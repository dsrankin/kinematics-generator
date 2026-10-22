// ROOT output and statistical comparison with the original myGenerator.h kinematics.
#include <cstdio>

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

#include "../myGenerator.h"
#include "helpers.h"
#include "kg/RootWriter.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

TEST(root_writer_new_format) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("randomized_tasks.yaml"));
  long long nObj = 0;
  {
    RootWriter w("kg_test_out.root");
    for (int it = 0; it < 3; ++it) {
      const TaskConfig t = sampleTask(p, 1, it);
      const DetectorModel det(t);
      for (int i = 0; i < 100; ++i) {
        const TruthEvent ev = generateTruthEvent(t, i);
        const ObservedEvent ob = det.apply(ev);
        nObj += ob.objects.size();
        w.fill(t, ev, ob);
      }
    }
  }
  TFile f("kg_test_out.root");
  TTree *ev = static_cast<TTree *>(f.Get("events"));
  TTree *tasks = static_cast<TTree *>(f.Get("tasks"));
  TTree *envs = static_cast<TTree *>(f.Get("environments"));
  CHECK(ev && tasks && envs);
  CHECK(ev->GetEntries() == 300 && tasks->GetEntries() == 3 && envs->GetEntries() == 3);
  std::vector<float> *pt = nullptr;
  std::vector<int> *pdg = nullptr;
  std::string *tid = nullptr;
  ev->SetBranchAddress("obj_pt", &pt);
  ev->SetBranchAddress("truth_pdg", &pdg);
  ev->SetBranchAddress("task_id", &tid);
  long long n = 0;
  for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
    ev->GetEntry(i);
    n += pt->size();
    CHECK(!pdg->empty());
  }
  CHECK(n == nObj);
  // The stored task YAML regenerates the stored events.
  std::string *yaml = nullptr;
  tasks->SetBranchAddress("yaml", &yaml);
  tasks->GetEntry(0);
  const TaskConfig back = taskFromYAML(*yaml);
  ev->GetEntry(0);
  CHECK(*tid == back.task_id);
  CHECK(pdg->size() == generateTruthEvent(back, 0).particles.size());
  std::remove("kg_test_out.root");
}

TEST(root_writer_legacy_format) {
  const PriorConfig p = priorFromLegacyCard(std::string(KG_SOURCE_DIR) + "/process.card");
  const TaskConfig t = sampleTask(p, 0, 0);
  {
    RootWriterOptions o;
    o.legacy_format = true;
    RootWriter w("kg_test_legacy.root", o);
    for (int i = 0; i < 200; ++i) {
      const TruthEvent ev = generateTruthEvent(t, i);
      w.fill(t, ev, DetectorModel(t).apply(ev));
    }
  }
  TFile f("kg_test_legacy.root");
  TTree *test = static_cast<TTree *>(f.Get("test"));
  TTree *info = static_cast<TTree *>(f.Get("info"));
  CHECK(test && info && test->GetEntries() == 200 && info->GetEntries() == 1);
  for (const char *b : {"P1Pt", "P1Eta", "P1Phi", "P1Id", "P2Pt", "P2Eta", "P2Phi", "P2Id", "T1Pt", "T1Eta", "T1Phi", "T1M",
                        "T2Pt", "T2Eta", "T2Phi", "T2M", "METPt", "METEta", "METPhi"})
    CHECK_MSG(test->GetBranch(b), b);
  std::vector<int> *p2id = nullptr;
  float t2m = 0, met = 0;
  test->SetBranchAddress("P2Id", &p2id);
  test->SetBranchAddress("T2M", &t2m);
  test->SetBranchAddress("METPt", &met);
  test->GetEntry(0);
  CHECK(p2id->size() == 3); // b, l, v
  CHECK(t2m >= 165 && t2m < 180);
  CHECK(met > 0);
  std::remove("kg_test_legacy.root");
}

TEST(matches_original_myGenerator_kinematics) {
  // Same physics through the old helpers and through the converted card:
  // compare sqrt(s), recoil pT and the decay-product pT spectrum.
  const PriorConfig p = kgtest::prior(convertLegacyCardString(R"(
production A B
particle A 400 400 1 0 1 q q
particle B 300 300 2 0 1 q l
particle q 0 0 1 1 1 x x
particle l 0 0 11 1 1 x x
)")
                                          .yaml);
  const TaskConfig t = sampleTask(p, 0, 0);
  Stats sNew, ptNew, pzNew;
  const int n = 40000;
  for (int i = 0; i < n; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    sNew.add(ev.sqrt_s_hat);
    P4 sys;
    for (int k : ev.producedParticles()) sys += ev.particles[k].p4;
    ptNew.add(sys.pt());
    for (int f : ev.finalState())
      if (ev.particles[f].name == "l") pzNew.add(std::fabs(ev.particles[f].p4.pt()));
  }
  delete gRandom;
  gRandom = new TRandom3(1);
  Stats sOld, ptOld, pzOld;
  for (int i = 0; i < n; ++i) {
    const tlv cm = myGenerator::getRandCM_softRad(400 + 300 + 1);
    sOld.add(cm.M());
    ptOld.add(cm.Pt());
    tlv a, b, q1, q2, q3, l;
    myGenerator::getRand2BodyDecay(a, b, cm, 400, 300);
    myGenerator::getRand2BodyDecay(q3, l, b, 0, 0);
    pzOld.add(l.Pt());
  }
  CHECK_CLOSE(sNew.mean(), sOld.mean(), 4 * std::hypot(sNew.err(), sOld.err()));
  CHECK_CLOSE(sNew.sd(), sOld.sd(), 0.03 * sOld.sd());
  CHECK_CLOSE(ptNew.mean(), ptOld.mean(), 4 * std::hypot(ptNew.err(), ptOld.err()));
  CHECK_CLOSE(pzNew.mean(), pzOld.mean(), 4 * std::hypot(pzNew.err(), pzOld.err()));
}
