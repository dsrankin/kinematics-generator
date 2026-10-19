// Observation model: smearing, efficiency, acceptance, fakes, misidentification, MET.
#include <map>
#include <set>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

namespace {

// pp -> X -> mu+ mu- + invisible, with a configurable detector block.
PriorConfig withDetector(const std::string &det, const std::string &extra = "") {
  return kgtest::prior(R"(
particles:
  X: {pdg: 32, mass: 400, decays: [[mu+, mu-, N]]}
  mu+: {pdg: -13, charge: 1}
  mu-: {pdg: 13, charge: -1}
  N: {pdg: 1000022, mass: 50, visible: false}
production: [X]
initial_state:
  sqrt_s: {model: exponential_pdf}
  recoil: {model: exponential, min_pt: 10, mean_pt: 40}
)" + extra + "detector:\n" + det);
}

} // namespace

TEST(perfect_detector_reproduces_truth_and_met) {
  const PriorConfig p = withDetector("  objects: {}\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  for (int i = 0; i < 1000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev);
    CHECK(ob.objects.size() == 2);
    for (const auto &o : ob.objects) {
      const TruthParticle &tp = ev.particles.at(o.truth_index);
      CHECK(o.type == "muon" && o.true_type == "muon" && !o.is_fake);
      CHECK_CLOSE(o.p4.pt(), tp.p4.pt(), 1e-9 * tp.p4.e);
      CHECK_CLOSE(o.p4.eta(), tp.p4.eta(), 1e-9);
      CHECK_CLOSE(o.charge, tp.charge, 0);
    }
    CHECK(ob.objects[0].pt() >= ob.objects[1].pt());
    const auto met = ev.trueMET();
    CHECK_CLOSE(ob.met_x, met.first, 1e-8 * ev.sqrt_s_hat);
    CHECK_CLOSE(ob.met_y, met.second, 1e-8 * ev.sqrt_s_hat);
  }
}

TEST(momentum_and_angular_smearing) {
  const PriorConfig p = withDetector("  objects: {muon: {resolution: {constant: 0.1}, sigma_eta: 0.05, sigma_phi: 0.02, scale: 1.02}}\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  Stats r, de, dp;
  for (int i = 0; i < 20000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    for (const auto &o : det.apply(ev).objects) {
      const P4 &tp = ev.particles[o.truth_index].p4;
      r.add(o.pt() / tp.pt());
      de.add(o.eta() - tp.eta());
      dp.add(deltaPhi(o.phi(), tp.phi()));
    }
  }
  CHECK_CLOSE(r.mean(), 1.02, 4 * r.err());
  CHECK_CLOSE(r.sd(), 0.102, 0.002);
  CHECK_CLOSE(de.sd(), 0.05, 0.001);
  CHECK_CLOSE(dp.sd(), 0.02, 0.0005);
}

TEST(calorimeter_resolution_formula) {
  Response r;
  r.stochastic = 0.5;
  r.constant = 0.03;
  r.noise = 2;
  CHECK_CLOSE(r.relativeResolution(100, 50), std::sqrt(0.25 / 100 + 0.0009 + 4.0 / 10000), 1e-12);
  r.efficiency = 0.9;
  r.eff_pt50 = 30;
  r.eff_width = 5;
  CHECK_CLOSE(r.efficiencyAt(30), 0.45, 1e-12);
  CHECK(r.efficiencyAt(100) > 0.8999);
  CHECK(r.efficiencyAt(5) < 1e-5);
}

TEST(efficiency_acceptance_and_thresholds) {
  const PriorConfig p = withDetector("  objects: {muon: {efficiency: 0.7, eta_max: 2.0, pt_min: 30}}\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  int nAcc = 0, nObs = 0;
  for (int i = 0; i < 20000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    for (int f : ev.finalState())
      if (ev.particles[f].visible() && std::fabs(ev.particles[f].p4.eta()) <= 2.0 && ev.particles[f].p4.pt() >= 30) ++nAcc;
    const ObservedEvent ob = det.apply(ev);
    for (const auto &o : ob.objects) {
      CHECK(std::fabs(o.eta()) <= 2.0 && o.pt() >= 30);
      ++nObs;
    }
    CHECK(ob.n_lost + int(ob.objects.size()) == 2);
  }
  CHECK_CLOSE(double(nObs) / nAcc, 0.7, 4 * std::sqrt(0.21 / nAcc));
}

TEST(fakes_and_misidentification) {
  const PriorConfig p = withDetector(R"(
  objects:
    muon: {misid: {electron: 0.2}}
  fakes: {rate: 1.5, types: {jet: 3, photon: 1}, pt_min: 15, pt_scale: 10, eta_max: 2.0}
)");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  Stats nf;
  int nMis = 0, nReal = 0, nJet = 0, nFake = 0;
  for (int i = 0; i < 20000; ++i) {
    const ObservedEvent ob = det.apply(generateTruthEvent(t, i));
    int fakes = 0;
    for (const auto &o : ob.objects) {
      if (o.is_fake) {
        ++fakes;
        ++nFake;
        CHECK(o.truth_index == -1 && o.true_type == "fake");
        CHECK(o.pt() >= 15 && std::fabs(o.eta()) <= 2.0);
        CHECK(o.type == "jet" || o.type == "photon");
        nJet += o.type == "jet";
        CHECK(o.type_code == t.typeCode(o.type) && o.type_code >= 0);
      } else {
        ++nReal;
        CHECK(o.true_type == "muon");
        if (o.type == "electron") ++nMis;
      }
    }
    CHECK(fakes == ob.n_fakes);
    nf.add(fakes);
  }
  CHECK_CLOSE(nf.mean(), 1.5, 4 * nf.err());
  CHECK_CLOSE(nf.sd() * nf.sd(), 1.5, 0.06); // Poisson
  CHECK_CLOSE(double(nMis) / nReal, 0.2, 4 * std::sqrt(0.16 / nReal));
  CHECK_CLOSE(double(nJet) / nFake, 0.75, 4 * std::sqrt(0.1875 / nFake));
}

TEST(met_reconstruction_and_smearing) {
  const PriorConfig p = withDetector("  objects: {}\n  met: {soft_sigma: 10}\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  Stats dx;
  for (int i = 0; i < 20000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev);
    dx.add(ob.met_x - ev.trueMET().first);
  }
  CHECK_CLOSE(dx.mean(), 0, 4 * dx.err());
  CHECK_CLOSE(dx.sd(), 10, 0.2);
  // Lost visible objects become missing momentum when lost_to_met = 1.
  const PriorConfig q = withDetector("  objects: {muon: {efficiency: 0}}\n");
  const TaskConfig u = sampleTask(q, 1, 0);
  for (int i = 0; i < 100; ++i) {
    const TruthEvent ev = generateTruthEvent(u, i);
    const ObservedEvent ob = DetectorModel(u).apply(ev);
    CHECK(ob.objects.empty() && ob.n_lost == 2);
    double vx = 0;
    for (int f : ev.finalState()) vx += ev.particles[f].p4.px; // visible + invisible
    CHECK_CLOSE(ob.met_x, vx, 1e-8 * ev.sqrt_s_hat);
  }
}

TEST(recoil_as_object_or_unclustered) {
  const PriorConfig p = kgtest::prior(R"(
particles:
  X: {pdg: 32, mass: 400, decays: [[mu+, mu-]]}
  mu+: {pdg: -13}
  mu-: {pdg: 13}
production: [X]
initial_state: {recoil: {model: exponential, min_pt: 50, mean_pt: 60, emit: true, object: jet}}
)");
  const TaskConfig t = sampleTask(p, 1, 0);
  for (int i = 0; i < 200; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = DetectorModel(t).apply(ev);
    CHECK(ob.objects.size() == 3);
    int njet = 0;
    for (const auto &o : ob.objects) njet += o.type == "jet";
    CHECK(njet == 1);
    CHECK_CLOSE(ob.met(), 0, 1e-6 * ev.sqrt_s_hat);
  }
}

TEST(environment_randomization_and_event_scope) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("detector_randomization.yaml"));
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  std::set<double> stoch;
  std::vector<Environment> envs;
  for (int k = 0; k < 10; ++k) {
    envs.push_back(sampleEnvironment(p.detector, 5, int64_t(k)));
    stoch.insert(envs.back().value("detector.jet.resolution.stochastic"));
    CHECK(envs.back().hasEventParams());
  }
  CHECK(stoch.size() == 10);
  const TruthEvent ev = generateTruthEvent(t, 3);
  std::set<double> mets;
  for (const auto &e : envs) {
    const ObservedEvent ob = det.apply(ev, e);
    CHECK(ob.env_id == e.env_id && ob.truth_seed == ev.seed);
    CHECK(ob.detector_params.size() == 1 && ob.detector_params[0].first == "detector.met.soft_sigma");
    mets.insert(ob.met_x);
  }
  CHECK(mets.size() == envs.size());
  // Event-scoped soft term changes between events of one environment.
  std::set<double> soft;
  for (int i = 0; i < 20; ++i) soft.insert(det.apply(generateTruthEvent(t, i), envs[0]).detector_params[0].second);
  CHECK(soft.size() == 20);
}
