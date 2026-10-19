#include "kg/RootWriter.h"

#include <set>
#include <stdexcept>
#include <vector>

#include "TFile.h"
#include "TTree.h"

namespace kg {

struct RootWriter::Impl {
  RootWriterOptions opt;
  TFile *file = nullptr;
  TTree *events = nullptr, *tasks = nullptr, *envs = nullptr;
  TTree *legacy = nullptr, *info = nullptr;
  std::set<std::string> seenTasks, seenEnvs;

  // events
  std::string task_id, env_id, process;
  Long64_t task_index = 0;
  ULong64_t event_index = 0, truth_seed = 0, detector_seed = 0;
  int component = 0, n_obj = 0, n_feat = 0, n_lost = 0, n_fakes = 0, n_misid = 0;
  std::vector<float> obj_px, obj_py, obj_pz, obj_e, obj_pt, obj_eta, obj_phi, obj_m, obj_charge, obj_feat;
  std::vector<int> obj_type, obj_truth, obj_fake;
  float met = 0, met_phi = 0, met_x = 0, met_y = 0, sum_et = 0, true_met_x = 0, true_met_y = 0, sqrt_s_hat = 0;
  std::vector<int> truth_pdg, truth_parent, truth_status, truth_species, truth_branch, truth_depth;
  std::vector<float> truth_px, truth_py, truth_pz, truth_e, truth_mass;

  // tasks
  std::string t_id, t_topology, t_yaml, t_prior;
  ULong64_t t_seed = 0, t_master = 0;
  Long64_t t_index = 0;
  std::vector<std::string> t_features, t_types, t_param_names, t_processes;
  std::vector<double> t_param_values, t_fractions;

  // environments
  std::string e_id, e_yaml;
  ULong64_t e_seed = 0;
  std::vector<std::string> e_param_names;
  std::vector<double> e_param_values;

  // legacy trees (test.cc layout)
  float T1Pt, T1Eta, T1Phi, T1M, T2Pt, T2Eta, T2Phi, T2M, METPt, METEta, METPhi;
  std::vector<float> P1Pt, P1Eta, P1Phi, P2Pt, P2Eta, P2Phi, truthM, truthMMin, truthMMax;
  std::vector<int> P1Id, P2Id, truthId;
};

RootWriter::RootWriter(const std::string &path, RootWriterOptions opt) : impl_(new Impl) {
  Impl &d = *impl_;
  d.opt = opt;
  d.file = TFile::Open(path.c_str(), "RECREATE");
  if (!d.file || d.file->IsZombie()) throw std::runtime_error("cannot create ROOT file '" + path + "'");
  if (opt.legacy_format) {
    d.legacy = new TTree("test", "test");
    d.info = new TTree("info", "info");
    d.legacy->Branch("P1Pt", &d.P1Pt);
    d.legacy->Branch("P1Eta", &d.P1Eta);
    d.legacy->Branch("P1Phi", &d.P1Phi);
    d.legacy->Branch("P1Id", &d.P1Id);
    d.legacy->Branch("P2Pt", &d.P2Pt);
    d.legacy->Branch("P2Eta", &d.P2Eta);
    d.legacy->Branch("P2Phi", &d.P2Phi);
    d.legacy->Branch("P2Id", &d.P2Id);
    d.legacy->Branch("T1Pt", &d.T1Pt, "T1Pt/F");
    d.legacy->Branch("T1Eta", &d.T1Eta, "T1Eta/F");
    d.legacy->Branch("T1Phi", &d.T1Phi, "T1Phi/F");
    d.legacy->Branch("T1M", &d.T1M, "T1M/F");
    d.legacy->Branch("T2Pt", &d.T2Pt, "T2Pt/F");
    d.legacy->Branch("T2Eta", &d.T2Eta, "T2Eta/F");
    d.legacy->Branch("T2Phi", &d.T2Phi, "T2Phi/F");
    d.legacy->Branch("T2M", &d.T2M, "T2M/F");
    d.legacy->Branch("METPt", &d.METPt, "METPt/F");
    d.legacy->Branch("METEta", &d.METEta, "METEta/F");
    d.legacy->Branch("METPhi", &d.METPhi, "METPhi/F");
    d.info->Branch("truthM", &d.truthM);
    d.info->Branch("truthMMin", &d.truthMMin);
    d.info->Branch("truthMMax", &d.truthMMax);
    d.info->Branch("truthId", &d.truthId);
    return;
  }
  d.events = new TTree("events", "kinegen observed events");
  d.events->Branch("task_id", &d.task_id);
  d.events->Branch("task_index", &d.task_index);
  d.events->Branch("env_id", &d.env_id);
  d.events->Branch("event_index", &d.event_index);
  d.events->Branch("truth_seed", &d.truth_seed);
  d.events->Branch("detector_seed", &d.detector_seed);
  d.events->Branch("component", &d.component);
  d.events->Branch("process", &d.process);
  d.events->Branch("n_obj", &d.n_obj);
  d.events->Branch("n_feat", &d.n_feat);
  d.events->Branch("obj_px", &d.obj_px);
  d.events->Branch("obj_py", &d.obj_py);
  d.events->Branch("obj_pz", &d.obj_pz);
  d.events->Branch("obj_e", &d.obj_e);
  d.events->Branch("obj_pt", &d.obj_pt);
  d.events->Branch("obj_eta", &d.obj_eta);
  d.events->Branch("obj_phi", &d.obj_phi);
  d.events->Branch("obj_m", &d.obj_m);
  d.events->Branch("obj_charge", &d.obj_charge);
  d.events->Branch("obj_type", &d.obj_type);
  d.events->Branch("obj_truth", &d.obj_truth);
  d.events->Branch("obj_fake", &d.obj_fake);
  d.events->Branch("obj_feat", &d.obj_feat);
  d.events->Branch("met", &d.met);
  d.events->Branch("met_phi", &d.met_phi);
  d.events->Branch("met_x", &d.met_x);
  d.events->Branch("met_y", &d.met_y);
  d.events->Branch("sum_et", &d.sum_et);
  d.events->Branch("n_lost", &d.n_lost);
  d.events->Branch("n_fakes", &d.n_fakes);
  d.events->Branch("n_misid", &d.n_misid);
  d.events->Branch("true_met_x", &d.true_met_x);
  d.events->Branch("true_met_y", &d.true_met_y);
  d.events->Branch("sqrt_s_hat", &d.sqrt_s_hat);
  if (opt.write_truth) {
    d.events->Branch("truth_pdg", &d.truth_pdg);
    d.events->Branch("truth_parent", &d.truth_parent);
    d.events->Branch("truth_status", &d.truth_status);
    d.events->Branch("truth_species", &d.truth_species);
    d.events->Branch("truth_branch", &d.truth_branch);
    d.events->Branch("truth_depth", &d.truth_depth);
    d.events->Branch("truth_px", &d.truth_px);
    d.events->Branch("truth_py", &d.truth_py);
    d.events->Branch("truth_pz", &d.truth_pz);
    d.events->Branch("truth_e", &d.truth_e);
    d.events->Branch("truth_mass", &d.truth_mass);
  }
  d.tasks = new TTree("tasks", "kinegen tasks");
  d.tasks->Branch("task_id", &d.t_id);
  d.tasks->Branch("task_index", &d.t_index);
  d.tasks->Branch("prior_name", &d.t_prior);
  d.tasks->Branch("seed", &d.t_seed);
  d.tasks->Branch("master_seed", &d.t_master);
  d.tasks->Branch("topology", &d.t_topology);
  d.tasks->Branch("processes", &d.t_processes);
  d.tasks->Branch("fractions", &d.t_fractions);
  d.tasks->Branch("feature_names", &d.t_features);
  d.tasks->Branch("object_types", &d.t_types);
  d.tasks->Branch("param_names", &d.t_param_names);
  d.tasks->Branch("param_values", &d.t_param_values);
  d.tasks->Branch("yaml", &d.t_yaml);
  d.envs = new TTree("environments", "kinegen detector environments");
  d.envs->Branch("env_id", &d.e_id);
  d.envs->Branch("seed", &d.e_seed);
  d.envs->Branch("param_names", &d.e_param_names);
  d.envs->Branch("param_values", &d.e_param_values);
  d.envs->Branch("yaml", &d.e_yaml);
}

RootWriter::~RootWriter() {
  try {
    close();
  } catch (...) {
  }
}

void RootWriter::writeTask(const TaskConfig &t) {
  Impl &d = *impl_;
  if (!d.seenTasks.insert(t.task_id).second) return;
  if (d.opt.legacy_format) {
    if (d.info->GetEntries() > 0) return;
    // test.cc wrote one entry with the minimum and maximum mass of every particle.
    for (const auto &p : t.particles) {
      const ParamDef &m = t.params.def(p.mass);
      double lo = t.values[p.mass], hi = lo;
      if (m.kind == ParamDef::Kind::Dist && m.dist.type == Distribution::Type::Uniform) {
        lo = m.dist.min;
        hi = m.dist.max;
      }
      d.truthM.push_back(lo);
      d.truthMMin.push_back(lo);
      d.truthMMax.push_back(hi);
      d.truthId.push_back(p.pdg);
    }
    d.info->Fill();
    return;
  }
  d.t_id = t.task_id;
  d.t_index = t.task_index;
  d.t_prior = t.prior_name;
  d.t_seed = t.seed;
  d.t_master = t.master_seed;
  d.t_topology = t.topology();
  d.t_processes.clear();
  d.t_fractions.clear();
  for (const auto &c : t.components) {
    d.t_processes.push_back(c.process);
    d.t_fractions.push_back(c.fraction);
  }
  d.t_features = *t.feature_names;
  d.t_types = t.object_types;
  d.t_param_names.clear();
  d.t_param_values.clear();
  for (const auto &kv : t.parameterValues()) {
    d.t_param_names.push_back(kv.first);
    d.t_param_values.push_back(kv.second);
  }
  d.t_yaml = taskToYAML(t);
  d.tasks->Fill();
  writeEnvironment(t.environment);
}

void RootWriter::writeEnvironment(const Environment &env) {
  Impl &d = *impl_;
  if (d.opt.legacy_format || !d.seenEnvs.insert(env.env_id).second) return;
  d.e_id = env.env_id;
  d.e_seed = env.seed;
  d.e_param_names.clear();
  d.e_param_values.clear();
  for (const auto &kv : env.parameterValues()) {
    d.e_param_names.push_back(kv.first);
    d.e_param_values.push_back(kv.second);
  }
  d.e_yaml = environmentToYAML(env);
  d.envs->Fill();
}

void RootWriter::fill(const TaskConfig &t, const TruthEvent &truth, const ObservedEvent &reco) {
  Impl &d = *impl_;
  writeTask(t);
  if (d.opt.legacy_format) {
    d.P1Pt.clear(); d.P1Eta.clear(); d.P1Phi.clear(); d.P1Id.clear();
    d.P2Pt.clear(); d.P2Eta.clear(); d.P2Phi.clear(); d.P2Id.clear();
    P4 met;
    for (const auto &p : truth.particles) {
      if (p.species < 0) continue;
      if (p.parent < 0 && p.branch <= 1) {
        (p.branch == 0 ? d.T1Pt : d.T2Pt) = p.p4.pt();
        (p.branch == 0 ? d.T1Eta : d.T2Eta) = p.p4.eta();
        (p.branch == 0 ? d.T1Phi : d.T2Phi) = p.p4.phi();
        (p.branch == 0 ? d.T1M : d.T2M) = p.mass;
      }
      if (!p.isFinal()) continue;
      if (p.object_type == "invisible") met += p.p4;
      if (p.branch == 0) {
        d.P1Pt.push_back(p.p4.pt()); d.P1Eta.push_back(p.p4.eta()); d.P1Phi.push_back(p.p4.phi()); d.P1Id.push_back(p.pdg);
      } else if (p.branch == 1) {
        d.P2Pt.push_back(p.p4.pt()); d.P2Eta.push_back(p.p4.eta()); d.P2Phi.push_back(p.p4.phi()); d.P2Id.push_back(p.pdg);
      }
    }
    d.METPt = met.pt();
    d.METEta = met.pt() > 0 ? met.eta() : 0;
    d.METPhi = met.phi();
    d.legacy->Fill();
    return;
  }
  d.task_id = truth.task_id;
  d.task_index = t.task_index;
  d.env_id = reco.env_id;
  d.event_index = reco.event_index;
  d.truth_seed = reco.truth_seed;
  d.detector_seed = reco.detector_seed;
  d.component = reco.component;
  d.process = reco.process;
  d.n_obj = static_cast<int>(reco.objects.size());
  d.n_feat = reco.feature_names ? static_cast<int>(reco.feature_names->size()) : 0;
  for (auto *v : {&d.obj_px, &d.obj_py, &d.obj_pz, &d.obj_e, &d.obj_pt, &d.obj_eta, &d.obj_phi, &d.obj_m, &d.obj_charge, &d.obj_feat}) v->clear();
  d.obj_type.clear();
  d.obj_truth.clear();
  d.obj_fake.clear();
  for (const auto &o : reco.objects) {
    d.obj_px.push_back(o.p4.px);
    d.obj_py.push_back(o.p4.py);
    d.obj_pz.push_back(o.p4.pz);
    d.obj_e.push_back(o.p4.e);
    d.obj_pt.push_back(o.pt());
    d.obj_eta.push_back(o.eta());
    d.obj_phi.push_back(o.phi());
    d.obj_m.push_back(o.m());
    d.obj_charge.push_back(o.charge);
    d.obj_type.push_back(o.type_code);
    d.obj_truth.push_back(o.truth_index);
    d.obj_fake.push_back(o.is_fake);
    for (double f : o.features) d.obj_feat.push_back(static_cast<float>(f));
  }
  d.met = reco.met();
  d.met_phi = reco.metPhi();
  d.met_x = reco.met_x;
  d.met_y = reco.met_y;
  d.sum_et = reco.sum_et;
  d.n_lost = reco.n_lost;
  d.n_fakes = reco.n_fakes;
  d.n_misid = reco.n_misid;
  const auto tm = truth.trueMET();
  d.true_met_x = tm.first;
  d.true_met_y = tm.second;
  d.sqrt_s_hat = truth.sqrt_s_hat;
  if (d.opt.write_truth) {
    for (auto *v : {&d.truth_pdg, &d.truth_parent, &d.truth_status, &d.truth_species, &d.truth_branch, &d.truth_depth}) v->clear();
    for (auto *v : {&d.truth_px, &d.truth_py, &d.truth_pz, &d.truth_e, &d.truth_mass}) v->clear();
    for (const auto &p : truth.particles) {
      d.truth_pdg.push_back(p.pdg);
      d.truth_parent.push_back(p.parent);
      d.truth_status.push_back(static_cast<int>(p.status));
      d.truth_species.push_back(p.species);
      d.truth_branch.push_back(p.branch);
      d.truth_depth.push_back(p.depth);
      d.truth_px.push_back(p.p4.px);
      d.truth_py.push_back(p.p4.py);
      d.truth_pz.push_back(p.p4.pz);
      d.truth_e.push_back(p.p4.e);
      d.truth_mass.push_back(p.mass);
    }
  }
  d.events->Fill();
}

void RootWriter::close() {
  Impl &d = *impl_;
  if (!d.file) return;
  d.file->cd();
  d.file->Write();
  d.file->Close();
  delete d.file;
  d.file = nullptr;
}

} // namespace kg
