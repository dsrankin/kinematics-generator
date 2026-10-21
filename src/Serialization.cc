#include "kg/Serialization.h"

#include <fstream>
#include <sstream>

#include <yaml-cpp/yaml.h>

#include "kg/Task.h"

namespace kg {

namespace {

const int kTaskFormat = 1;

void emitAngular(YAML::Emitter &out, const AngularModel &a) {
  out << YAML::Flow << YAML::BeginMap;
  out << YAML::Key << "type" << YAML::Value << (a.type == AngularModel::Type::Isotropic ? "isotropic" : "polynomial");
  out << YAML::Key << "axis" << YAML::Value << axisName(a.axis);
  if (a.type == AngularModel::Type::Polynomial) out << YAML::Key << "coeffs" << YAML::Value << YAML::Flow << a.coeffs;
  out << YAML::EndMap;
}

AngularModel parseAngular(const YAML::Node &n) {
  AngularModel a;
  a.type = n["type"].as<std::string>() == "polynomial" ? AngularModel::Type::Polynomial : AngularModel::Type::Isotropic;
  a.axis = axisFromString(n["axis"].as<std::string>());
  if (n["coeffs"]) a.coeffs = n["coeffs"].as<std::vector<double>>();
  a.prepare();
  return a;
}

void emitP4(YAML::Emitter &out, const P4 &p) { out << YAML::Flow << std::vector<double>{p.px, p.py, p.pz, p.e}; }

P4 parseP4(const YAML::Node &n) {
  const auto v = n.as<std::vector<double>>();
  if (v.size() != 4) throw std::invalid_argument("four-vector must have 4 entries");
  return {v[0], v[1], v[2], v[3]};
}

const char *statusName(TruthStatus s) {
  switch (s) {
    case TruthStatus::Incoming: return "incoming";
    case TruthStatus::Recoil: return "recoil";
    case TruthStatus::Decayed: return "decayed";
    case TruthStatus::Final: return "final";
  }
  return "final";
}

TruthStatus statusFromString(const std::string &s) {
  if (s == "incoming") return TruthStatus::Incoming;
  if (s == "recoil") return TruthStatus::Recoil;
  if (s == "decayed") return TruthStatus::Decayed;
  if (s == "final") return TruthStatus::Final;
  throw std::invalid_argument("unknown truth status '" + s + "'");
}

} // namespace

// ---------------------------------------------------------------------------
// Tasks

std::string taskToYAML(const TaskConfig &t) {
  YAML::Emitter out;
  out.SetDoublePrecision(17);
  out << YAML::BeginMap;
  out << YAML::Key << "kinegen_task" << YAML::Value << kTaskFormat;
  out << YAML::Key << "task_id" << YAML::Value << t.task_id;
  out << YAML::Key << "prior_name" << YAML::Value << t.prior_name;
  out << YAML::Key << "seed" << YAML::Value << t.seed;
  out << YAML::Key << "master_seed" << YAML::Value << t.master_seed;
  out << YAML::Key << "task_index" << YAML::Value << t.task_index;
  out << YAML::Key << "attempts" << YAML::Value << t.attempts;
  out << YAML::Key << "topology" << YAML::Value << t.topology();
  out << YAML::Key << "bw_cutoff" << YAML::Value << t.bw_cutoff;
  out << YAML::Key << "parameters" << YAML::Value;
  t.params.emit(out);
  out << YAML::Key << "values" << YAML::Value << YAML::Flow << t.values;
  out << YAML::Key << "particles" << YAML::Value << YAML::BeginSeq;
  for (const auto &p : t.particles) {
    out << YAML::BeginMap;
    out << YAML::Key << "name" << YAML::Value << p.name;
    out << YAML::Key << "pdg" << YAML::Value << p.pdg;
    out << YAML::Key << "charge" << YAML::Value << p.charge;
    out << YAML::Key << "object" << YAML::Value << p.object;
    if (p.fixed_mode >= 0) out << YAML::Key << "fixed_mode" << YAML::Value << p.fixed_mode;
    if (!p.modes.empty()) {
      out << YAML::Key << "decays" << YAML::Value << YAML::BeginSeq;
      for (const auto &m : p.modes) {
        std::vector<std::string> d;
        for (int s : m.daughters) d.push_back(t.particles[s].name);
        out << YAML::Flow << YAML::BeginMap;
        out << YAML::Key << "daughters" << YAML::Value << YAML::Flow << d;
        out << YAML::Key << "br" << YAML::Value << m.br;
        out << YAML::Key << "model" << YAML::Value << threeBodyModelName(m.model);
        if (!m.reconstruct_as.empty()) out << YAML::Key << "reconstruct_as" << YAML::Value << m.reconstruct_as;
        out << YAML::Key << "angular" << YAML::Value;
        emitAngular(out, m.angular);
        out << YAML::EndMap;
      }
      out << YAML::EndSeq;
    }
    out << YAML::EndMap;
  }
  out << YAML::EndSeq;
  out << YAML::Key << "components" << YAML::Value << YAML::BeginSeq;
  for (const auto &c : t.components) {
    std::vector<std::string> d;
    for (int s : c.produced) d.push_back(t.particles[s].name);
    out << YAML::BeginMap;
    out << YAML::Key << "process" << YAML::Value << c.process;
    out << YAML::Key << "produced" << YAML::Value << YAML::Flow << d;
    out << YAML::Key << "fraction" << YAML::Value << c.fraction;
    out << YAML::Key << "angular" << YAML::Value;
    emitAngular(out, c.angular);
    out << YAML::EndMap;
  }
  out << YAML::EndSeq;
  out << YAML::Key << "initial_state" << YAML::Value;
  emitInitialStateFlags(out, t.initial);
  out << YAML::Key << "object_types" << YAML::Value << YAML::Flow << t.object_types;
  out << YAML::Key << "features" << YAML::Value << YAML::BeginSeq;
  for (const auto &f : t.features) emitFeatureSpec(out, f);
  out << YAML::EndSeq;
  out << YAML::Key << "environment" << YAML::Value;
  emitEnvironment(out, t.environment);
  out << YAML::EndMap;
  return std::string(out.c_str()) + "\n";
}

TaskConfig taskFromYAML(const std::string &yaml) {
  const YAML::Node n = YAML::Load(yaml);
  if (!n["kinegen_task"]) throw std::invalid_argument("not a kinegen task (missing 'kinegen_task')");
  if (n["kinegen_task"].as<int>() != kTaskFormat) throw std::invalid_argument("unsupported task format version");
  TaskConfig t;
  t.task_id = n["task_id"].as<std::string>();
  t.prior_name = n["prior_name"].as<std::string>();
  t.seed = n["seed"].as<uint64_t>();
  t.master_seed = n["master_seed"].as<uint64_t>();
  t.task_index = n["task_index"].as<int64_t>();
  t.attempts = n["attempts"].as<int>();
  t.bw_cutoff = n["bw_cutoff"].as<double>();
  t.params = ParamSet::fromYAML(n["parameters"]);
  t.values = n["values"].as<std::vector<double>>();
  if (t.values.size() != t.params.size()) throw std::invalid_argument("task: 'values' and 'parameters' differ in length");
  const YAML::Node parts = n["particles"];
  for (size_t i = 0; i < parts.size(); ++i) t.particle_index[parts[i]["name"].as<std::string>()] = static_cast<int>(i);
  for (const auto &pn : parts) {
    TaskParticle p;
    p.name = pn["name"].as<std::string>();
    p.pdg = pn["pdg"].as<int>();
    p.charge = pn["charge"].as<double>();
    p.object = pn["object"].as<std::string>();
    p.mass = t.params.at(p.name + ".mass");
    p.width = t.params.at(p.name + ".width");
    if (pn["fixed_mode"]) p.fixed_mode = pn["fixed_mode"].as<int>();
    if (pn["decays"])
      for (const auto &mn : pn["decays"]) {
        TaskDecayMode m;
        for (const auto &d : mn["daughters"]) m.daughters.push_back(t.species(d.as<std::string>()));
        m.br = mn["br"].as<double>();
        m.model = threeBodyModelFromString(mn["model"].as<std::string>());
        m.angular = parseAngular(mn["angular"]);
        if (mn["reconstruct_as"]) m.reconstruct_as = mn["reconstruct_as"].as<std::string>();
        p.br_weights.push_back(m.br);
        p.modes.push_back(std::move(m));
      }
    p.stable = p.modes.empty();
    t.particles.push_back(std::move(p));
  }
  for (const auto &cn : n["components"]) {
    TaskComponent c;
    c.process = cn["process"].as<std::string>();
    for (const auto &d : cn["produced"]) c.produced.push_back(t.species(d.as<std::string>()));
    c.fraction = cn["fraction"].as<double>();
    c.angular = parseAngular(cn["angular"]);
    t.components.push_back(std::move(c));
  }
  const YAML::Node is = n["initial_state"];
  const std::string model = is["sqrt_s_model"].as<std::string>();
  for (SqrtSModel m : {SqrtSModel::ExponentialPdf, SqrtSModel::ThresholdExponential, SqrtSModel::Fixed, SqrtSModel::Uniform})
    if (sqrtSModelName(m) == model) t.initial.model = m;
  t.initial.rapidity_from_pdf = is["rapidity_from_pdf"].as<bool>();
  const std::string rec = is["recoil_model"].as<std::string>();
  t.initial.recoil = rec == "none" ? RecoilModel::None : (rec == "exponential" ? RecoilModel::Exponential : RecoilModel::Param);
  t.initial.emit_recoil = is["emit_recoil"].as<bool>();
  t.initial.recoil_object = is["recoil_object"].as<std::string>();
  t.initial.max_attempts = is["max_attempts"].as<int>();
  bindInitialStateParams(t.initial, t.params);
  t.object_types = n["object_types"].as<std::vector<std::string>>();
  auto names = std::make_shared<std::vector<std::string>>();
  for (const auto &fn : n["features"]) {
    t.features.push_back(parseFeatureSpec(fn));
    names->push_back(t.features.back().name);
  }
  t.feature_names = names;
  t.environment = parseEnvironment(n["environment"]);
  validateTask(t);
  return t;
}

void saveTask(const TaskConfig &task, const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write '" + path + "'");
  out << taskToYAML(task);
}

TaskConfig loadTask(const std::string &path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open task file '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return taskFromYAML(ss.str());
}

// ---------------------------------------------------------------------------
// Events

std::string truthEventToYAML(const TruthEvent &ev) {
  YAML::Emitter out;
  out.SetDoublePrecision(17);
  out << YAML::BeginMap;
  out << YAML::Key << "task_id" << YAML::Value << ev.task_id;
  out << YAML::Key << "task_seed" << YAML::Value << ev.task_seed;
  out << YAML::Key << "event_index" << YAML::Value << ev.event_index;
  out << YAML::Key << "seed" << YAML::Value << ev.seed;
  out << YAML::Key << "component" << YAML::Value << ev.component;
  out << YAML::Key << "process" << YAML::Value << ev.process;
  out << YAML::Key << "sqrt_s_hat" << YAML::Value << ev.sqrt_s_hat;
  out << YAML::Key << "x1" << YAML::Value << ev.x1;
  out << YAML::Key << "x2" << YAML::Value << ev.x2;
  out << YAML::Key << "event_params" << YAML::Value << YAML::Flow << YAML::BeginMap;
  for (const auto &kv : ev.event_params) out << YAML::Key << kv.first << YAML::Value << kv.second;
  out << YAML::EndMap;
  out << YAML::Key << "particles" << YAML::Value << YAML::BeginSeq;
  for (const auto &p : ev.particles) {
    out << YAML::Flow << YAML::BeginMap;
    out << YAML::Key << "i" << YAML::Value << p.index;
    out << YAML::Key << "name" << YAML::Value << p.name;
    out << YAML::Key << "pdg" << YAML::Value << p.pdg;
    out << YAML::Key << "charge" << YAML::Value << p.charge;
    out << YAML::Key << "p4" << YAML::Value;
    emitP4(out, p.p4);
    out << YAML::Key << "mass" << YAML::Value << p.mass;
    out << YAML::Key << "parent" << YAML::Value << p.parent;
    out << YAML::Key << "children" << YAML::Value << YAML::Flow << p.children;
    out << YAML::Key << "status" << YAML::Value << statusName(p.status);
    out << YAML::Key << "species" << YAML::Value << p.species;
    out << YAML::Key << "branch" << YAML::Value << p.branch;
    out << YAML::Key << "depth" << YAML::Value << p.depth;
    out << YAML::Key << "decay_mode" << YAML::Value << p.decay_mode;
    out << YAML::Key << "object" << YAML::Value << p.object_type;
    if (!p.reconstruct_as.empty()) out << YAML::Key << "reconstruct_as" << YAML::Value << p.reconstruct_as;
    out << YAML::EndMap;
  }
  out << YAML::EndSeq << YAML::EndMap;
  return std::string(out.c_str()) + "\n";
}

TruthEvent truthEventFromYAML(const std::string &yaml) {
  const YAML::Node n = YAML::Load(yaml);
  TruthEvent ev;
  ev.task_id = n["task_id"].as<std::string>();
  ev.task_seed = n["task_seed"].as<uint64_t>();
  ev.event_index = n["event_index"].as<uint64_t>();
  ev.seed = n["seed"].as<uint64_t>();
  ev.component = n["component"].as<int>();
  ev.process = n["process"].as<std::string>();
  ev.sqrt_s_hat = n["sqrt_s_hat"].as<double>();
  ev.x1 = n["x1"].as<double>();
  ev.x2 = n["x2"].as<double>();
  for (const auto &kv : n["event_params"]) ev.event_params.emplace_back(kv.first.as<std::string>(), kv.second.as<double>());
  for (const auto &pn : n["particles"]) {
    TruthParticle p;
    p.index = pn["i"].as<int>();
    p.name = pn["name"].as<std::string>();
    p.pdg = pn["pdg"].as<int>();
    p.charge = pn["charge"].as<double>();
    p.p4 = parseP4(pn["p4"]);
    p.mass = pn["mass"].as<double>();
    p.parent = pn["parent"].as<int>();
    p.children = pn["children"].as<std::vector<int>>();
    p.status = statusFromString(pn["status"].as<std::string>());
    p.species = pn["species"].as<int>();
    p.branch = pn["branch"].as<int>();
    p.depth = pn["depth"].as<int>();
    p.decay_mode = pn["decay_mode"].as<int>();
    p.object_type = pn["object"].as<std::string>();
    if (pn["reconstruct_as"]) p.reconstruct_as = pn["reconstruct_as"].as<std::string>();
    ev.particles.push_back(std::move(p));
  }
  return ev;
}

std::string observedEventToYAML(const ObservedEvent &ev) {
  YAML::Emitter out;
  out.SetDoublePrecision(17);
  out << YAML::BeginMap;
  out << YAML::Key << "task_id" << YAML::Value << ev.task_id;
  out << YAML::Key << "env_id" << YAML::Value << ev.env_id;
  out << YAML::Key << "event_index" << YAML::Value << ev.event_index;
  out << YAML::Key << "truth_seed" << YAML::Value << ev.truth_seed;
  out << YAML::Key << "detector_seed" << YAML::Value << ev.detector_seed;
  out << YAML::Key << "component" << YAML::Value << ev.component;
  out << YAML::Key << "process" << YAML::Value << ev.process;
  out << YAML::Key << "met" << YAML::Value << YAML::Flow << std::vector<double>{ev.met_x, ev.met_y};
  out << YAML::Key << "sum_et" << YAML::Value << ev.sum_et;
  out << YAML::Key << "n_lost" << YAML::Value << ev.n_lost;
  out << YAML::Key << "n_fakes" << YAML::Value << ev.n_fakes;
  out << YAML::Key << "n_misid" << YAML::Value << ev.n_misid;
  out << YAML::Key << "feature_names" << YAML::Value << YAML::Flow
      << (ev.feature_names ? *ev.feature_names : std::vector<std::string>{});
  out << YAML::Key << "detector_params" << YAML::Value << YAML::Flow << YAML::BeginMap;
  for (const auto &kv : ev.detector_params) out << YAML::Key << kv.first << YAML::Value << kv.second;
  out << YAML::EndMap;
  out << YAML::Key << "objects" << YAML::Value << YAML::BeginSeq;
  for (const auto &o : ev.objects) {
    out << YAML::Flow << YAML::BeginMap;
    out << YAML::Key << "type" << YAML::Value << o.type;
    out << YAML::Key << "type_code" << YAML::Value << o.type_code;
    out << YAML::Key << "p4" << YAML::Value;
    emitP4(out, o.p4);
    out << YAML::Key << "charge" << YAML::Value << o.charge;
    out << YAML::Key << "features" << YAML::Value << YAML::Flow << o.features;
    out << YAML::Key << "truth_index" << YAML::Value << o.truth_index;
    out << YAML::Key << "true_type" << YAML::Value << o.true_type;
    out << YAML::Key << "is_fake" << YAML::Value << o.is_fake;
    out << YAML::Key << "n_constituents" << YAML::Value << o.n_constituents;
    out << YAML::Key << "n_charged" << YAML::Value << o.n_charged;
    out << YAML::EndMap;
  }
  out << YAML::EndSeq << YAML::EndMap;
  return std::string(out.c_str()) + "\n";
}

ObservedEvent observedEventFromYAML(const std::string &yaml) {
  const YAML::Node n = YAML::Load(yaml);
  ObservedEvent ev;
  ev.task_id = n["task_id"].as<std::string>();
  ev.env_id = n["env_id"].as<std::string>();
  ev.event_index = n["event_index"].as<uint64_t>();
  ev.truth_seed = n["truth_seed"].as<uint64_t>();
  ev.detector_seed = n["detector_seed"].as<uint64_t>();
  ev.component = n["component"].as<int>();
  ev.process = n["process"].as<std::string>();
  const auto met = n["met"].as<std::vector<double>>();
  ev.met_x = met.at(0);
  ev.met_y = met.at(1);
  ev.sum_et = n["sum_et"].as<double>();
  ev.n_lost = n["n_lost"].as<int>();
  ev.n_fakes = n["n_fakes"].as<int>();
  ev.n_misid = n["n_misid"].as<int>();
  ev.feature_names = std::make_shared<std::vector<std::string>>(n["feature_names"].as<std::vector<std::string>>());
  for (const auto &kv : n["detector_params"]) ev.detector_params.emplace_back(kv.first.as<std::string>(), kv.second.as<double>());
  for (const auto &on : n["objects"]) {
    ObservedObject o;
    o.type = on["type"].as<std::string>();
    o.type_code = on["type_code"].as<int>();
    o.p4 = parseP4(on["p4"]);
    o.charge = on["charge"].as<double>();
    o.features = on["features"].as<std::vector<double>>();
    o.truth_index = on["truth_index"].as<int>();
    o.true_type = on["true_type"].as<std::string>();
    o.is_fake = on["is_fake"].as<bool>();
    if (on["n_constituents"]) o.n_constituents = on["n_constituents"].as<int>();
    if (on["n_charged"]) o.n_charged = on["n_charged"].as<int>();
    ev.objects.push_back(std::move(o));
  }
  return ev;
}

} // namespace kg
