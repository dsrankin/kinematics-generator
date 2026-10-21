#include "kg/Prior.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <set>
#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace kg {

namespace {

void requireKeys(const YAML::Node &n, const std::set<std::string> &allowed, const std::string &ctx) {
  if (!n.IsMap()) return;
  for (const auto &kv : n) {
    const std::string k = kv.first.as<std::string>();
    if (!allowed.count(k)) throw std::invalid_argument(ctx + ": unknown key '" + k + "'");
  }
}

} // namespace

std::string defaultObjectType(int pdg, bool visible) {
  if (!visible) return "invisible";
  const int a = std::abs(pdg);
  if ((a >= 1 && a <= 4) || a == 21) return "jet";
  if (a == 5) return "bjet";
  if (a == 11) return "electron";
  if (a == 13) return "muon";
  if (a == 15) return "tau";
  if (a == 22) return "photon";
  if (a == 12 || a == 14 || a == 16) return "invisible";
  return "";
}

AngularModel::Axis axisFromString(const std::string &s) {
  if (s == "helicity") return AngularModel::Axis::Helicity;
  if (s == "beam") return AngularModel::Axis::Beam;
  throw std::invalid_argument("unknown angular axis '" + s + "' (helicity|beam)");
}
std::string axisName(AngularModel::Axis a) { return a == AngularModel::Axis::Beam ? "beam" : "helicity"; }

ThreeBodyModel threeBodyModelFromString(const std::string &s) {
  if (s == "phase_space") return ThreeBodyModel::PhaseSpace;
  if (s == "legacy_sequential") return ThreeBodyModel::LegacySequential;
  throw std::invalid_argument("unknown three-body model '" + s + "' (phase_space|legacy_sequential)");
}
std::string threeBodyModelName(ThreeBodyModel m) {
  return m == ThreeBodyModel::LegacySequential ? "legacy_sequential" : "phase_space";
}

std::string sqrtSModelName(SqrtSModel m) {
  switch (m) {
    case SqrtSModel::ExponentialPdf: return "exponential_pdf";
    case SqrtSModel::ThresholdExponential: return "threshold_exponential";
    case SqrtSModel::Fixed: return "fixed";
    case SqrtSModel::Uniform: return "uniform";
  }
  return "exponential_pdf";
}

AngularDef parseAngularDef(const YAML::Node &n, ParamSet &ps, const std::string &prefix, AngularModel::Axis def_axis) {
  AngularDef a;
  a.axis = def_axis;
  if (!n) return a;
  std::string type;
  if (n.IsScalar()) {
    type = n.as<std::string>();
  } else {
    requireKeys(n, {"type", "coeffs", "axis"}, prefix);
    type = n["type"] ? n["type"].as<std::string>() : "isotropic";
    if (n["axis"]) a.axis = axisFromString(n["axis"].as<std::string>());
  }
  if (type == "isotropic" || type == "flat") {
    a.type = AngularModel::Type::Isotropic;
  } else if (type == "polynomial") {
    a.type = AngularModel::Type::Polynomial;
    if (!n.IsMap() || !n["coeffs"] || !n["coeffs"].IsSequence() || n["coeffs"].size() == 0)
      throw std::invalid_argument(prefix + ": polynomial angular model needs a non-empty 'coeffs' list");
    for (size_t k = 0; k < n["coeffs"].size(); ++k) a.coeffs.push_back(ps.parse(n["coeffs"][k], prefix + ".c" + std::to_string(k)));
  } else if (type == "W" || type == "w_leptonic") {
    // myGenerator's "W" mode: f(cos) = 1 - cos about the rest-frame z axis.
    a.type = AngularModel::Type::Polynomial;
    a.axis = AngularModel::Axis::Beam;
    a.coeffs.push_back(ps.addFixed(prefix + ".c0", 1.0));
    a.coeffs.push_back(ps.addFixed(prefix + ".c1", -1.0));
  } else {
    throw std::invalid_argument(prefix + ": unknown angular model '" + type + "'");
  }
  return a;
}

void parseInitialStateFlags(const YAML::Node &n, InitialStateDef &is) {
  if (!n) return;
  requireKeys(n, {"sqrt_s", "rapidity", "recoil", "max_attempts"}, "initial_state");
  const YAML::Node s = n["sqrt_s"];
  if (s) {
    const std::string m = s["model"] ? s["model"].as<std::string>() : "exponential_pdf";
    if (m == "exponential_pdf") is.model = SqrtSModel::ExponentialPdf;
    else if (m == "threshold_exponential") is.model = SqrtSModel::ThresholdExponential;
    else if (m == "fixed") is.model = SqrtSModel::Fixed;
    else if (m == "uniform") is.model = SqrtSModel::Uniform;
    else throw std::invalid_argument("initial_state.sqrt_s: unknown model '" + m + "'");
    requireKeys(s, {"model", "ebeam", "tau", "tau_slope", "threshold_offset", "scale", "value", "min", "max"}, "initial_state.sqrt_s");
  }
  const YAML::Node y = n["rapidity"];
  if (y && y.IsScalar() && y.as<std::string>() == "from_pdf") is.rapidity_from_pdf = true;
  else if (y) is.rapidity_from_pdf = false;
  else is.rapidity_from_pdf = (is.model == SqrtSModel::ExponentialPdf);
  if (is.rapidity_from_pdf && is.model != SqrtSModel::ExponentialPdf)
    throw std::invalid_argument("initial_state.rapidity: 'from_pdf' requires sqrt_s model exponential_pdf");
  const YAML::Node r = n["recoil"];
  if (r) {
    if (r.IsScalar()) {
      if (r.as<std::string>() != "none") throw std::invalid_argument("initial_state.recoil: expected 'none' or a map");
      is.recoil = RecoilModel::None;
    } else {
      requireKeys(r, {"model", "min_pt", "mean_pt", "pt", "n_jets", "emit", "object", "rapidity"}, "initial_state.recoil");
      const std::string m = r["model"] ? r["model"].as<std::string>() : "exponential";
      if (m == "none") is.recoil = RecoilModel::None;
      else if (m == "exponential") is.recoil = RecoilModel::Exponential;
      else if (m == "param") is.recoil = RecoilModel::Param;
      else throw std::invalid_argument("initial_state.recoil: unknown model '" + m + "'");
      if (r["emit"]) is.emit_recoil = r["emit"].as<bool>();
      if (r["object"]) is.recoil_object = r["object"].as<std::string>();
      if (is.recoil == RecoilModel::Param && !r["pt"]) throw std::invalid_argument("initial_state.recoil: model 'param' needs 'pt'");
    }
  }
  if (n["max_attempts"]) is.max_attempts = n["max_attempts"].as<int>();
}

void registerInitialStateParams(const YAML::Node &n, ParamSet &ps) {
  const YAML::Node s = n ? n["sqrt_s"] : YAML::Node(YAML::NodeType::Undefined);
  auto f = [&](const YAML::Node &parent, const char *key, const std::string &name, double def) {
    ps.parseOr(parent ? parent[key] : YAML::Node(YAML::NodeType::Undefined), "initial_state." + name, def);
  };
  f(s, "ebeam", "ebeam", 6500.0);
  f(s, "tau", "tau", 0.05);
  f(s, "tau_slope", "tau_slope", 0.02);
  f(s, "threshold_offset", "threshold_offset", 0.0);
  f(s, "scale", "scale", 100.0);
  f(s, "value", "value", 1000.0);
  f(s, "min", "min", 0.0);
  f(s, "max", "max", 13000.0);
  const YAML::Node y = n ? n["rapidity"] : YAML::Node(YAML::NodeType::Undefined);
  if (y && !(y.IsScalar() && y.as<std::string>() == "from_pdf")) ps.parse(y, "initial_state.rapidity");
  else ps.addFixed("initial_state.rapidity", 0.0);
  const YAML::Node r = (n && n["recoil"] && n["recoil"].IsMap()) ? n["recoil"] : YAML::Node(YAML::NodeType::Undefined);
  f(r, "min_pt", "recoil_min_pt", 10.0);
  f(r, "mean_pt", "recoil_mean_pt", 15.0);
  f(r, "pt", "recoil_pt", 0.0);
  f(r, "n_jets", "recoil_n_jets", 1.0);
  if (r && r["rapidity"]) ps.parse(r["rapidity"], "initial_state.recoil_rapidity");
  else if (!r) ps.addFixed("initial_state.recoil_rapidity", 0.0); // unused without a recoil model
  else {
    ParamDef d;
    d.name = "initial_state.recoil_rapidity";
    d.kind = ParamDef::Kind::Dist;
    d.dist = Distribution::uniform(-2.5, 2.5);
    d.scope = Scope::Event;
    ps.add(d);
  }
}

void bindInitialStateParams(InitialStateDef &is, const ParamSet &ps) {
  is.ebeam = ps.at("initial_state.ebeam");
  is.tau = ps.at("initial_state.tau");
  is.tau_slope = ps.at("initial_state.tau_slope");
  is.offset = ps.at("initial_state.threshold_offset");
  is.scale = ps.at("initial_state.scale");
  is.value = ps.at("initial_state.value");
  is.min = ps.at("initial_state.min");
  is.max = ps.at("initial_state.max");
  is.rapidity = ps.at("initial_state.rapidity");
  is.recoil_min = ps.at("initial_state.recoil_min_pt");
  is.recoil_mean = ps.at("initial_state.recoil_mean_pt");
  is.recoil_pt = ps.at("initial_state.recoil_pt");
  is.recoil_rapidity = ps.at("initial_state.recoil_rapidity");
  is.recoil_n_jets = ps.at("initial_state.recoil_n_jets");
}

void emitInitialStateFlags(YAML::Emitter &out, const InitialStateDef &is) {
  out << YAML::BeginMap;
  out << YAML::Key << "sqrt_s_model" << YAML::Value << sqrtSModelName(is.model);
  out << YAML::Key << "rapidity_from_pdf" << YAML::Value << is.rapidity_from_pdf;
  const char *rec = is.recoil == RecoilModel::None ? "none" : (is.recoil == RecoilModel::Exponential ? "exponential" : "param");
  out << YAML::Key << "recoil_model" << YAML::Value << rec;
  out << YAML::Key << "emit_recoil" << YAML::Value << is.emit_recoil;
  out << YAML::Key << "recoil_object" << YAML::Value << is.recoil_object;
  out << YAML::Key << "max_attempts" << YAML::Value << is.max_attempts;
  out << YAML::EndMap;
}

const ParticleDef &PriorConfig::particle(const std::string &n) const {
  auto it = particle_index.find(n);
  if (it == particle_index.end()) throw std::invalid_argument("unknown particle '" + n + "'");
  return particles[it->second];
}

std::vector<std::string> PriorConfig::objectTypes() const {
  std::set<std::string> s;
  for (const auto &p : particles) {
    if (p.object != "invisible" && p.stable) s.insert(p.object);
    for (const auto &m : p.modes)
      if (!m.reconstruct_as.empty()) s.insert(m.reconstruct_as);
  }
  if (initial.emit_recoil) s.insert(initial.recoil_object);
  for (const auto &t : detector.objectTypes()) s.insert(t);
  return {s.begin(), s.end()};
}

void PriorConfig::validate() const {
  if (processes.empty()) throw std::invalid_argument("prior '" + name + "': no production processes defined");
  for (const auto &p : particles) {
    for (size_t im = 0; im < p.modes.size(); ++im) {
      const auto &m = p.modes[im];
      const std::string ctx = "particle '" + p.name + "' decay " + std::to_string(im);
      if (m.daughters.size() != 2 && m.daughters.size() != 3)
        throw std::invalid_argument(ctx + ": only two- and three-body decays are supported");
      for (const auto &d : m.daughters)
        if (!particle_index.count(d)) throw std::invalid_argument(ctx + ": unknown daughter '" + d + "'");
      if (params.isEvent(m.br)) throw std::invalid_argument(ctx + ": branching fraction cannot be event-scoped");
      for (int c : m.angular.coeffs)
        if (params.isEvent(c)) throw std::invalid_argument(ctx + ": angular coefficients cannot be event-scoped");
      if (m.model == ThreeBodyModel::LegacySequential && m.daughters.size() != 3)
        throw std::invalid_argument(ctx + ": legacy_sequential applies to three-body decays only");
    }
    // Branching fractions known now must already be consistent.
    bool allFixed = true;
    double sum = 0;
    for (const auto &m : p.modes) {
      const ParamDef &d = params.def(m.br);
      if (d.kind != ParamDef::Kind::Fixed) allFixed = false;
      else {
        if (d.value < 0 || (!p.normalize_br && d.value > 1))
          throw std::invalid_argument("particle '" + p.name + "': branching fraction " + std::to_string(d.value) + " outside [0, 1]");
        sum += d.value;
      }
    }
    if (allFixed && !p.modes.empty() && !p.normalize_br && std::fabs(sum - 1) > 1e-6) {
      std::ostringstream os;
      os << "particle '" << p.name << "': branching fractions sum to " << sum << " (set normalize_br: true to rescale explicitly)";
      throw std::invalid_argument(os.str());
    }
    if (allFixed && p.normalize_br && !p.modes.empty() && !(sum > 0))
      throw std::invalid_argument("particle '" + p.name + "': branching fractions sum to zero");
  }
  // Decay graph must be acyclic.
  std::vector<int> state(particles.size(), 0);
  std::function<void(int, std::vector<std::string> &)> visit = [&](int i, std::vector<std::string> &path) {
    path.push_back(particles[i].name);
    if (state[i] == 1) {
      std::string s;
      for (size_t k = 0; k < path.size(); ++k) s += (k ? " -> " : "") + path[k];
      throw std::invalid_argument("cyclic decay chain: " + s);
    }
    if (state[i] == 0) {
      state[i] = 1;
      for (const auto &m : particles[i].modes)
        for (const auto &d : m.daughters) visit(particle_index.at(d), path);
      state[i] = 2;
    }
    path.pop_back();
  };
  for (size_t i = 0; i < particles.size(); ++i) {
    std::vector<std::string> path;
    visit(static_cast<int>(i), path);
  }
  for (const auto &pr : processes) {
    const std::string ctx = "process '" + pr.name + "'";
    if (!(pr.weight > 0)) throw std::invalid_argument(ctx + ": weight must be > 0");
    if (pr.slots.empty()) throw std::invalid_argument(ctx + ": empty production");
    if (pr.multiplicity.empty()) {
      if (pr.slots.size() > 3) throw std::invalid_argument(ctx + ": at most three produced particles are supported");
    } else {
      if (pr.slots.size() != 1) throw std::invalid_argument(ctx + ": 'multiplicity' requires exactly one production entry");
      for (int m : pr.multiplicity)
        if (m < 1 || m > 3) throw std::invalid_argument(ctx + ": multiplicity values must be 1, 2 or 3");
    }
    for (const auto &s : pr.slots) {
      if (s.choices.empty()) throw std::invalid_argument(ctx + ": empty choice");
      if (!s.weights.empty() && s.weights.size() != s.choices.size())
        throw std::invalid_argument(ctx + ": choice weights and particles differ in length");
      for (const auto &c : s.choices)
        if (!particle_index.count(c)) throw std::invalid_argument(ctx + ": unknown particle '" + c + "'");
    }
    for (int c : pr.angular.coeffs)
      if (params.isEvent(c)) throw std::invalid_argument(ctx + ": angular coefficients cannot be event-scoped");
  }
  for (const auto &f : features.defs)
    if (f.kind == FeatureKind::Environment && !detector.params.has(f.env_param))
      throw std::invalid_argument("feature '" + f.name + "': unknown detector parameter '" + f.env_param + "'");
  task.n_components.validate("task.components");
  if (task.distinct_processes) {
    const Distribution &nc = task.n_components;
    double nmax = nc.value;
    if (nc.type == Distribution::Type::Categorical) nmax = *std::max_element(nc.values.begin(), nc.values.end());
    else if (nc.type != Distribution::Type::Fixed) nmax = nc.max;
    if (std::lround(nmax) > static_cast<long>(processes.size()))
      throw std::invalid_argument("task.components can reach " + std::to_string(std::lround(nmax)) +
                                  ", which exceeds the number of distinct processes (" + std::to_string(processes.size()) + ")");
  }
  if (task.fractions == "fixed" && task.n_components.type != Distribution::Type::Fixed)
    throw std::invalid_argument("task.fractions: a fixed fraction list requires a fixed number of components");
  if (task.fractions == "fixed" && task.fixed_fractions.size() != static_cast<size_t>(task.n_components.value))
    throw std::invalid_argument("task.fractions: list length differs from the number of components");
}

PriorConfig parsePrior(const YAML::Node &root) {
  if (!root || !root.IsMap()) throw std::invalid_argument("prior: top level must be a map");
  requireKeys(root, {"name", "seed", "parameters", "particles", "processes", "production", "initial_state", "task",
                     "bw_cutoff", "detector", "features"},
              "prior");
  PriorConfig p;
  if (root["name"]) p.name = root["name"].as<std::string>();
  if (root["seed"]) p.seed = root["seed"].as<uint64_t>();
  if (root["bw_cutoff"]) p.bw_cutoff = root["bw_cutoff"].as<double>();
  ParamSet &ps = p.params;
  if (root["parameters"]) {
    if (!root["parameters"].IsMap()) throw std::invalid_argument("'parameters' must be a map");
    for (const auto &kv : root["parameters"]) ps.parse(kv.second, kv.first.as<std::string>());
  }
  // Particles.
  const YAML::Node parts = root["particles"];
  if (!parts || !parts.IsMap()) throw std::invalid_argument("prior: 'particles' map is required");
  for (const auto &kv : parts) {
    const std::string name = kv.first.as<std::string>();
    const YAML::Node n = kv.second;
    const std::string ctx = "particle '" + name + "'";
    if (name.empty()) throw std::invalid_argument("particle with empty name");
    if (p.particle_index.count(name)) throw std::invalid_argument(ctx + " defined twice");
    requireKeys(n, {"pdg", "charge", "mass", "width", "stable", "visible", "object", "decays", "normalize_br", "decay_choice"}, ctx);
    ParticleDef d;
    d.name = name;
    d.pdg = n["pdg"] ? n["pdg"].as<int>() : 0;
    d.charge = n["charge"] ? n["charge"].as<double>() : 0.0;
    d.mass = ps.parseOr(n["mass"], name + ".mass", 0.0);
    d.width = ps.parseOr(n["width"], name + ".width", 0.0);
    const bool visible = n["visible"] ? n["visible"].as<bool>() : true;
    d.object = n["object"] ? n["object"].as<std::string>() : defaultObjectType(d.pdg, visible);
    if (!visible) d.object = "invisible";
    d.normalize_br = n["normalize_br"] && n["normalize_br"].as<bool>();
    if (n["decay_choice"]) {
      const std::string c = n["decay_choice"].as<std::string>();
      if (c != "event" && c != "task") throw std::invalid_argument(ctx + ": decay_choice must be 'event' or 'task'");
      d.mode_per_task = (c == "task");
    }
    const YAML::Node dec = n["decays"];
    if (dec) {
      if (!dec.IsSequence() || dec.size() == 0) throw std::invalid_argument(ctx + ": 'decays' must be a non-empty list");
      for (size_t i = 0; i < dec.size(); ++i) {
        const std::string pre = name + ".decay" + std::to_string(i);
        DecayModeDef m;
        if (dec[i].IsSequence()) {
          m.daughters = dec[i].as<std::vector<std::string>>();
          if (dec.size() != 1) throw std::invalid_argument(ctx + ": several decay modes need explicit 'br' values");
          m.br = ps.addFixed(pre + ".br", 1.0);
        } else {
          requireKeys(dec[i], {"daughters", "br", "angular", "model", "reconstruct_as"}, ctx + " decay " + std::to_string(i));
          if (!dec[i]["daughters"]) throw std::invalid_argument(ctx + ": decay without 'daughters'");
          m.daughters = dec[i]["daughters"].as<std::vector<std::string>>();
          if (!dec[i]["br"] && dec.size() != 1) throw std::invalid_argument(ctx + ": several decay modes need explicit 'br' values");
          m.br = ps.parseOr(dec[i]["br"], pre + ".br", 1.0);
          m.angular = parseAngularDef(dec[i]["angular"], ps, pre + ".angular", AngularModel::Axis::Helicity);
          if (dec[i]["model"]) m.model = threeBodyModelFromString(dec[i]["model"].as<std::string>());
          if (dec[i]["reconstruct_as"]) {
            m.reconstruct_as = dec[i]["reconstruct_as"].as<std::string>();
            if (m.reconstruct_as.empty() || m.reconstruct_as == "invisible" || m.reconstruct_as == "default")
              throw std::invalid_argument(ctx + ": reconstruct_as must name an object type");
          }
        }
        d.modes.push_back(std::move(m));
      }
    }
    d.stable = d.modes.empty();
    if (n["stable"] && n["stable"].as<bool>() != d.stable)
      throw std::invalid_argument(ctx + ": 'stable' contradicts the presence/absence of 'decays'");
    if (d.stable && d.object.empty())
      throw std::invalid_argument(ctx + ": cannot infer object type from pdg " + std::to_string(d.pdg) + "; set 'object' or 'visible: false'");
    if (!d.stable && d.object.empty()) d.object = "invisible";
    p.particle_index[name] = static_cast<int>(p.particles.size());
    p.particles.push_back(std::move(d));
  }
  // Processes.
  auto parseProcess = [&](const YAML::Node &n, const std::string &defname) {
    ProcessDef pr;
    pr.name = n["name"] ? n["name"].as<std::string>() : defname;
    const std::string ctx = "process '" + pr.name + "'";
    requireKeys(n, {"name", "weight", "production", "multiplicity", "angular"}, ctx);
    if (n["weight"]) pr.weight = n["weight"].as<double>();
    if (!n["production"] || !n["production"].IsSequence()) throw std::invalid_argument(ctx + ": 'production' list is required");
    for (const auto &e : n["production"]) {
      SlotDef s;
      if (e.IsScalar()) s.choices.push_back(e.as<std::string>());
      else {
        requireKeys(e, {"choice", "weights"}, ctx);
        if (!e["choice"]) throw std::invalid_argument(ctx + ": production entry must be a name or {choice: [...]}");
        s.choices = e["choice"].as<std::vector<std::string>>();
        if (e["weights"]) s.weights = e["weights"].as<std::vector<double>>();
        for (double w : s.weights)
          if (!(w >= 0)) throw std::invalid_argument(ctx + ": choice weights must be >= 0");
      }
      pr.slots.push_back(std::move(s));
    }
    if (n["multiplicity"]) {
      const YAML::Node m = n["multiplicity"];
      if (m.IsScalar()) pr.multiplicity = {m.as<int>()};
      else {
        pr.multiplicity = m["values"].as<std::vector<int>>();
        if (m["weights"]) pr.multiplicity_weights = m["weights"].as<std::vector<double>>();
        if (!pr.multiplicity_weights.empty() && pr.multiplicity_weights.size() != pr.multiplicity.size())
          throw std::invalid_argument(ctx + ": multiplicity weights and values differ in length");
      }
    }
    pr.angular = parseAngularDef(n["angular"], ps, "process." + pr.name + ".angular", AngularModel::Axis::Beam);
    return pr;
  };
  if (root["processes"] && root["production"]) throw std::invalid_argument("prior: give either 'processes' or 'production', not both");
  if (root["processes"]) {
    if (!root["processes"].IsSequence()) throw std::invalid_argument("'processes' must be a list");
    std::set<std::string> names;
    int i = 0;
    for (const auto &n : root["processes"]) {
      ProcessDef pr = parseProcess(n, "process" + std::to_string(i++));
      if (!names.insert(pr.name).second) throw std::invalid_argument("process '" + pr.name + "' defined twice");
      p.processes.push_back(std::move(pr));
    }
  } else if (root["production"]) {
    YAML::Node n;
    n["name"] = "default";
    n["production"] = root["production"];
    p.processes.push_back(parseProcess(n, "default"));
  }
  // Initial state.
  parseInitialStateFlags(root["initial_state"], p.initial);
  registerInitialStateParams(root["initial_state"], ps);
  // Task assembly.
  const YAML::Node t = root["task"];
  if (t) {
    requireKeys(t, {"components", "fractions", "distinct_processes", "validation", "max_attempts"}, "task");
    if (t["components"]) {
      if (t["components"].IsScalar()) p.task.n_components = Distribution::fixed(t["components"].as<double>());
      else p.task.n_components = Distribution::fromYAML(t["components"], "task.components");
    }
    if (t["fractions"]) {
      const YAML::Node f = t["fractions"];
      if (f.IsScalar()) {
        p.task.fractions = f.as<std::string>();
        if (p.task.fractions != "equal") throw std::invalid_argument("task.fractions: expected 'equal', a list or {dirichlet: alpha}");
      } else if (f.IsSequence()) {
        p.task.fractions = "fixed";
        p.task.fixed_fractions = f.as<std::vector<double>>();
      } else if (f["dirichlet"]) {
        p.task.fractions = "dirichlet";
        p.task.dirichlet_alpha = f["dirichlet"].as<double>();
        if (!(p.task.dirichlet_alpha > 0)) throw std::invalid_argument("task.fractions: dirichlet alpha must be > 0");
      } else {
        throw std::invalid_argument("task.fractions: expected 'equal', a list or {dirichlet: alpha}");
      }
    }
    if (t["distinct_processes"]) p.task.distinct_processes = t["distinct_processes"].as<bool>();
    if (t["validation"]) {
      const std::string v = t["validation"].as<std::string>();
      if (v != "reject" && v != "error") throw std::invalid_argument("task.validation must be 'reject' or 'error'");
      p.task.reject_invalid = (v == "reject");
    }
    if (t["max_attempts"]) p.task.max_attempts = t["max_attempts"].as<int>();
  }
  ps.compile();
  bindInitialStateParams(p.initial, ps);
  p.detector = parseDetectorPrior(root["detector"]);
  p.features = parseFeaturePrior(root["features"]);
  p.validate();
  return p;
}

PriorConfig parsePriorString(const std::string &yaml) { return parsePrior(YAML::Load(yaml)); }

PriorConfig loadPrior(const std::string &path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open prior file '" + path + "'");
  try {
    return parsePrior(YAML::LoadFile(path));
  } catch (const YAML::Exception &e) {
    throw std::invalid_argument("prior '" + path + "': " + e.what());
  }
}

} // namespace kg
