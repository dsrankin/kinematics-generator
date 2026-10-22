#include "kg/Detector.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace kg {

namespace {

const double kInf = std::numeric_limits<double>::infinity();

std::shared_ptr<const ResponseMap> parseResponseMap(const YAML::Node &n, const std::string &ctx) {
  auto m = std::make_shared<ResponseMap>();
  for (const auto &kv : n)
    if (!std::set<std::string>{"eta_edges", "pt_edges", "scale", "resolution", "efficiency"}.count(kv.first.as<std::string>()))
      throw std::invalid_argument(ctx + ": unknown key '" + kv.first.as<std::string>() + "'");
  if (!n["eta_edges"] || !n["pt_edges"] || !n["scale"] || !n["resolution"])
    throw std::invalid_argument(ctx + ": needs eta_edges, pt_edges, scale and resolution");
  m->eta_edges = n["eta_edges"].as<std::vector<double>>();
  m->pt_edges = n["pt_edges"].as<std::vector<double>>();
  auto increasing = [&](const std::vector<double> &e, const char *what) {
    if (e.size() < 2) throw std::invalid_argument(ctx + ": " + what + " needs at least two edges");
    for (size_t i = 1; i < e.size(); ++i)
      if (!(e[i] > e[i - 1])) throw std::invalid_argument(ctx + ": " + what + " must increase");
  };
  increasing(m->eta_edges, "eta_edges");
  increasing(m->pt_edges, "pt_edges");
  const size_t ne = m->eta_edges.size() - 1, np = m->pt_edges.size() - 1;
  auto table = [&](const char *key, double lo, double hi, bool hiInclusive) {
    auto t = n[key].as<std::vector<std::vector<double>>>();
    if (t.size() != ne) throw std::invalid_argument(ctx + ": " + key + " needs " + std::to_string(ne) + " eta rows");
    for (const auto &row : t) {
      if (row.size() != np) throw std::invalid_argument(ctx + ": " + key + " needs " + std::to_string(np) + " pT columns");
      for (double v : row)
        if (!(v >= lo) || (hiInclusive ? !(v <= hi) : !(v < hi)))
          throw std::invalid_argument(ctx + ": " + key + " value " + std::to_string(v) + " out of range");
    }
    return t;
  };
  m->scale = table("scale", 1e-6, kInf, false);
  m->resolution = table("resolution", 0.0, kInf, false);
  if (n["efficiency"]) m->efficiency = table("efficiency", 0.0, 1.0, true);
  return m;
}

std::shared_ptr<const MassResponseMap> parseMassResponse(const YAML::Node &n, const std::string &ctx) {
  for (const auto &kv : n)
    if (!std::set<std::string>{"pt_edges", "scale", "resolution"}.count(kv.first.as<std::string>()))
      throw std::invalid_argument(ctx + ": unknown key '" + kv.first.as<std::string>() + "'");
  if (!n["pt_edges"] || !n["scale"] || !n["resolution"]) throw std::invalid_argument(ctx + ": needs pt_edges, scale and resolution");
  auto m = std::make_shared<MassResponseMap>();
  m->pt_edges = n["pt_edges"].as<std::vector<double>>();
  m->scale = n["scale"].as<std::vector<double>>();
  m->resolution = n["resolution"].as<std::vector<double>>();
  if (m->pt_edges.size() < 2) throw std::invalid_argument(ctx + ": pt_edges needs at least two edges");
  for (size_t i = 1; i < m->pt_edges.size(); ++i)
    if (!(m->pt_edges[i] > m->pt_edges[i - 1])) throw std::invalid_argument(ctx + ": pt_edges must increase");
  const size_t nb = m->pt_edges.size() - 1;
  if (m->scale.size() != nb || m->resolution.size() != nb)
    throw std::invalid_argument(ctx + ": scale and resolution need " + std::to_string(nb) + " entries");
  for (size_t i = 0; i < nb; ++i)
    if (!(m->scale[i] > 0) || !(m->resolution[i] >= 0)) throw std::invalid_argument(ctx + ": scale must be > 0, resolution >= 0");
  return m;
}

void onlyKeys(const YAML::Node &n, const std::set<std::string> &allowed, const std::string &ctx) {
  if (!n || !n.IsMap()) return;
  for (const auto &kv : n)
    if (!allowed.count(kv.first.as<std::string>()))
      throw std::invalid_argument(ctx + ": unknown key '" + kv.first.as<std::string>() + "'");
}

ResponseDef parseResponse(const YAML::Node &n, const std::string &type, ParamSet &ps) {
  const std::string p = "detector." + type + ".";
  ResponseDef r;
  const YAML::Node eff = n["efficiency"];
  if (eff && eff.IsMap() && !eff["dist"] && !eff["expr"] && !eff["value"]) {
    onlyKeys(eff, {"plateau", "pt50", "width"}, "detector." + type + ".efficiency");
    r.efficiency = ps.parseOr(eff["plateau"], p + "efficiency", 1.0);
    r.eff_pt50 = ps.parseOr(eff["pt50"], p + "eff_pt50", 0.0);
    r.eff_width = ps.parseOr(eff["width"], p + "eff_width", 0.0);
  } else {
    r.efficiency = ps.parseOr(eff, p + "efficiency", 1.0);
    r.eff_pt50 = ps.parseOr(n["eff_pt50"], p + "eff_pt50", 0.0);
    r.eff_width = ps.parseOr(n["eff_width"], p + "eff_width", 0.0);
  }
  r.pt_min = ps.parseOr(n["pt_min"], p + "pt_min", 0.0);
  r.eta_max = ps.parseOr(n["eta_max"], p + "eta_max", kInf);
  const YAML::Node res = n["resolution"];
  onlyKeys(res, {"stochastic", "constant", "noise", "pt_term"}, "detector." + type + ".resolution");
  r.stochastic = ps.parseOr(res ? res["stochastic"] : YAML::Node(YAML::NodeType::Undefined), p + "resolution.stochastic", 0.0);
  r.constant = ps.parseOr(res ? res["constant"] : YAML::Node(YAML::NodeType::Undefined), p + "resolution.constant", 0.0);
  r.noise = ps.parseOr(res ? res["noise"] : YAML::Node(YAML::NodeType::Undefined), p + "resolution.noise", 0.0);
  r.pt_term = ps.parseOr(res ? res["pt_term"] : YAML::Node(YAML::NodeType::Undefined), p + "resolution.pt_term", 0.0);
  r.scale = ps.parseOr(n["scale"], p + "scale", 1.0);
  r.sigma_eta = ps.parseOr(n["sigma_eta"], p + "sigma_eta", 0.0);
  r.sigma_phi = ps.parseOr(n["sigma_phi"], p + "sigma_phi", 0.0);
  r.resolution_scale = ps.parseOr(n["resolution_scale"], p + "resolution_scale", 1.0);
  if (n["mass_response"]) r.mass_map = parseMassResponse(n["mass_response"], "detector." + type + ".mass_response");
  if (n["mass_scale"] || n["mass_resolution"] || r.mass_map) {
    r.mass_scale = ps.parseOr(n["mass_scale"], p + "mass_scale", 1.0);
    // With a mass_response map, mass_resolution multiplies the map's resolution.
    r.mass_resolution = ps.parseOr(n["mass_resolution"], p + "mass_resolution", r.mass_map ? 1.0 : 0.0);
  }
  if (n["response_map"]) r.map = parseResponseMap(n["response_map"], "detector." + type + ".response_map");
  if (n["misid"]) {
    if (!n["misid"].IsMap()) throw std::invalid_argument("detector." + type + ".misid must be a map type -> probability");
    for (const auto &kv : n["misid"]) {
      const std::string t = kv.first.as<std::string>();
      if (t == type) throw std::invalid_argument("detector." + type + ".misid lists its own type");
      r.misid.emplace_back(t, ps.parse(kv.second, p + "misid." + t));
    }
  }
  static const std::set<std::string> known = {"efficiency", "eff_pt50", "eff_width", "pt_min", "eta_max", "resolution",
                                              "scale", "sigma_eta", "sigma_phi", "misid", "resolution_scale", "response_map",
                                              "mass_scale", "mass_resolution", "mass_response"};
  for (const auto &kv : n)
    if (!known.count(kv.first.as<std::string>()))
      throw std::invalid_argument("detector." + type + ": unknown key '" + kv.first.as<std::string>() + "'");
  return r;
}

void emitP(YAML::Emitter &out, const ParamSet &ps, const char *key, int idx) {
  out << YAML::Key << key << YAML::Value;
  emitParamDef(out, ps.def(idx));
}

void check(bool ok, const std::string &ctx, const std::string &msg) {
  if (!ok) throw std::invalid_argument(ctx + ": " + msg);
}

} // namespace

std::vector<std::string> DetectorPrior::objectTypes() const {
  std::set<std::string> s;
  for (const auto &kv : objects) {
    if (kv.first != "default") s.insert(kv.first);
    for (const auto &m : kv.second.misid) s.insert(m.first);
  }
  for (const auto &t : fakes.types) s.insert(t.first);
  for (const auto &j : jets) s.insert(j.type);
  return {s.begin(), s.end()};
}

DetectorPrior parseDetectorPrior(const YAML::Node &n) {
  DetectorPrior d;
  ParamSet &ps = d.params;
  if (n && !n.IsMap()) throw std::invalid_argument("'detector' must be a map");
  if (n && n["objects"]) {
    for (const auto &kv : n["objects"]) {
      const std::string type = kv.first.as<std::string>();
      d.objects[type] = parseResponse(kv.second, type, ps);
    }
  }
  if (n && n["jets"]) {
    if (!n["jets"].IsSequence()) throw std::invalid_argument("detector.jets must be a list of collections");
    std::set<std::string> types, consumed;
    for (const auto &jn : n["jets"]) {
      for (const auto &k : jn)
        if (!std::set<std::string>{"type", "inputs", "radius", "in_met"}.count(k.first.as<std::string>()))
          throw std::invalid_argument("detector.jets: unknown key '" + k.first.as<std::string>() + "'");
      ClusteringDef c;
      if (!jn["type"]) throw std::invalid_argument("detector.jets: each collection needs a 'type'");
      c.type = jn["type"].as<std::string>();
      const std::string ctx = "detector.jets." + c.type;
      if (!types.insert(c.type).second) throw std::invalid_argument(ctx + ": duplicate collection");
      c.inputs = jn["inputs"] ? jn["inputs"].as<std::vector<std::string>>() : std::vector<std::string>{"jet"};
      if (c.inputs.empty()) throw std::invalid_argument(ctx + ": 'inputs' must not be empty");
      if (jn["in_met"]) c.in_met = jn["in_met"].as<bool>();
      if (c.in_met)
        for (const auto &t : c.inputs)
          if (!consumed.insert(t).second) throw std::invalid_argument(ctx + ": input '" + t + "' is already consumed by another in_met collection");
      c.radius = ps.parseOr(jn["radius"], ctx + ".radius", 0.4);
      d.jets.push_back(std::move(c));
    }
  }
  const YAML::Node f = n ? n["fakes"] : YAML::Node(YAML::NodeType::Undefined);
  onlyKeys(f, {"rate", "types", "pt_min", "pt_scale", "eta_max"}, "detector.fakes");
  d.fakes.rate = ps.parseOr(f ? f["rate"] : YAML::Node(YAML::NodeType::Undefined), "detector.fakes.rate", 0.0);
  d.fakes.pt_min = ps.parseOr(f ? f["pt_min"] : YAML::Node(YAML::NodeType::Undefined), "detector.fakes.pt_min", 20.0);
  d.fakes.pt_scale = ps.parseOr(f ? f["pt_scale"] : YAML::Node(YAML::NodeType::Undefined), "detector.fakes.pt_scale", 10.0);
  d.fakes.eta_max = ps.parseOr(f ? f["eta_max"] : YAML::Node(YAML::NodeType::Undefined), "detector.fakes.eta_max", 2.5);
  if (f && f["types"]) {
    for (const auto &kv : f["types"]) {
      const std::string t = kv.first.as<std::string>();
      d.fakes.types.emplace_back(t, ps.parse(kv.second, "detector.fakes.types." + t));
    }
  }
  const YAML::Node m = n ? n["met"] : YAML::Node(YAML::NodeType::Undefined);
  onlyKeys(m, {"soft_sigma", "sumet_coeff", "recoil_scale", "recoil_resolution", "lost_to_met"}, "detector.met");
  d.met.soft_sigma = ps.parseOr(m ? m["soft_sigma"] : YAML::Node(YAML::NodeType::Undefined), "detector.met.soft_sigma", 0.0);
  d.met.sumet_coeff = ps.parseOr(m ? m["sumet_coeff"] : YAML::Node(YAML::NodeType::Undefined), "detector.met.sumet_coeff", 0.0);
  d.met.recoil_scale = ps.parseOr(m ? m["recoil_scale"] : YAML::Node(YAML::NodeType::Undefined), "detector.met.recoil_scale", 1.0);
  d.met.recoil_resolution = ps.parseOr(m ? m["recoil_resolution"] : YAML::Node(YAML::NodeType::Undefined), "detector.met.recoil_resolution", 0.0);
  d.met.lost_to_met = ps.parseOr(m ? m["lost_to_met"] : YAML::Node(YAML::NodeType::Undefined), "detector.met.lost_to_met", 1.0);
  if (n && n["correlations"]) {
    if (!n["correlations"].IsMap()) throw std::invalid_argument("detector.correlations must be a map name -> group");
    std::set<std::string> used;
    bool unclustered = false;
    for (const auto &kv : n["correlations"]) {
      CorrelationDef c;
      c.name = kv.first.as<std::string>();
      const YAML::Node g = kv.second;
      const std::string ctx = "detector.correlations." + c.name;
      for (const auto &k : g)
        if (!std::set<std::string>{"types", "scale_sigma", "resolution_correlation", "unclustered"}.count(k.first.as<std::string>()))
          throw std::invalid_argument(ctx + ": unknown key '" + k.first.as<std::string>() + "'");
      if (!g["types"] || !g["types"].IsSequence() || g["types"].size() == 0)
        throw std::invalid_argument(ctx + ": needs a non-empty 'types' list");
      c.types = g["types"].as<std::vector<std::string>>();
      for (const auto &t : c.types)
        if (!used.insert(t).second) throw std::invalid_argument(ctx + ": type '" + t + "' belongs to more than one correlation group");
      c.unclustered = g["unclustered"] && g["unclustered"].as<bool>();
      if (c.unclustered && unclustered) throw std::invalid_argument(ctx + ": only one correlation group may include the unclustered recoil");
      unclustered = unclustered || c.unclustered;
      c.scale_sigma = ps.parseOr(g["scale_sigma"], ctx + ".scale_sigma", 0.0);
      c.resolution_correlation = ps.parseOr(g["resolution_correlation"], ctx + ".resolution_correlation", 0.0);
      d.correlations.push_back(std::move(c));
    }
  }
  if (n && n["sort_by_pt"]) d.sort_by_pt = n["sort_by_pt"].as<bool>();
  if (n)
    for (const auto &kv : n)
      if (!std::set<std::string>{"objects", "jets", "fakes", "met", "correlations", "sort_by_pt"}.count(kv.first.as<std::string>()))
        throw std::invalid_argument("detector: unknown key '" + kv.first.as<std::string>() + "'");
  ps.compile();
  return d;
}

void emitDetectorPrior(YAML::Emitter &out, const DetectorPrior &d) {
  const ParamSet &ps = d.params;
  out << YAML::BeginMap;
  out << YAML::Key << "objects" << YAML::Value << YAML::BeginMap;
  for (const auto &kv : d.objects) {
    const ResponseDef &r = kv.second;
    out << YAML::Key << kv.first << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "efficiency" << YAML::Value << YAML::BeginMap;
    emitP(out, ps, "plateau", r.efficiency);
    emitP(out, ps, "pt50", r.eff_pt50);
    emitP(out, ps, "width", r.eff_width);
    out << YAML::EndMap;
    emitP(out, ps, "pt_min", r.pt_min);
    emitP(out, ps, "eta_max", r.eta_max);
    out << YAML::Key << "resolution" << YAML::Value << YAML::BeginMap;
    emitP(out, ps, "stochastic", r.stochastic);
    emitP(out, ps, "constant", r.constant);
    emitP(out, ps, "noise", r.noise);
    emitP(out, ps, "pt_term", r.pt_term);
    out << YAML::EndMap;
    emitP(out, ps, "scale", r.scale);
    emitP(out, ps, "sigma_eta", r.sigma_eta);
    emitP(out, ps, "sigma_phi", r.sigma_phi);
    emitP(out, ps, "resolution_scale", r.resolution_scale);
    if (r.mass_scale >= 0) emitP(out, ps, "mass_scale", r.mass_scale);
    if (r.mass_resolution >= 0) emitP(out, ps, "mass_resolution", r.mass_resolution);
    if (r.mass_map) {
      out << YAML::Key << "mass_response" << YAML::Value << YAML::BeginMap;
      out << YAML::Key << "pt_edges" << YAML::Value << YAML::Flow << r.mass_map->pt_edges;
      out << YAML::Key << "scale" << YAML::Value << YAML::Flow << r.mass_map->scale;
      out << YAML::Key << "resolution" << YAML::Value << YAML::Flow << r.mass_map->resolution;
      out << YAML::EndMap;
    }
    if (r.map) {
      out << YAML::Key << "response_map" << YAML::Value << YAML::BeginMap;
      out << YAML::Key << "eta_edges" << YAML::Value << YAML::Flow << r.map->eta_edges;
      out << YAML::Key << "pt_edges" << YAML::Value << YAML::Flow << r.map->pt_edges;
      auto tab = [&](const char *k, const std::vector<std::vector<double>> &t) {
        out << YAML::Key << k << YAML::Value << YAML::BeginSeq;
        for (const auto &row : t) out << YAML::Flow << row;
        out << YAML::EndSeq;
      };
      tab("scale", r.map->scale);
      tab("resolution", r.map->resolution);
      if (!r.map->efficiency.empty()) tab("efficiency", r.map->efficiency);
      out << YAML::EndMap;
    }
    if (!r.misid.empty()) {
      out << YAML::Key << "misid" << YAML::Value << YAML::BeginMap;
      for (const auto &m : r.misid) emitP(out, ps, m.first.c_str(), m.second);
      out << YAML::EndMap;
    }
    out << YAML::EndMap;
  }
  out << YAML::EndMap;
  if (!d.jets.empty()) {
    out << YAML::Key << "jets" << YAML::Value << YAML::BeginSeq;
    for (const auto &c : d.jets) {
      out << YAML::BeginMap;
      out << YAML::Key << "type" << YAML::Value << c.type;
      out << YAML::Key << "inputs" << YAML::Value << YAML::Flow << c.inputs;
      emitP(out, ps, "radius", c.radius);
      out << YAML::Key << "in_met" << YAML::Value << c.in_met;
      out << YAML::EndMap;
    }
    out << YAML::EndSeq;
  }
  out << YAML::Key << "fakes" << YAML::Value << YAML::BeginMap;
  emitP(out, ps, "rate", d.fakes.rate);
  emitP(out, ps, "pt_min", d.fakes.pt_min);
  emitP(out, ps, "pt_scale", d.fakes.pt_scale);
  emitP(out, ps, "eta_max", d.fakes.eta_max);
  if (!d.fakes.types.empty()) {
    out << YAML::Key << "types" << YAML::Value << YAML::BeginMap;
    for (const auto &t : d.fakes.types) emitP(out, ps, t.first.c_str(), t.second);
    out << YAML::EndMap;
  }
  out << YAML::EndMap;
  out << YAML::Key << "met" << YAML::Value << YAML::BeginMap;
  emitP(out, ps, "soft_sigma", d.met.soft_sigma);
  emitP(out, ps, "sumet_coeff", d.met.sumet_coeff);
  emitP(out, ps, "recoil_scale", d.met.recoil_scale);
  emitP(out, ps, "recoil_resolution", d.met.recoil_resolution);
  emitP(out, ps, "lost_to_met", d.met.lost_to_met);
  out << YAML::EndMap;
  if (!d.correlations.empty()) {
    out << YAML::Key << "correlations" << YAML::Value << YAML::BeginMap;
    for (const auto &c : d.correlations) {
      out << YAML::Key << c.name << YAML::Value << YAML::BeginMap;
      out << YAML::Key << "types" << YAML::Value << YAML::Flow << c.types;
      emitP(out, ps, "scale_sigma", c.scale_sigma);
      emitP(out, ps, "resolution_correlation", c.resolution_correlation);
      if (c.unclustered) out << YAML::Key << "unclustered" << YAML::Value << true;
      out << YAML::EndMap;
    }
    out << YAML::EndMap;
  }
  out << YAML::Key << "sort_by_pt" << YAML::Value << d.sort_by_pt;
  out << YAML::EndMap;
}

double Response::efficiencyAt(double pt) const {
  if (eff_width <= 0) return pt >= eff_pt50 ? efficiency : 0.0;
  return efficiency * 0.5 * (1.0 + std::erf((pt - eff_pt50) / (std::sqrt(2.0) * eff_width)));
}

size_t MassResponseMap::ptBin(double pt) const {
  const size_t i = std::upper_bound(pt_edges.begin(), pt_edges.end(), pt) - pt_edges.begin();
  return std::min(std::max<size_t>(i, 1), pt_edges.size() - 1) - 1;
}

size_t ResponseMap::etaBin(double a) const {
  const size_t i = std::upper_bound(eta_edges.begin(), eta_edges.end(), a) - eta_edges.begin();
  return std::min(std::max<size_t>(i, 1), eta_edges.size() - 1) - 1;
}

size_t ResponseMap::ptBin(double pt) const {
  const size_t i = std::upper_bound(pt_edges.begin(), pt_edges.end(), pt) - pt_edges.begin();
  return std::min(std::max<size_t>(i, 1), pt_edges.size() - 1) - 1;
}

void Response::lookup(double e, double pt, double eta, double &eff, double &mean, double &sigma) const {
  if (!map) {
    eff = efficiencyAt(pt);
    mean = 1.0;
    sigma = relativeResolution(e, pt);
    return;
  }
  const size_t ie = map->etaBin(std::fabs(eta)), ip = map->ptBin(pt);
  eff = map->efficiency.empty() ? efficiencyAt(pt) : efficiency * map->efficiency.at(ie).at(ip);
  mean = map->scale.at(ie).at(ip);
  sigma = resolution_scale * map->resolution.at(ie).at(ip);
}

double Response::relativeResolution(double e, double pt) const {
  double s2 = constant * constant + pt_term * pt_term * pt * pt;
  if (e > 0) s2 += stochastic * stochastic / e + noise * noise / (e * e);
  return std::sqrt(s2);
}

DetectorState buildDetectorState(const DetectorPrior &d, const std::vector<double> &v, const std::string &ctx) {
  DetectorState s;
  auto build = [&](const std::string &type, const ResponseDef &r) {
    const std::string c = ctx + " detector." + type;
    Response o;
    o.efficiency = v[r.efficiency];
    o.eff_pt50 = v[r.eff_pt50];
    o.eff_width = v[r.eff_width];
    o.pt_min = v[r.pt_min];
    o.eta_max = v[r.eta_max];
    o.stochastic = v[r.stochastic];
    o.constant = v[r.constant];
    o.noise = v[r.noise];
    o.pt_term = v[r.pt_term];
    o.scale = v[r.scale];
    o.sigma_eta = v[r.sigma_eta];
    o.sigma_phi = v[r.sigma_phi];
    o.resolution_scale = v[r.resolution_scale];
    o.map = r.map;
    o.mass_map = r.mass_map;
    if (r.mass_scale >= 0 || r.mass_resolution >= 0) {
      o.smear_mass = true;
      if (r.mass_scale >= 0) o.mass_scale = v[r.mass_scale];
      if (r.mass_resolution >= 0) o.mass_resolution = v[r.mass_resolution];
      check(o.mass_scale > 0 && o.mass_resolution >= 0, c, "mass_scale must be > 0 and mass_resolution >= 0");
    }
    check(o.resolution_scale >= 0, c, "resolution_scale must be >= 0");
    check(o.efficiency >= 0 && o.efficiency <= 1, c, "efficiency must be in [0, 1], got " + std::to_string(o.efficiency));
    check(o.eff_width >= 0, c, "efficiency width must be >= 0");
    check(o.eta_max > 0, c, "eta_max must be > 0");
    check(o.stochastic >= 0 && o.constant >= 0 && o.noise >= 0 && o.pt_term >= 0, c, "resolution terms must be >= 0");
    check(o.scale > 0, c, "scale must be > 0");
    check(o.sigma_eta >= 0 && o.sigma_phi >= 0, c, "angular resolutions must be >= 0");
    double tot = 0;
    for (const auto &m : r.misid) {
      const double p = v[m.second];
      check(p >= 0 && p <= 1, c, "misid probability to '" + m.first + "' must be in [0, 1]");
      tot += p;
      o.misid.emplace_back(m.first, p);
    }
    check(tot <= 1 + 1e-12, c, "misid probabilities sum to more than 1");
    return o;
  };
  for (const auto &kv : d.objects) {
    if (kv.first == "default") s.fallback = build(kv.first, kv.second);
    else s.objects[kv.first] = build(kv.first, kv.second);
  }
  for (const auto &c : d.jets) {
    Clustering j;
    j.type = c.type;
    j.inputs = c.inputs;
    j.radius = v[c.radius];
    j.in_met = c.in_met;
    check(j.radius > 0 && j.radius < 3, ctx + " detector.jets." + c.type, "radius must be in (0, 3)");
    s.jets.push_back(std::move(j));
  }
  s.fake_rate = v[d.fakes.rate];
  s.fake_pt_min = v[d.fakes.pt_min];
  s.fake_pt_scale = v[d.fakes.pt_scale];
  s.fake_eta_max = v[d.fakes.eta_max];
  check(s.fake_rate >= 0, ctx + " detector.fakes", "rate must be >= 0");
  check(s.fake_pt_scale > 0, ctx + " detector.fakes", "pt_scale must be > 0");
  double wsum = 0;
  for (const auto &t : d.fakes.types) {
    const double w = v[t.second];
    check(w >= 0, ctx + " detector.fakes", "type weights must be >= 0");
    wsum += w;
    s.fake_types.push_back(t.first);
    s.fake_weights.push_back(w);
  }
  check(s.fake_rate == 0 || wsum > 0, ctx + " detector.fakes", "rate > 0 requires 'types' with positive weights");
  s.soft_sigma = v[d.met.soft_sigma];
  s.sumet_coeff = v[d.met.sumet_coeff];
  s.recoil_scale = v[d.met.recoil_scale];
  s.recoil_resolution = v[d.met.recoil_resolution];
  s.lost_to_met = v[d.met.lost_to_met];
  check(s.soft_sigma >= 0 && s.sumet_coeff >= 0 && s.recoil_resolution >= 0, ctx + " detector.met", "resolutions must be >= 0");
  check(s.lost_to_met >= 0 && s.lost_to_met <= 1, ctx + " detector.met", "lost_to_met must be in [0, 1]");
  for (const auto &c : d.correlations) {
    CorrelationGroup g;
    g.name = c.name;
    g.types = c.types;
    g.unclustered = c.unclustered;
    g.scale_sigma = v[c.scale_sigma];
    g.rho = v[c.resolution_correlation];
    const std::string cc = ctx + " detector.correlations." + c.name;
    check(g.scale_sigma >= 0, cc, "scale_sigma must be >= 0");
    check(g.rho >= 0 && g.rho <= 1, cc, "resolution_correlation must be in [0, 1]");
    for (const auto &t : g.types) s.group_of[t] = static_cast<int>(s.groups.size());
    s.groups.push_back(std::move(g));
  }
  s.sort_by_pt = d.sort_by_pt;
  return s;
}

std::map<std::string, double> Environment::parameterValues() const {
  std::map<std::string, double> out;
  for (size_t i = 0; i < values.size(); ++i) out[detector->params.def(static_cast<int>(i)).name] = values[i];
  return out;
}

Environment sampleEnvironment(const DetectorPrior &prior, uint64_t seed, const std::string &env_id) {
  Environment env;
  env.env_id = env_id;
  env.seed = seed;
  Rng rng(mixSeed(seed, kStreamEnvironment));
  env.values = prior.params.sample(rng);
  auto frozen = std::make_shared<DetectorPrior>(prior);
  frozen->params = prior.params.frozen(env.values);
  env.detector = frozen;
  env.state = buildDetectorState(*frozen, env.values, "environment '" + env_id + "'");
  return env;
}

Environment sampleEnvironment(const DetectorPrior &prior, uint64_t master_seed, int64_t index) {
  return sampleEnvironment(prior, mixSeed(master_seed, kStreamEnvironment, static_cast<uint64_t>(index)),
                           "env_" + std::to_string(index));
}

void emitEnvironment(YAML::Emitter &out, const Environment &env) {
  out << YAML::BeginMap;
  out << YAML::Key << "env_id" << YAML::Value << env.env_id;
  out << YAML::Key << "seed" << YAML::Value << env.seed;
  out << YAML::Key << "values" << YAML::Value << YAML::Flow << env.parameterValues();
  out << YAML::Key << "detector" << YAML::Value;
  emitDetectorPrior(out, *env.detector);
  out << YAML::EndMap;
}

Environment parseEnvironment(const YAML::Node &n) {
  // The serialized detector holds fixed values (and event-scoped definitions),
  // so re-sampling with the stored seed reproduces the environment exactly.
  const DetectorPrior d = parseDetectorPrior(n["detector"]);
  Environment env = sampleEnvironment(d, n["seed"].as<uint64_t>(), n["env_id"].as<std::string>());
  // Keep the original initial draws of event-scoped parameters (metadata only).
  if (n["values"]) {
    for (const auto &kv : n["values"].as<std::map<std::string, double>>()) env.values.at(env.detector->params.at(kv.first)) = kv.second;
  }
  return env;
}

std::string environmentToYAML(const Environment &env) {
  YAML::Emitter out;
  out.SetDoublePrecision(17);
  emitEnvironment(out, env);
  return out.c_str();
}

Environment environmentFromYAML(const std::string &yaml) { return parseEnvironment(YAML::Load(yaml)); }

} // namespace kg
