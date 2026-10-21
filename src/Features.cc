#include "kg/Features.h"

#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace kg {

FeatureKind featureKindFromString(const std::string &s) {
  if (s == "gaussian") return FeatureKind::Gaussian;
  if (s == "uniform") return FeatureKind::Uniform;
  if (s == "truth_proxy") return FeatureKind::TruthProxy;
  if (s == "observed_proxy") return FeatureKind::ObservedProxy;
  if (s == "label_correlated") return FeatureKind::LabelCorrelated;
  if (s == "environment") return FeatureKind::Environment;
  throw std::invalid_argument("unknown feature kind '" + s + "'");
}

std::string featureKindName(FeatureKind k) {
  switch (k) {
    case FeatureKind::Gaussian: return "gaussian";
    case FeatureKind::Uniform: return "uniform";
    case FeatureKind::TruthProxy: return "truth_proxy";
    case FeatureKind::ObservedProxy: return "observed_proxy";
    case FeatureKind::LabelCorrelated: return "label_correlated";
    case FeatureKind::Environment: return "environment";
  }
  return "gaussian";
}

namespace {

// Numeric fields accepted per kind, with defaults.
const std::map<FeatureKind, std::map<std::string, double>> &numericFields() {
  static const std::map<FeatureKind, std::map<std::string, double>> f = {
      {FeatureKind::Gaussian, {{"mean", 0}, {"sigma", 1}}},
      {FeatureKind::Uniform, {{"min", 0}, {"max", 1}}},
      {FeatureKind::TruthProxy, {{"rel_sigma", 0}, {"abs_sigma", 0}, {"bias", 0}}},
      {FeatureKind::ObservedProxy, {{"rel_sigma", 0}, {"abs_sigma", 0}, {"bias", 0}}},
      {FeatureKind::LabelCorrelated, {{"base_mean", 0}, {"sigma", 1}}},
      {FeatureKind::Environment, {{"coeff", 1}, {"offset", 0}, {"sigma", 0}}},
  };
  return f;
}

const std::set<std::string> &truthQuantities() {
  static const std::set<std::string> q = {"pt", "eta", "phi", "m", "e", "px", "py", "pz", "charge", "pdg",
                                          "is_fake", "depth", "parent_pdg", "parent_mass", "branch"};
  return q;
}
const std::set<std::string> &observedQuantities() {
  static const std::set<std::string> q = {"pt", "eta", "phi", "m", "e", "px", "py", "pz", "charge", "n_constituents", "n_charged"};
  return q;
}

} // namespace

FeaturePrior parseFeaturePrior(const YAML::Node &node) {
  FeaturePrior fp;
  if (!node) {
    fp.params.compile();
    return fp;
  }
  if (!node.IsSequence()) throw std::invalid_argument("'features' must be a list");
  std::set<std::string> names;
  for (const auto &n : node) {
    FeatureDef d;
    if (!n["name"]) throw std::invalid_argument("feature without 'name'");
    d.name = n["name"].as<std::string>();
    if (!names.insert(d.name).second) throw std::invalid_argument("feature '" + d.name + "' defined twice");
    const std::string ctx = "feature '" + d.name + "'";
    if (!n["kind"]) throw std::invalid_argument(ctx + ": missing 'kind'");
    d.kind = featureKindFromString(n["kind"].as<std::string>());
    if (n["applies_to"]) d.applies_to = n["applies_to"].as<std::vector<std::string>>();
    const std::string pre = "feature." + d.name + ".";
    for (const auto &kv : numericFields().at(d.kind)) d.num[kv.first] = fp.params.parseOr(n[kv.first], pre + kv.first, kv.second);
    d.num["missing_prob"] = fp.params.parseOr(n["missing_prob"], pre + "missing_prob", 0.0);
    d.num["presence_prob"] = fp.params.parseOr(n["presence_prob"], pre + "presence_prob", 1.0);
    std::set<std::string> allowed = {"name", "kind", "applies_to", "missing_prob", "presence_prob"};
    for (const auto &kv : numericFields().at(d.kind)) allowed.insert(kv.first);
    if (d.kind == FeatureKind::TruthProxy || d.kind == FeatureKind::ObservedProxy) {
      if (!n["quantity"]) throw std::invalid_argument(ctx + ": missing 'quantity'");
      d.quantity = n["quantity"].as<std::string>();
      const auto &q = d.kind == FeatureKind::TruthProxy ? truthQuantities() : observedQuantities();
      if (!q.count(d.quantity)) throw std::invalid_argument(ctx + ": unknown quantity '" + d.quantity + "'");
      allowed.insert("quantity");
      if (d.kind == FeatureKind::TruthProxy) {
        if (n["fake_policy"]) d.fake_policy = n["fake_policy"].as<std::string>();
        if (d.fake_policy != "missing" && d.fake_policy != "observed")
          throw std::invalid_argument(ctx + ": fake_policy must be 'missing' or 'observed'");
        if (d.fake_policy == "observed" && !observedQuantities().count(d.quantity))
          throw std::invalid_argument(ctx + ": fake_policy 'observed' needs an observable quantity");
        allowed.insert("fake_policy");
      }
    }
    if (d.kind == FeatureKind::LabelCorrelated) {
      if (!n["label"]) throw std::invalid_argument(ctx + ": missing 'label'");
      d.label = n["label"].as<std::string>();
      static const std::set<std::string> labels = {"component", "process", "truth_type", "parent", "is_fake"};
      if (!labels.count(d.label)) throw std::invalid_argument(ctx + ": unknown label '" + d.label + "'");
      if (!n["shift"]) throw std::invalid_argument(ctx + ": missing 'shift' distribution");
      d.shift = n["shift"].IsScalar() ? Distribution::fixed(n["shift"].as<double>()) : Distribution::fromYAML(n["shift"], ctx + " shift");
      allowed.insert("label");
      allowed.insert("shift");
    }
    if (d.kind == FeatureKind::Environment) {
      if (!n["param"]) throw std::invalid_argument(ctx + ": missing 'param' (detector parameter name)");
      d.env_param = n["param"].as<std::string>();
      allowed.insert("param");
    }
    for (const auto &kv : n)
      if (!allowed.count(kv.first.as<std::string>()))
        throw std::invalid_argument(ctx + ": unknown key '" + kv.first.as<std::string>() + "'");
    fp.defs.push_back(std::move(d));
  }
  fp.params.compile();
  for (size_t i = 0; i < fp.params.size(); ++i)
    if (fp.params.isEvent(static_cast<int>(i)))
      throw std::invalid_argument("feature parameter '" + fp.params.def(static_cast<int>(i)).name + "' cannot be event-scoped");
  return fp;
}

void emitFeatureSpec(YAML::Emitter &out, const FeatureSpec &f) {
  out << YAML::BeginMap;
  out << YAML::Key << "name" << YAML::Value << f.name;
  out << YAML::Key << "kind" << YAML::Value << featureKindName(f.kind);
  if (!f.applies_to.empty())
    out << YAML::Key << "applies_to" << YAML::Value << YAML::Flow << std::vector<std::string>(f.applies_to.begin(), f.applies_to.end());
  if (!f.quantity.empty()) out << YAML::Key << "quantity" << YAML::Value << f.quantity;
  if (!f.label.empty()) out << YAML::Key << "label" << YAML::Value << f.label;
  if (f.kind == FeatureKind::TruthProxy) out << YAML::Key << "fake_policy" << YAML::Value << f.fake_policy;
  if (!f.env_param.empty()) out << YAML::Key << "param" << YAML::Value << f.env_param;
  out << YAML::Key << "values" << YAML::Value << YAML::Flow << f.p;
  if (!f.shifts.empty()) out << YAML::Key << "shifts" << YAML::Value << YAML::Flow << f.shifts;
  out << YAML::EndMap;
}

FeatureSpec parseFeatureSpec(const YAML::Node &n) {
  FeatureSpec f;
  f.name = n["name"].as<std::string>();
  f.kind = featureKindFromString(n["kind"].as<std::string>());
  if (n["applies_to"])
    for (const auto &t : n["applies_to"]) f.applies_to.insert(t.as<std::string>());
  if (n["quantity"]) f.quantity = n["quantity"].as<std::string>();
  if (n["label"]) f.label = n["label"].as<std::string>();
  if (n["fake_policy"]) f.fake_policy = n["fake_policy"].as<std::string>();
  if (n["param"]) f.env_param = n["param"].as<std::string>();
  if (n["values"]) f.p = n["values"].as<std::map<std::string, double>>();
  if (n["shifts"]) f.shifts = n["shifts"].as<std::map<std::string, double>>();
  return f;
}

} // namespace kg
