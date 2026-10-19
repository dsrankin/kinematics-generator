// Auxiliary per-object features.
//
// A feature is a named scalar attached to every observed object. Feature
// generators can be pure noise, noisy proxies for truth or observed
// quantities, or deliberately (spuriously) correlated with event labels,
// truth categories or detector conditions. Each task samples its own schema:
// which features are present (presence_prob) and their parameters.
#ifndef KG_FEATURES_H
#define KG_FEATURES_H

#include <map>
#include <set>
#include <string>
#include <vector>

#include "kg/Param.h"

namespace kg {

enum class FeatureKind {
  Gaussian,        // mean, sigma
  Uniform,         // min, max
  TruthProxy,      // quantity of the matched truth particle, with noise
  ObservedProxy,   // quantity of the observed object, with noise
  LabelCorrelated, // base_mean + shift[category] + sigma * N(0,1)
  Environment      // offset + coeff * <detector parameter> + sigma * N(0,1)
};

FeatureKind featureKindFromString(const std::string &s);
std::string featureKindName(FeatureKind k);

// Feature generator as written in a prior (numeric fields are parameters).
struct FeatureDef {
  std::string name;
  FeatureKind kind = FeatureKind::Gaussian;
  std::vector<std::string> applies_to; // object types; empty = all
  std::string quantity;                // proxies
  std::string label;                   // label_correlated: component | process | truth_type | parent | is_fake
  std::string fake_policy = "missing"; // truth_proxy for fakes: missing | observed
  std::string env_param;               // environment feature: detector parameter name
  Distribution shift;                  // label_correlated: per-category shift, drawn per task
  std::map<std::string, int> num;      // field -> parameter index
};

struct FeaturePrior {
  ParamSet params;
  std::vector<FeatureDef> defs;
};

FeaturePrior parseFeaturePrior(const YAML::Node &node);

// Feature generator as realised in a task (all numbers fixed).
struct FeatureSpec {
  std::string name;
  FeatureKind kind = FeatureKind::Gaussian;
  std::set<std::string> applies_to;
  std::string quantity, label, fake_policy = "missing", env_param;
  std::map<std::string, double> p;      // numeric fields
  std::map<std::string, double> shifts; // label_correlated: category -> shift

  double get(const std::string &k, double def = 0) const {
    auto it = p.find(k);
    return it == p.end() ? def : it->second;
  }
  bool applies(const std::string &type) const { return applies_to.empty() || applies_to.count(type) > 0; }
};

void emitFeatureSpec(YAML::Emitter &out, const FeatureSpec &f);
FeatureSpec parseFeatureSpec(const YAML::Node &node);

} // namespace kg

#endif
