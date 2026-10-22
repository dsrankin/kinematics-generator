// Detector priors and sampled detector environments.
//
// A DetectorPrior describes, per reconstructed object type, efficiency,
// acceptance, resolution, scale and mis-identification, plus fakes and
// missing-transverse-momentum reconstruction. Every number is a parameter
// (see Param.h), so it can be fixed, drawn once per environment
// (scope: task/environment, the default) or redrawn per event (scope: event).
//
// Sampling a DetectorPrior gives an Environment: a concrete detector that can
// be applied to any number of truth events (DetectorModel::apply).
#ifndef KG_DETECTOR_H
#define KG_DETECTOR_H

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "kg/Param.h"

namespace kg {

// Binned response in (|eta|, pT) of the true object, e.g. measured from a full simulation:
// mean response (reco/true pT), relative resolution and, optionally, efficiency per bin.
// Values outside the edges use the nearest bin. Arrays are [eta bin][pT bin].
struct ResponseMap {
  std::vector<double> eta_edges, pt_edges;
  std::vector<std::vector<double>> scale, resolution, efficiency; // efficiency may be empty

  size_t etaBin(double abs_eta) const;
  size_t ptBin(double pt) const;
};

// Binned jet-mass response in true pT: median reco/true mass and relative mass resolution.
struct MassResponseMap {
  std::vector<double> pt_edges, scale, resolution;
  size_t ptBin(double pt) const;
};

// Parameter indices describing one object type's response.
struct ResponseDef {
  int efficiency = -1, eff_pt50 = -1, eff_width = -1;
  int pt_min = -1, eta_max = -1;
  int stochastic = -1, constant = -1, noise = -1, pt_term = -1;
  int scale = -1, sigma_eta = -1, sigma_phi = -1;
  int resolution_scale = -1; // multiplies the map's resolution (domain randomisation around a measurement)
  int mass_scale = -1, mass_resolution = -1; // optional extra mass response (-1: mass scales with pT only)
  std::vector<std::pair<std::string, int>> misid; // reconstructed-as type -> probability
  std::shared_ptr<const ResponseMap> map;         // optional binned response
  std::shared_ptr<const MassResponseMap> mass_map; // optional binned mass response
};

// One jet collection: anti-kT clustering of the reconstruction candidates of the `inputs` types
// (after the parton shower) into objects of `type`. A collection with in_met consumes its inputs
// (they are replaced by the jets and enter MET); otherwise it is an extra view of the same
// candidates (e.g. large-radius jets) that is reconstructed independently and not added to MET.
struct ClusteringDef {
  std::string type;
  std::vector<std::string> inputs; // default: {"jet"}
  int radius = -1;                 // parameter "detector.jets.<type>.radius"
  bool in_met = true;
};

struct FakeDef {
  int rate = -1, pt_min = -1, pt_scale = -1, eta_max = -1;
  std::vector<std::pair<std::string, int>> types; // type -> relative weight
};

struct MetDef {
  int soft_sigma = -1, sumet_coeff = -1, recoil_scale = -1, recoil_resolution = -1, lost_to_met = -1;
};

// Correlated response of several object types within one event: a common
// relative scale shift N(0, scale_sigma) and a shared component of the
// resolution fluctuation, g = sqrt(rho) g_common + sqrt(1 - rho) g_object.
struct CorrelationDef {
  std::string name;
  std::vector<std::string> types;
  bool unclustered = false; // also shift the unclustered recoil in MET
  int scale_sigma = -1, resolution_correlation = -1;
};

struct DetectorPrior {
  ParamSet params;
  std::vector<CorrelationDef> correlations;
  std::map<std::string, ResponseDef> objects; // may contain "default"
  std::vector<ClusteringDef> jets;            // empty: every candidate is its own object
  FakeDef fakes;
  MetDef met;
  bool sort_by_pt = true;

  // Object types this detector can output (excluding "default").
  std::vector<std::string> objectTypes() const;
};

DetectorPrior parseDetectorPrior(const YAML::Node &node);
void emitDetectorPrior(YAML::Emitter &out, const DetectorPrior &d);

// Concrete response of one object type.
struct Response {
  double efficiency = 1, eff_pt50 = 0, eff_width = 0;
  double pt_min = 0, eta_max = 1e300;
  double stochastic = 0, constant = 0, noise = 0, pt_term = 0;
  double scale = 1, sigma_eta = 0, sigma_phi = 0;
  double resolution_scale = 1;
  bool smear_mass = false; // mass_scale / mass_resolution given
  double mass_scale = 1, mass_resolution = 0;
  std::vector<std::pair<std::string, double>> misid;
  std::shared_ptr<const ResponseMap> map;
  std::shared_ptr<const MassResponseMap> mass_map;

  double efficiencyAt(double pt) const;
  double relativeResolution(double e, double pt) const;
  // With a map: efficiency = plateau x map efficiency (when given), mean = map scale,
  // resolution = resolution_scale x map resolution; otherwise the parametric formulas.
  void lookup(double e, double pt, double eta, double &eff, double &mean, double &sigma) const;
};

struct CorrelationGroup {
  std::string name;
  std::vector<std::string> types;
  bool unclustered = false;
  double scale_sigma = 0, rho = 0;
};

struct Clustering {
  std::string type;
  std::vector<std::string> inputs;
  double radius = 0.4;
  bool in_met = true;
};

struct DetectorState {
  std::map<std::string, Response> objects;
  std::vector<Clustering> jets;
  std::vector<CorrelationGroup> groups;
  std::map<std::string, int> group_of; // object type -> index into groups
  Response fallback; // "default" entry, or a perfect detector
  double fake_rate = 0, fake_pt_min = 20, fake_pt_scale = 10, fake_eta_max = 2.5;
  std::vector<std::string> fake_types;
  std::vector<double> fake_weights;
  double soft_sigma = 0, sumet_coeff = 0, recoil_scale = 1, recoil_resolution = 0, lost_to_met = 1;
  bool sort_by_pt = true;

  const Response &response(const std::string &type) const {
    auto it = objects.find(type);
    return it == objects.end() ? fallback : it->second;
  }
};

DetectorState buildDetectorState(const DetectorPrior &d, const std::vector<double> &values, const std::string &context);

struct Environment {
  std::string env_id;
  uint64_t seed = 0;
  std::shared_ptr<const DetectorPrior> detector; // frozen: only event-scoped parameters remain random
  std::vector<double> values;
  DetectorState state;

  double value(const std::string &param) const { return values.at(detector->params.at(param)); }
  bool hasEventParams() const { return detector && detector->params.hasEventParams(); }
  // Parameter name -> value (event-scoped entries show this environment's initial draw).
  std::map<std::string, double> parameterValues() const;
};

// Draw an environment. `seed` fully determines it.
Environment sampleEnvironment(const DetectorPrior &prior, uint64_t seed, const std::string &env_id);
// Convenience: environment `index` of a domain-randomization sequence started from `master_seed`.
Environment sampleEnvironment(const DetectorPrior &prior, uint64_t master_seed, int64_t index);

std::string environmentToYAML(const Environment &env);
Environment environmentFromYAML(const std::string &yaml);
void emitEnvironment(YAML::Emitter &out, const Environment &env);
Environment parseEnvironment(const YAML::Node &node);

} // namespace kg

#endif
