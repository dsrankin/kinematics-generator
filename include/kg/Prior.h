// PriorConfig: a distribution over physics tasks and detector environments.
//
//   Prior --sampleTask--> Task --generateTruthEvent--> TruthEvent --DetectorModel::apply--> ObservedEvent
//
// A prior holds a catalogue of particles (string names, explicit decay graphs
// with branching fractions and angular models), a catalogue of production
// processes (1-3 produced particles, possibly with choices and multiplicity
// distributions), the initial-state model, how tasks are assembled from
// processes, a detector prior and a list of auxiliary feature generators.
// See DESIGN.md for the YAML schema.
#ifndef KG_PRIOR_H
#define KG_PRIOR_H

#include <map>
#include <string>
#include <vector>

#include "kg/Detector.h"
#include "kg/Features.h"
#include "kg/Kinematics.h"
#include "kg/Param.h"

namespace kg {

struct AngularDef {
  AngularModel::Type type = AngularModel::Type::Isotropic;
  AngularModel::Axis axis = AngularModel::Axis::Helicity;
  std::vector<int> coeffs; // parameter indices
};

struct DecayModeDef {
  std::vector<std::string> daughters; // 2 or 3
  int br = -1;                        // parameter index
  AngularDef angular;
  ThreeBodyModel model = ThreeBodyModel::PhaseSpace;
  std::string reconstruct_as; // non-empty: visible descendants are reconstructed as one object of this type
};

struct ParticleDef {
  std::string name;
  int pdg = 0;
  double charge = 0;
  int mass = -1, width = -1; // parameter indices ("<name>.mass", "<name>.width")
  bool stable = true;
  std::string object; // reconstructed object type, or "invisible"
  std::vector<DecayModeDef> modes;
  bool normalize_br = false; // explicitly allow branching fractions that do not sum to one
  bool mode_per_task = false; // decay_choice: task -> one mode is fixed per task
};

struct SlotDef {
  std::vector<std::string> choices; // one entry = fixed particle
  std::vector<double> weights;      // empty = equal
};

// One production alternative of a process (e.g. W, W + q, W + q q).
struct SubprocessDef {
  std::vector<SlotDef> slots;
  std::vector<int> multiplicity;              // optional: number of copies of the single slot
  std::vector<double> multiplicity_weights;
  AngularDef angular;                         // production angle (axis defaults to beam)
  int weight = -1;                            // parameter index of the per-event weight (-1: weight 1)
};

struct ProcessDef {
  std::string name;
  double weight = 1;                          // task-level sampling weight of the process
  std::vector<SubprocessDef> subprocesses;    // >= 1; chosen per event by weight
};

enum class SqrtSModel { ExponentialPdf, ThresholdExponential, Fixed, Uniform, Pdf };
enum class RecoilModel { None, Exponential, Param };

// Parton-level initial state. All numeric fields are parameters named
// "initial_state.<field>" so that they can be fixed, randomised per task or per event.
struct InitialStateDef {
  SqrtSModel model = SqrtSModel::ExponentialPdf;
  bool rapidity_from_pdf = true;
  RecoilModel recoil = RecoilModel::None;
  bool emit_recoil = false;          // add the recoil jets as visible final-state objects
  std::string recoil_object = "jet";
  int max_attempts = 10000000;
  // parameter indices
  int ebeam = -1, tau = -1, tau_slope = -1, offset = -1, scale = -1, value = -1, min = -1, max = -1;
  int rapidity = -1, recoil_min = -1, recoil_mean = -1, recoil_pt = -1, recoil_rapidity = -1;
  int pdf_a = -1, pdf_b = -1; // pdf model: effective parton density f(x) ~ x^-a (1-x)^b
  int recoil_n_jets = -1; // number of recoil (ISR) jets; each draws its own pT, rapidity and azimuth
};

// Final-state parton shower (see Jets.h). Parameters "parton_shower.<field>".
struct ShowerDef {
  bool enabled = false;
  int alpha_s = -1, kt_min = -1, theta_max = -1, hadron_mass = -1, hadron_shift = -1; // parameter indices
};

struct TaskSamplingDef {
  Distribution n_components = Distribution::fixed(1);
  std::string fractions = "equal"; // equal | dirichlet | fixed
  double dirichlet_alpha = 1;
  std::vector<double> fixed_fractions;
  bool distinct_processes = true;
  bool reject_invalid = true; // false: an invalid sampled task is an error
  int max_attempts = 1000;
};

struct PriorConfig {
  std::string name = "prior";
  uint64_t seed = 0; // default master seed
  ParamSet params;
  std::vector<ParticleDef> particles;
  std::map<std::string, int> particle_index;
  std::vector<ProcessDef> processes;
  InitialStateDef initial;
  ShowerDef shower;
  TaskSamplingDef task;
  double bw_cutoff = 15; // truncate Breit-Wigner masses at +-bw_cutoff*width (MadGraph default; 0: physical limits only)
  DetectorPrior detector;
  FeaturePrior features;

  const ParticleDef &particle(const std::string &n) const;
  // All reconstructed object types (particles, detector, fakes), sorted, without "invisible".
  std::vector<std::string> objectTypes() const;
  // Static graph checks: unknown daughters, multiplicity, cycles, scopes, BR sums (when fixed).
  void validate() const;
};

PriorConfig loadPrior(const std::string &path);
PriorConfig parsePrior(const YAML::Node &node);
PriorConfig parsePriorString(const std::string &yaml);

// Helpers shared with task (de)serialisation.
std::string defaultObjectType(int pdg, bool visible);
AngularDef parseAngularDef(const YAML::Node &node, ParamSet &ps, const std::string &prefix, AngularModel::Axis def_axis);
void parseInitialStateFlags(const YAML::Node &node, InitialStateDef &is);
void registerInitialStateParams(const YAML::Node &node, ParamSet &ps);
void bindInitialStateParams(InitialStateDef &is, const ParamSet &ps);
void emitInitialStateFlags(YAML::Emitter &out, const InitialStateDef &is);
// `parton_shower: {alpha_s, kt_min, theta_max}` (or false / none): registers the parameters.
bool registerShowerParams(const YAML::Node &node, ParamSet &ps);
void bindShowerParams(ShowerDef &s, const ParamSet &ps); // disabled when the parameters are absent
std::string sqrtSModelName(SqrtSModel m);
std::string threeBodyModelName(ThreeBodyModel m);
ThreeBodyModel threeBodyModelFromString(const std::string &s);
AngularModel::Axis axisFromString(const std::string &s);
std::string axisName(AngularModel::Axis a);

} // namespace kg

#endif
