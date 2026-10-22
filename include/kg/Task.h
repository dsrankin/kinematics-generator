// TaskConfig: one fully specified physics task sampled from a PriorConfig.
//
// All task-level randomness (masses, widths, branching fractions, process
// choice, multiplicities, feature schema, default detector environment) is
// resolved. Event-scoped parameters keep their distributions and are redrawn
// per event. A task is self-contained and can be serialised to YAML and
// reconstructed exactly (taskToYAML / taskFromYAML).
#ifndef KG_TASK_H
#define KG_TASK_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "kg/Prior.h"

namespace kg {

// Thrown when a sampled task violates kinematic or consistency constraints.
// sampleTask() rejects and redraws such tasks when task.validation is "reject".
class InvalidTask : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

struct TaskDecayMode {
  std::vector<int> daughters; // species indices
  double br = 1;
  AngularModel angular;
  ThreeBodyModel model = ThreeBodyModel::PhaseSpace;
  std::string reconstruct_as; // composite reconstruction type (empty: none)
};

struct TaskParticle {
  std::string name;
  int pdg = 0;
  double charge = 0;
  int mass = -1, width = -1; // parameter indices
  bool stable = true;
  std::string object;
  std::vector<TaskDecayMode> modes;
  std::vector<double> br_weights; // == modes[i].br
  int fixed_mode = -1;            // >= 0 when the decay mode is fixed for the task
};

struct TaskSubprocess {
  std::vector<int> produced; // species indices (1-3)
  double weight = 1;         // relative per-event probability within the component
  AngularModel angular;      // production angle in the hard-system rest frame
};

struct TaskComponent {
  std::string process;
  double fraction = 1;
  std::vector<TaskSubprocess> subprocesses; // >= 1

  std::vector<double> weights() const {
    std::vector<double> w;
    for (const auto &s : subprocesses) w.push_back(s.weight);
    return w;
  }
};

struct TaskConfig {
  std::string task_id;
  std::string prior_name;
  uint64_t seed = 0;        // task seed: all event seeds derive from it
  uint64_t master_seed = 0; // seed and index this task was sampled with
  int64_t task_index = 0;
  int attempts = 1;         // number of draws needed to obtain a valid task

  ParamSet params;            // frozen: only event-scoped parameters are random
  std::vector<double> values; // task-level values (initial draw for event-scoped ones)

  std::vector<TaskParticle> particles;
  std::map<std::string, int> particle_index;
  std::vector<TaskComponent> components;
  InitialStateDef initial;
  ShowerDef shower;
  double bw_cutoff = 0;

  std::vector<FeatureSpec> features;
  std::shared_ptr<const std::vector<std::string>> feature_names;
  std::vector<std::string> object_types; // vocabulary for ObservedObject::type_code

  Environment environment; // default detector environment for this task

  uint64_t next_event = 0; // used by generateTruthEvent(TaskConfig&)

  int species(const std::string &name) const;
  double mass(int species) const { return values[particles[species].mass]; }
  double width(int species) const { return values[particles[species].width]; }
  double param(const std::string &name) const { return values.at(params.at(name)); }
  std::map<std::string, double> parameterValues() const;
  int typeCode(const std::string &type) const;
  // Possible decay topologies, e.g. "ttbar: t>(W>{(l,v)|(q,q)},b) + t>(...)".
  std::string topology() const;
};

// Sample task `index` from `prior`, deterministically in (master_seed, index).
TaskConfig sampleTask(const PriorConfig &prior, uint64_t master_seed, int64_t index);
// Task 0 for the prior's own seed.
TaskConfig sampleTask(const PriorConfig &prior);

// Re-check a task's kinematic consistency (throws InvalidTask).
void validateTask(const TaskConfig &task);

// Sequential sampler: task i of the stream is sampleTask(prior, seed, i).
class TaskSampler {
 public:
  TaskSampler(const PriorConfig &prior, uint64_t master_seed, int64_t first_index = 0)
      : prior_(prior), seed_(master_seed), next_(first_index) {}
  TaskConfig next() { return sampleTask(prior_, seed_, next_++); }
  TaskConfig at(int64_t index) const { return sampleTask(prior_, seed_, index); }

 private:
  const PriorConfig &prior_;
  uint64_t seed_;
  int64_t next_;
};

std::string taskToYAML(const TaskConfig &task);
TaskConfig taskFromYAML(const std::string &yaml);
void saveTask(const TaskConfig &task, const std::string &path);
TaskConfig loadTask(const std::string &path);

} // namespace kg

#endif
