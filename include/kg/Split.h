// Deterministic dataset-split assignment for tasks and environments.
//
// Hash-based modes give every task with the same key (task id, topology,
// process set or environment id) the same split, independent of generation
// order. Parameter mode bins a task-level parameter (e.g. a mass) so that
// held-out mass ranges never appear in training.
#ifndef KG_SPLIT_H
#define KG_SPLIT_H

#include <string>
#include <utility>
#include <vector>

#include "kg/Task.h"

namespace kg {

struct SplitSpec {
  enum class By { Task, Topology, Process, Environment, Parameter };
  By by = By::Task;
  std::vector<std::pair<std::string, double>> fractions = {{"train", 0.8}, {"val", 0.1}, {"test", 0.1}};
  uint64_t salt = 0;
  // Parameter mode: labels[i] for edges[i-1] <= value < edges[i].
  std::string parameter;
  std::vector<double> edges;
  std::vector<std::string> labels;
};

SplitSpec::By splitByFromString(const std::string &s);

// Split label of `task` (with `env` for By::Environment; defaults to the task's own environment).
std::string assignSplit(const TaskConfig &task, const SplitSpec &spec, const Environment *env = nullptr);
std::string assignSplit(const Environment &env, const SplitSpec &spec);

} // namespace kg

#endif
