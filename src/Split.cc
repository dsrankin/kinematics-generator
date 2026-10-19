#include "kg/Split.h"

#include <algorithm>
#include <set>

namespace kg {

SplitSpec::By splitByFromString(const std::string &s) {
  if (s == "task") return SplitSpec::By::Task;
  if (s == "topology") return SplitSpec::By::Topology;
  if (s == "process") return SplitSpec::By::Process;
  if (s == "environment") return SplitSpec::By::Environment;
  if (s == "parameter") return SplitSpec::By::Parameter;
  throw std::invalid_argument("unknown split mode '" + s + "'");
}

namespace {

std::string byHash(const std::string &key, const SplitSpec &spec) {
  double tot = 0;
  for (const auto &f : spec.fractions) {
    if (!(f.second >= 0)) throw std::invalid_argument("split fractions must be >= 0");
    tot += f.second;
  }
  if (!(tot > 0)) throw std::invalid_argument("split fractions must have a positive sum");
  const double u = (mixSeed(hashString(key), spec.salt) >> 11) * 0x1.0p-53 * tot;
  double acc = 0;
  for (const auto &f : spec.fractions) {
    acc += f.second;
    if (u < acc) return f.first;
  }
  return spec.fractions.back().first;
}

} // namespace

std::string assignSplit(const TaskConfig &task, const SplitSpec &spec, const Environment *env) {
  switch (spec.by) {
    case SplitSpec::By::Task: return byHash(task.task_id + "#" + std::to_string(task.seed), spec);
    case SplitSpec::By::Topology: return byHash(task.topology(), spec);
    case SplitSpec::By::Process: {
      std::set<std::string> s;
      for (const auto &c : task.components) s.insert(c.process);
      std::string key;
      for (const auto &p : s) key += p + ";";
      return byHash(key, spec);
    }
    case SplitSpec::By::Environment: return assignSplit(env ? *env : task.environment, spec);
    case SplitSpec::By::Parameter: {
      if (spec.labels.size() != spec.edges.size() + 1)
        throw std::invalid_argument("parameter split needs labels.size() == edges.size() + 1");
      if (!std::is_sorted(spec.edges.begin(), spec.edges.end())) throw std::invalid_argument("parameter split edges must be sorted");
      const double x = task.param(spec.parameter);
      const size_t bin = std::upper_bound(spec.edges.begin(), spec.edges.end(), x) - spec.edges.begin();
      return spec.labels[bin];
    }
  }
  return "";
}

std::string assignSplit(const Environment &env, const SplitSpec &spec) { return byHash(env.env_id + "#" + std::to_string(env.seed), spec); }

} // namespace kg
