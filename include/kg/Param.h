// Named, possibly random and possibly dependent numeric parameters.
//
// Every numeric quantity in a prior (masses, widths, branching fractions,
// beam energy, detector resolutions, ...) is a ParamDef in a ParamSet. A
// ParamDef is one of
//   * Fixed:  a number,
//   * Dist:   a draw from a Distribution, optionally multiplied by another
//             parameter (`ratio_of`, used for daughter/parent mass ratios),
//   * Expr:   an arithmetic expression of other parameters.
// and has a scope: `task` (sampled once per task/environment) or `event`
// (re-sampled for every event). Anything depending on an event-scoped
// parameter is itself event-scoped.
#ifndef KG_PARAM_H
#define KG_PARAM_H

#include <string>
#include <unordered_map>
#include <vector>

#include "kg/Distribution.h"
#include "kg/Expression.h"

namespace kg {

enum class Scope { Task, Event };

struct ParamDef {
  enum class Kind { Fixed, Dist, Expr };
  std::string name;
  Kind kind = Kind::Fixed;
  double value = 0;
  Distribution dist;
  std::string ratio_of; // Dist only: multiply sample by this parameter
  Expression expr;
  Scope scope = Scope::Task;
};

class ParamSet {
 public:
  // Add a definition; names must be unique. Returns its index.
  int add(ParamDef def);
  int addFixed(const std::string &name, double value);

  // Parse a YAML value (number, expression string, or map with dist/expr/value)
  // and add it under `name`.
  int parse(const YAML::Node &node, const std::string &name);
  // As parse(), but use `def` when the node is missing.
  int parseOr(const YAML::Node &node, const std::string &name, double def);

  int find(const std::string &name) const; // -1 if absent
  int at(const std::string &name) const;   // throws if absent
  bool has(const std::string &name) const { return find(name) >= 0; }
  size_t size() const { return defs_.size(); }
  const ParamDef &def(int i) const { return defs_.at(i); }
  const std::vector<ParamDef> &defs() const { return defs_; }

  // Resolve references, order dependencies and propagate event scope.
  // Throws on unknown references or dependency cycles.
  void compile();
  bool compiled() const { return compiled_; }

  // Effective scope after propagation.
  bool isEvent(int i) const { return isEvent_.at(i); }
  bool hasEventParams() const { return !eventOrder_.empty(); }

  // Draw all parameters (task-level and an initial draw of event-level ones).
  std::vector<double> sample(Rng &rng) const;
  // Re-draw event-scoped parameters in place.
  void resampleEvent(std::vector<double> &values, Rng &rng) const;

  // Fresh value of parameter i given the other values (a new draw for distributions;
  // used for per-object draws such as individual ISR jets).
  double redraw(int i, const std::vector<double> &values, Rng &rng) const { return evalOne(i, values, rng); }

  // Copy in which every task-scoped parameter is replaced by its value in `values`.
  ParamSet frozen(const std::vector<double> &values) const;

  // name -> definition map (frozen task values emit as plain numbers).
  void emit(YAML::Emitter &out) const;
  static ParamSet fromYAML(const YAML::Node &node);

 private:
  double evalOne(int i, const std::vector<double> &v, Rng &rng) const;
  std::vector<ParamDef> defs_;
  std::unordered_map<std::string, int> index_;
  std::vector<int> refIdx_;
  std::vector<int> order_;
  std::vector<int> eventOrder_;
  std::vector<bool> isEvent_;
  bool compiled_ = false;
};

void emitParamDef(YAML::Emitter &out, const ParamDef &d);

} // namespace kg

#endif
