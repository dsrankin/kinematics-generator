#include "kg/Param.h"

#include <cmath>
#include <functional>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace kg {

namespace {

Scope parseScope(const YAML::Node &n, const std::string &name) {
  if (!n) return Scope::Task;
  const std::string s = n.as<std::string>();
  if (s == "task" || s == "environment" || s == "env") return Scope::Task;
  if (s == "event") return Scope::Event;
  throw std::invalid_argument("parameter '" + name + "': unknown scope '" + s + "'");
}

bool scalarAsDouble(const YAML::Node &n, double &out) {
  try {
    out = n.as<double>();
    return true;
  } catch (const YAML::Exception &) {
    return false;
  }
}

} // namespace

int ParamSet::add(ParamDef def) {
  if (def.name.empty()) throw std::invalid_argument("parameter with empty name");
  if (index_.count(def.name)) throw std::invalid_argument("parameter '" + def.name + "' defined twice");
  if (def.kind == ParamDef::Kind::Dist) def.dist.validate("parameter '" + def.name + "'");
  const int i = static_cast<int>(defs_.size());
  index_[def.name] = i;
  defs_.push_back(std::move(def));
  compiled_ = false;
  return i;
}

int ParamSet::addFixed(const std::string &name, double value) {
  ParamDef d;
  d.name = name;
  d.kind = ParamDef::Kind::Fixed;
  d.value = value;
  return add(d);
}

int ParamSet::parse(const YAML::Node &n, const std::string &name) {
  if (!n || n.IsNull()) throw std::invalid_argument("parameter '" + name + "' is missing");
  ParamDef d;
  d.name = name;
  if (n.IsScalar()) {
    double v;
    if (scalarAsDouble(n, v)) {
      d.kind = ParamDef::Kind::Fixed;
      d.value = v;
    } else {
      d.kind = ParamDef::Kind::Expr;
      d.expr = Expression(n.as<std::string>());
    }
  } else if (n.IsMap()) {
    d.scope = parseScope(n["scope"], name);
    if (n["dist"]) {
      d.kind = ParamDef::Kind::Dist;
      d.dist = Distribution::fromYAML(n, "parameter '" + name + "'", true); // also rejects unknown keys
      if (n["ratio_of"]) d.ratio_of = n["ratio_of"].as<std::string>();
    } else if (n["expr"] || n["value"]) {
      const char *k = n["expr"] ? "expr" : "value";
      for (const auto &kv : n) {
        const std::string key = kv.first.as<std::string>();
        if (key != k && key != "scope") throw std::invalid_argument("parameter '" + name + "': unknown key '" + key + "'");
      }
      if (n["expr"]) {
        d.kind = ParamDef::Kind::Expr;
        d.expr = Expression(n["expr"].as<std::string>());
      } else {
        d.kind = ParamDef::Kind::Fixed;
        d.value = n["value"].as<double>();
      }
    } else {
      throw std::invalid_argument("parameter '" + name + "': map must contain 'dist', 'expr' or 'value'");
    }
  } else {
    throw std::invalid_argument("parameter '" + name + "': expected a number, expression or map");
  }
  return add(std::move(d));
}

int ParamSet::parseOr(const YAML::Node &n, const std::string &name, double def) {
  if (!n || n.IsNull()) return addFixed(name, def);
  return parse(n, name);
}

int ParamSet::find(const std::string &name) const {
  auto it = index_.find(name);
  return it == index_.end() ? -1 : it->second;
}

int ParamSet::at(const std::string &name) const {
  const int i = find(name);
  if (i < 0) throw std::invalid_argument("unknown parameter '" + name + "'");
  return i;
}

void ParamSet::compile() {
  const size_t n = defs_.size();
  auto resolve = [&](const std::string &ref) {
    int i = find(ref);
    if (i < 0) i = find(ref + ".mass"); // bare particle name means its mass
    return i;
  };
  std::vector<std::vector<int>> deps(n);
  refIdx_.assign(n, -1);
  for (size_t i = 0; i < n; ++i) {
    ParamDef &d = defs_[i];
    if (d.kind == ParamDef::Kind::Dist && !d.ratio_of.empty()) {
      refIdx_[i] = resolve(d.ratio_of);
      if (refIdx_[i] < 0) throw std::invalid_argument("parameter '" + d.name + "': ratio_of unknown parameter '" + d.ratio_of + "'");
      deps[i].push_back(refIdx_[i]);
    } else if (d.kind == ParamDef::Kind::Expr) {
      d.expr.bind([&](const std::string &v) {
        const int j = resolve(v);
        if (j >= 0) deps[i].push_back(j);
        return j;
      });
    }
  }
  // Topological order (DFS) with cycle detection.
  std::vector<int> state(n, 0);
  order_.clear();
  std::function<void(int, std::vector<int> &)> visit = [&](int i, std::vector<int> &stack) {
    if (state[i] == 2) return;
    stack.push_back(i);
    if (state[i] == 1) {
      std::string cyc;
      for (size_t k = 0; k < stack.size(); ++k) cyc += (k ? " -> " : "") + defs_[stack[k]].name;
      throw std::invalid_argument("dependency cycle among parameters: " + cyc);
    }
    state[i] = 1;
    for (int j : deps[i]) visit(j, stack);
    state[i] = 2;
    stack.pop_back();
    order_.push_back(i);
  };
  for (size_t i = 0; i < n; ++i) {
    std::vector<int> stack;
    visit(static_cast<int>(i), stack);
  }
  isEvent_.assign(n, false);
  eventOrder_.clear();
  for (int i : order_) {
    bool ev = defs_[i].scope == Scope::Event && defs_[i].kind != ParamDef::Kind::Fixed;
    for (int j : deps[i]) ev = ev || isEvent_[j];
    isEvent_[i] = ev;
    if (ev) eventOrder_.push_back(i);
  }
  compiled_ = true;
}

double ParamSet::evalOne(int i, const std::vector<double> &v, Rng &rng) const {
  const ParamDef &d = defs_[i];
  switch (d.kind) {
    case ParamDef::Kind::Fixed: return d.value;
    case ParamDef::Kind::Dist: {
      double x = d.dist.sample(rng);
      if (refIdx_[i] >= 0) x *= v[refIdx_[i]];
      return x;
    }
    case ParamDef::Kind::Expr: return d.expr.eval(v);
  }
  return 0;
}

std::vector<double> ParamSet::sample(Rng &rng) const {
  if (!compiled_) throw std::logic_error("ParamSet::sample called before compile()");
  std::vector<double> v(defs_.size(), 0.0);
  for (int i : order_) v[i] = evalOne(i, v, rng);
  return v;
}

void ParamSet::resampleEvent(std::vector<double> &v, Rng &rng) const {
  for (int i : eventOrder_) v[i] = evalOne(i, v, rng);
}

ParamSet ParamSet::frozen(const std::vector<double> &values) const {
  ParamSet out;
  for (size_t i = 0; i < defs_.size(); ++i) {
    ParamDef d = defs_[i];
    if (!isEvent_[i]) {
      d.kind = ParamDef::Kind::Fixed;
      d.value = values[i];
      d.ratio_of.clear();
      d.expr = Expression();
      d.scope = Scope::Task;
    }
    out.add(std::move(d));
  }
  out.compile();
  return out;
}

void emitParamDef(YAML::Emitter &out, const ParamDef &d) {
  switch (d.kind) {
    case ParamDef::Kind::Fixed: out << d.value; break;
    case ParamDef::Kind::Dist:
      out << YAML::Flow << YAML::BeginMap;
      d.dist.emitFields(out);
      if (!d.ratio_of.empty()) out << YAML::Key << "ratio_of" << YAML::Value << d.ratio_of;
      if (d.scope == Scope::Event) out << YAML::Key << "scope" << YAML::Value << "event";
      out << YAML::EndMap;
      break;
    case ParamDef::Kind::Expr:
      out << YAML::Flow << YAML::BeginMap << YAML::Key << "expr" << YAML::Value << d.expr.text();
      if (d.scope == Scope::Event) out << YAML::Key << "scope" << YAML::Value << "event";
      out << YAML::EndMap;
      break;
  }
}

void ParamSet::emit(YAML::Emitter &out) const {
  out << YAML::BeginMap;
  for (const auto &d : defs_) {
    out << YAML::Key << d.name << YAML::Value;
    emitParamDef(out, d);
  }
  out << YAML::EndMap;
}

ParamSet ParamSet::fromYAML(const YAML::Node &node) {
  ParamSet out;
  if (node) {
    if (!node.IsMap()) throw std::invalid_argument("'parameters' must be a map");
    for (const auto &kv : node) out.parse(kv.second, kv.first.as<std::string>());
  }
  out.compile();
  return out;
}

} // namespace kg
