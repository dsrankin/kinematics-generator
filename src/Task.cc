#include "kg/Task.h"

#include <cmath>
#include <functional>
#include <set>
#include <sstream>

namespace kg {

int TaskConfig::species(const std::string &name) const {
  auto it = particle_index.find(name);
  if (it == particle_index.end()) throw std::invalid_argument("task '" + task_id + "': unknown particle '" + name + "'");
  return it->second;
}

std::map<std::string, double> TaskConfig::parameterValues() const {
  std::map<std::string, double> out;
  for (size_t i = 0; i < values.size(); ++i) out[params.def(static_cast<int>(i)).name] = values[i];
  return out;
}

int TaskConfig::typeCode(const std::string &type) const {
  for (size_t i = 0; i < object_types.size(); ++i)
    if (object_types[i] == type) return static_cast<int>(i);
  return -1;
}

std::string TaskConfig::topology() const {
  std::function<std::string(int)> sig = [&](int s) -> std::string {
    const TaskParticle &p = particles[s];
    if (p.stable) return p.name;
    auto mode = [&](const TaskDecayMode &m) {
      std::string o = "(";
      for (size_t i = 0; i < m.daughters.size(); ++i) o += (i ? "," : "") + sig(m.daughters[i]);
      return o + ")";
    };
    if (p.fixed_mode >= 0) return p.name + ">" + mode(p.modes[p.fixed_mode]);
    std::vector<std::string> alts;
    for (const auto &m : p.modes)
      if (m.br > 0) alts.push_back(mode(m));
    if (alts.size() == 1) return p.name + ">" + alts[0];
    std::string o = p.name + ">{";
    for (size_t i = 0; i < alts.size(); ++i) o += (i ? "|" : "") + alts[i];
    return o + "}";
  };
  std::string out;
  for (size_t c = 0; c < components.size(); ++c) {
    if (c) out += " ; ";
    out += components[c].process + ": ";
    for (size_t i = 0; i < components[c].produced.size(); ++i) out += (i ? " + " : "") + sig(components[c].produced[i]);
  }
  return out;
}

namespace {

void fail(const std::string &msg) { throw InvalidTask(msg); }

std::string num(double x) {
  std::ostringstream os;
  os << x;
  return os.str();
}

AngularModel resolveAngular(const AngularDef &d, const std::vector<double> &v, const std::string &ctx) {
  AngularModel a;
  a.type = d.type;
  a.axis = d.axis;
  for (int c : d.coeffs) a.coeffs.push_back(v[c]);
  try {
    a.prepare();
  } catch (const std::invalid_argument &e) {
    fail(ctx + ": " + e.what());
  }
  return a;
}

// Lowest mass a species can take in an event (worst case over its decay modes).
double minimumMass(const TaskConfig &t, int s, std::vector<double> &memo) {
  if (memo[s] >= 0) return memo[s];
  const TaskParticle &p = t.particles[s];
  const double M = t.mass(s), G = t.width(s);
  const double cut = t.bw_cutoff;
  double lo;
  if (p.stable) {
    lo = G > 0 ? (cut > 0 ? std::max(0.0, M - cut * G) : 0.0) : M;
  } else {
    double thrMax = 0;
    for (size_t im = 0; im < p.modes.size(); ++im) {
      const auto &m = p.modes[im];
      if (m.br <= 0 && static_cast<int>(im) != p.fixed_mode) continue;
      if (p.fixed_mode >= 0 && static_cast<int>(im) != p.fixed_mode) continue;
      double thr = 0;
      for (int d : m.daughters) thr += minimumMass(t, d, memo);
      std::string ds;
      for (int d : m.daughters) ds += " " + t.particles[d].name;
      if (G == 0 && M + 1e-9 * std::max(1.0, M) < thr)
        fail("decay " + p.name + " ->" + ds + " is kinematically forbidden: m(" + p.name + ")=" + num(M) + " < " + num(thr));
      if (G > 0 && cut > 0 && thr > M + cut * G)
        fail("decay " + p.name + " ->" + ds + " lies beyond the Breit-Wigner cutoff: threshold " + num(thr) + " > " +
             num(M + cut * G));
      thrMax = std::max(thrMax, thr);
    }
    lo = G > 0 ? std::max(thrMax, cut > 0 ? M - cut * G : 0.0) : M;
  }
  memo[s] = lo;
  return lo;
}

} // namespace

void validateTask(const TaskConfig &t) {
  std::set<int> reach;
  std::function<void(int)> walk = [&](int s) {
    if (!reach.insert(s).second) return;
    for (const auto &m : t.particles[s].modes)
      for (int d : m.daughters) walk(d);
  };
  for (const auto &c : t.components)
    for (int s : c.produced) walk(s);
  for (int s : reach) {
    const TaskParticle &p = t.particles[s];
    if (!(t.mass(s) >= 0)) fail("particle '" + p.name + "' has negative mass " + num(t.mass(s)));
    if (!(t.width(s) >= 0)) fail("particle '" + p.name + "' has negative width " + num(t.width(s)));
    if (!p.stable) {
      double sum = 0;
      for (const auto &m : p.modes) {
        if (!(m.br >= 0)) fail("particle '" + p.name + "' has a negative branching fraction");
        sum += m.br;
      }
      if (std::fabs(sum - 1) > 1e-6) fail("particle '" + p.name + "': branching fractions sum to " + num(sum));
    }
  }
  std::vector<double> memo(t.particles.size(), -1.0);
  for (int s : reach) minimumMass(t, s, memo);
  const auto &v = t.values;
  const InitialStateDef &is = t.initial;
  for (const auto &c : t.components) {
    if (c.produced.empty() || c.produced.size() > 3) fail("component '" + c.process + "' must produce 1-3 particles");
    double thr = 0;
    for (int s : c.produced) thr += memo[s];
    if (c.produced.size() > 1) thr += v[is.offset];
    switch (is.model) {
      case SqrtSModel::ExponentialPdf:
        if (!(v[is.ebeam] > 0) || !(v[is.tau] > 0)) fail("initial_state: ebeam and tau must be > 0");
        if (thr >= 2 * v[is.ebeam]) fail("process '" + c.process + "': threshold " + num(thr) + " exceeds the collider energy");
        break;
      case SqrtSModel::ThresholdExponential:
        if (!(v[is.scale] > 0)) fail("initial_state: scale must be > 0");
        break;
      case SqrtSModel::Fixed:
        if (c.produced.size() > 1 && v[is.value] < thr)
          fail("process '" + c.process + "': fixed sqrt(s)=" + num(v[is.value]) + " below threshold " + num(thr));
        break;
      case SqrtSModel::Uniform:
        if (c.produced.size() > 1 && v[is.max] < thr)
          fail("process '" + c.process + "': sqrt(s) range max " + num(v[is.max]) + " below threshold " + num(thr));
        break;
    }
    if (!(c.fraction >= 0)) fail("component '" + c.process + "' has a negative fraction");
  }
  if (is.recoil == RecoilModel::Exponential && !(v[is.recoil_mean] > v[is.recoil_min] && v[is.recoil_min] >= 0))
    fail("initial_state.recoil: requires mean_pt > min_pt >= 0");
}

TaskConfig sampleTask(const PriorConfig &prior) { return sampleTask(prior, prior.seed, 0); }

TaskConfig sampleTask(const PriorConfig &prior, uint64_t master_seed, int64_t index) {
  const uint64_t task_seed = mixSeed(master_seed, kStreamTask, static_cast<uint64_t>(index));
  Rng rng(task_seed);
  std::ostringstream id;
  id << "task_" << index << "_" << std::hex << (task_seed & 0xffffffffULL);
  std::string lastError;
  const int maxAttempts = prior.task.reject_invalid ? std::max(1, prior.task.max_attempts) : 1;
  for (int attempt = 1; attempt <= maxAttempts; ++attempt) {
    TaskConfig t;
    t.task_id = id.str();
    t.prior_name = prior.name;
    t.seed = task_seed;
    t.master_seed = master_seed;
    t.task_index = index;
    t.attempts = attempt;
    t.initial = prior.initial;
    t.bw_cutoff = prior.bw_cutoff;
    t.object_types = prior.objectTypes();
    try {
      t.values = prior.params.sample(rng);
      const auto &v = t.values;
      // Particles.
      t.particle_index = prior.particle_index;
      for (const auto &pd : prior.particles) {
        TaskParticle p;
        p.name = pd.name;
        p.pdg = pd.pdg;
        p.charge = pd.charge;
        p.mass = pd.mass;
        p.width = pd.width;
        p.stable = pd.stable;
        p.object = pd.object;
        double sum = 0;
        for (size_t im = 0; im < pd.modes.size(); ++im) {
          const auto &md = pd.modes[im];
          TaskDecayMode m;
          for (const auto &d : md.daughters) m.daughters.push_back(prior.particle_index.at(d));
          m.br = v[md.br];
          if (!(m.br >= 0) || (!pd.normalize_br && m.br > 1))
            fail("particle '" + pd.name + "': branching fraction " + num(m.br) + " outside [0, 1]");
          sum += m.br;
          m.angular = resolveAngular(md.angular, v, "particle '" + pd.name + "' decay " + std::to_string(im));
          m.model = md.model;
          p.modes.push_back(std::move(m));
        }
        if (!p.modes.empty()) {
          if (pd.normalize_br) {
            if (!(sum > 0)) fail("particle '" + pd.name + "': branching fractions sum to zero");
            for (auto &m : p.modes) m.br /= sum;
          } else if (std::fabs(sum - 1) > 1e-6) {
            fail("particle '" + pd.name + "': branching fractions sum to " + num(sum) + " (set normalize_br: true to rescale)");
          }
        }
        for (const auto &m : p.modes) p.br_weights.push_back(m.br);
        t.particles.push_back(std::move(p));
      }
      for (size_t i = 0; i < prior.particles.size(); ++i)
        if (prior.particles[i].mode_per_task && !t.particles[i].modes.empty())
          t.particles[i].fixed_mode = static_cast<int>(rng.categorical(t.particles[i].br_weights));
      // Components.
      const double nraw = prior.task.n_components.sample(rng);
      const int n = static_cast<int>(std::lround(nraw));
      if (n < 1) fail("task.components sampled " + num(nraw) + " < 1");
      if (prior.task.distinct_processes && n > static_cast<int>(prior.processes.size()))
        fail("task.components=" + std::to_string(n) + " exceeds the number of distinct processes");
      std::vector<double> pw;
      for (const auto &pr : prior.processes) pw.push_back(pr.weight);
      for (int c = 0; c < n; ++c) {
        const size_t ip = rng.categorical(pw);
        if (prior.task.distinct_processes) pw[ip] = 0;
        const ProcessDef &pr = prior.processes[ip];
        TaskComponent comp;
        comp.process = pr.name;
        std::vector<const SlotDef *> slots;
        if (pr.multiplicity.empty()) {
          for (const auto &s : pr.slots) slots.push_back(&s);
        } else {
          const int k = pr.multiplicity_weights.empty() ? pr.multiplicity[static_cast<size_t>(rng.uniform() * pr.multiplicity.size()) % pr.multiplicity.size()]
                                                        : pr.multiplicity[rng.categorical(pr.multiplicity_weights)];
          for (int j = 0; j < k; ++j) slots.push_back(&pr.slots[0]);
        }
        for (const SlotDef *s : slots) {
          size_t ic = 0;
          if (s->choices.size() > 1)
            ic = s->weights.empty() ? static_cast<size_t>(rng.uniform() * s->choices.size()) % s->choices.size() : rng.categorical(s->weights);
          comp.produced.push_back(prior.particle_index.at(s->choices[ic]));
        }
        comp.angular = resolveAngular(pr.angular, v, "process '" + pr.name + "'");
        t.components.push_back(std::move(comp));
      }
      // Mixture fractions.
      if (prior.task.fractions == "fixed") {
        double tot = 0;
        for (double f : prior.task.fixed_fractions) tot += f;
        if (!(tot > 0)) fail("task.fractions sum to zero");
        for (int c = 0; c < n; ++c) t.components[c].fraction = prior.task.fixed_fractions[c] / tot;
      } else if (prior.task.fractions == "dirichlet") {
        double tot = 0;
        for (auto &c : t.components) tot += (c.fraction = rng.gamma(prior.task.dirichlet_alpha));
        for (auto &c : t.components) c.fraction /= tot;
      } else {
        for (auto &c : t.components) c.fraction = 1.0 / n;
      }
      t.params = prior.params.frozen(t.values);
      validateTask(t);
      // Feature schema.
      const std::vector<double> fv = prior.features.params.sample(rng);
      auto names = std::make_shared<std::vector<std::string>>();
      for (const auto &fd : prior.features.defs) {
        const double presence = fv[fd.num.at("presence_prob")];
        if (!(presence >= 0 && presence <= 1)) fail("feature '" + fd.name + "': presence_prob outside [0, 1]");
        if (!(rng.uniform() < presence)) continue;
        FeatureSpec f;
        f.name = fd.name;
        f.kind = fd.kind;
        f.applies_to.insert(fd.applies_to.begin(), fd.applies_to.end());
        f.quantity = fd.quantity;
        f.label = fd.label;
        f.fake_policy = fd.fake_policy;
        f.env_param = fd.env_param;
        for (const auto &kv : fd.num)
          if (kv.first != "presence_prob") f.p[kv.first] = fv[kv.second];
        const double miss = f.p["missing_prob"];
        if (!(miss >= 0 && miss <= 1)) fail("feature '" + fd.name + "': missing_prob outside [0, 1]");
        if (f.kind == FeatureKind::LabelCorrelated) {
          std::vector<std::string> cats;
          if (f.label == "component") {
            for (int c = 0; c < n; ++c) cats.push_back(std::to_string(c));
          } else if (f.label == "process") {
            for (const auto &c : t.components) cats.push_back(c.process);
          } else if (f.label == "is_fake") {
            cats = {"0", "1"};
          } else {
            for (const auto &p : t.particles) cats.push_back(p.name);
            cats.push_back("fake");
            if (f.label == "parent") cats.push_back("none");
          }
          for (const auto &c : cats)
            if (!f.shifts.count(c)) f.shifts[c] = fd.shift.sample(rng);
        }
        names->push_back(f.name);
        t.features.push_back(std::move(f));
      }
      t.feature_names = names;
      // Default detector environment for this task.
      try {
        t.environment = sampleEnvironment(prior.detector, mixSeed(task_seed, kStreamEnvironment), t.task_id + "/env");
      } catch (const std::invalid_argument &e) {
        fail(e.what());
      }
      return t;
    } catch (const InvalidTask &e) {
      lastError = e.what();
      if (!prior.task.reject_invalid) throw;
    }
  }
  throw InvalidTask("prior '" + prior.name + "': no valid task after " + std::to_string(maxAttempts) +
                    " attempts; last problem: " + lastError);
}

} // namespace kg
