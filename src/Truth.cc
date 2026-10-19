#include <algorithm>
#include <cmath>
#include <functional>
#include <sstream>

#include "kg/Generator.h"

namespace kg {

// ---------------------------------------------------------------------------
// TruthEvent helpers

std::vector<int> TruthEvent::finalState() const {
  std::vector<int> out;
  for (const auto &p : particles)
    if (p.isFinal()) out.push_back(p.index);
  return out;
}

std::vector<int> TruthEvent::producedParticles() const {
  std::vector<int> out;
  for (const auto &p : particles)
    if (p.species >= 0 && p.parent < 0) out.push_back(p.index);
  return out;
}

std::vector<int> TruthEvent::ancestors(int i) const {
  std::vector<int> out;
  for (int p = particles.at(i).parent; p >= 0; p = particles[p].parent) out.push_back(p);
  return out;
}

bool TruthEvent::isDescendant(int i, int ancestor) const {
  for (int p = particles.at(i).parent; p >= 0; p = particles[p].parent)
    if (p == ancestor) return true;
  return false;
}

P4 TruthEvent::incomingSum() const {
  P4 s;
  for (const auto &p : particles)
    if (p.status == TruthStatus::Incoming) s += p.p4;
  return s;
}

P4 TruthEvent::finalSum() const {
  P4 s;
  for (const auto &p : particles)
    if (p.status == TruthStatus::Final || p.status == TruthStatus::Recoil) s += p.p4;
  return s;
}

std::pair<double, double> TruthEvent::trueMET() const {
  double x = 0, y = 0;
  for (const auto &p : particles)
    if (p.isFinal() && p.object_type == "invisible") {
      x += p.p4.px;
      y += p.p4.py;
    }
  return {x, y};
}

std::string TruthEvent::topology() const {
  std::function<std::string(int)> sig = [&](int i) -> std::string {
    const TruthParticle &p = particles[i];
    if (p.children.empty()) return p.name;
    std::string o = p.name + ">(";
    for (size_t k = 0; k < p.children.size(); ++k) o += (k ? "," : "") + sig(p.children[k]);
    return o + ")";
  };
  std::string out;
  for (int i : producedParticles()) out += (out.empty() ? "" : " + ") + sig(i);
  return out;
}

// ---------------------------------------------------------------------------
// Generation

namespace {

struct Node {
  int species;
  int parent;
  int mode = -1;
  std::vector<int> children;
  double mass = 0, minMass = 0;
  int branch = 0, depth = 0;
  P4 p4;
};

std::string num(double x) {
  std::ostringstream os;
  os << x;
  return os.str();
}

double bwUpper(const TaskConfig &t, int s, const std::vector<double> &v) {
  const double M = v[t.particles[s].mass], G = v[t.particles[s].width];
  return (G > 0 && t.bw_cutoff > 0) ? M + t.bw_cutoff * G : (G > 0 ? std::numeric_limits<double>::infinity() : M);
}

double sampleMass(const TaskConfig &t, int s, const std::vector<double> &v, double lo, double hi, Rng &rng) {
  const double M = v[t.particles[s].mass], G = v[t.particles[s].width];
  if (G <= 0) return M;
  return sampleRelativisticBW(rng, M, G, lo, hi);
}

} // namespace

TruthEvent generateTruthEvent(TaskConfig &task) { return generateTruthEvent(task, task.next_event++); }

TruthEvent generateTruthEvent(const TaskConfig &t, uint64_t event_index) {
  TruthEvent ev;
  ev.task_id = t.task_id;
  ev.task_seed = t.seed;
  ev.event_index = event_index;
  ev.seed = mixSeed(t.seed, kStreamEvent, event_index);
  Rng rng(ev.seed);
  const std::string ctx = "task '" + t.task_id + "' event " + std::to_string(event_index);

  // Event-scoped parameters.
  std::vector<double> v = t.values;
  if (t.params.hasEventParams()) {
    t.params.resampleEvent(v, rng);
    for (size_t i = 0; i < v.size(); ++i)
      if (t.params.isEvent(static_cast<int>(i))) ev.event_params.emplace_back(t.params.def(static_cast<int>(i)).name, v[i]);
  }

  // Component (process) of this event.
  if (t.components.size() > 1) {
    std::vector<double> w;
    for (const auto &c : t.components) w.push_back(c.fraction);
    ev.component = static_cast<int>(rng.categorical(w));
  }
  const TaskComponent &comp = t.components[ev.component];
  ev.process = comp.process;

  // 1) Decay tree structure: choose decay modes top-down.
  std::vector<Node> nodes;
  for (size_t b = 0; b < comp.produced.size(); ++b) {
    Node n;
    n.species = comp.produced[b];
    n.parent = -1;
    n.branch = static_cast<int>(b);
    nodes.push_back(n);
  }
  for (size_t i = 0; i < nodes.size(); ++i) {
    const TaskParticle &p = t.particles[nodes[i].species];
    if (p.stable) continue;
    const int mode = p.fixed_mode >= 0 ? p.fixed_mode : static_cast<int>(rng.categorical(p.br_weights));
    nodes[i].mode = mode;
    for (int d : p.modes[mode].daughters) {
      Node c;
      c.species = d;
      c.parent = static_cast<int>(i);
      c.branch = nodes[i].branch;
      c.depth = nodes[i].depth + 1;
      nodes[i].children.push_back(static_cast<int>(nodes.size()));
      nodes.push_back(c);
    }
    if (nodes.size() > 100000) throw std::runtime_error(ctx + ": decay tree exceeds 100000 particles");
  }

  // 2) Mass thresholds bottom-up (children always follow their parent).
  for (size_t k = nodes.size(); k-- > 0;) {
    Node &n = nodes[k];
    const double M = v[t.particles[n.species].mass], G = v[t.particles[n.species].width];
    double thr = 0;
    for (int c : n.children) thr += nodes[c].minMass;
    if (G > 0) {
      n.minMass = std::max(thr, t.bw_cutoff > 0 ? M - t.bw_cutoff * G : 0.0);
      if (n.minMass > bwUpper(t, n.species, v))
        throw std::runtime_error(ctx + ": " + t.particles[n.species].name + " threshold " + num(thr) + " beyond its Breit-Wigner cutoff");
    } else {
      if (M + 1e-9 * std::max(1.0, M) < thr)
        throw std::runtime_error(ctx + ": decay of " + t.particles[n.species].name + " (m=" + num(M) + ") below threshold " + num(thr));
      n.minMass = M;
    }
  }

  // 3) Masses top-down.
  // Produced particles: at most the collider energy minus the other produced particles' minima.
  const InitialStateDef &is = t.initial;
  double rootMin = 0;
  for (auto &n : nodes)
    if (n.parent < 0) rootMin += n.minMass;
  const double sMax = is.model == SqrtSModel::ExponentialPdf ? 2 * v[is.ebeam] - (comp.produced.size() > 1 ? v[is.offset] : 0.0)
                    : is.model == SqrtSModel::Fixed && comp.produced.size() > 1 ? v[is.value] - v[is.offset]
                    : is.model == SqrtSModel::Uniform && comp.produced.size() > 1 ? v[is.max] - v[is.offset]
                    : std::numeric_limits<double>::infinity();
  for (auto &n : nodes)
    if (n.parent < 0) {
      const double hi = std::min(bwUpper(t, n.species, v), sMax - (rootMin - n.minMass));
      if (hi < n.minMass) throw std::runtime_error(ctx + ": produced particles exceed the available energy");
      n.mass = sampleMass(t, n.species, v, n.minMass, hi, rng);
    }
  for (auto &n : nodes) {
    if (n.children.empty()) continue;
    double sumMin = 0;
    for (int c : n.children) sumMin += nodes[c].minMass;
    bool ok = false;
    for (int attempt = 0; attempt < 10000 && !ok; ++attempt) {
      double sum = 0;
      for (int c : n.children) {
        Node &ch = nodes[c];
        const double hi = std::min(bwUpper(t, ch.species, v), n.mass - (sumMin - ch.minMass));
        if (hi < ch.minMass - 1e-9)
          throw std::runtime_error(ctx + ": no phase space for " + t.particles[ch.species].name + " in decay of " +
                                   t.particles[n.species].name);
        ch.mass = sampleMass(t, ch.species, v, ch.minMass, std::max(hi, ch.minMass), rng);
        sum += ch.mass;
      }
      ok = sum <= n.mass * (1 + 1e-12);
    }
    if (!ok) throw std::runtime_error(ctx + ": could not sample daughter masses of " + t.particles[n.species].name);
  }

  // 4) Initial state: hard-system mass, rapidity and recoil.
  double msum = 0;
  for (auto &n : nodes)
    if (n.parent < 0) msum += n.mass;
  const bool single = comp.produced.size() == 1;
  const double ebeam = v[is.ebeam];
  double sqrtS = 0, y = 0;
  if (single) {
    sqrtS = msum;
    if (is.rapidity_from_pdf) {
      // x1 x2 = s/(4 Eb^2) fixed; density in y ~ exp(-2 sqrt(x1x2) cosh(y)/tau).
      const double tau = v[is.tau] + v[is.tau_slope] * (sqrtS / 1000.0);
      const double sq = sqrtS / (2 * ebeam);
      if (!(sq < 1)) throw std::runtime_error(ctx + ": resonance mass exceeds collider energy");
      const double ymax = -std::log(sq);
      bool ok = false;
      for (int a = 0; a < is.max_attempts && !ok; ++a) {
        y = rng.uniform(-ymax, ymax);
        ok = rng.uniform() < std::exp(-2 * sq * (std::cosh(y) - 1) / tau);
      }
      if (!ok) throw std::runtime_error(ctx + ": rapidity sampling failed");
    } else {
      y = v[is.rapidity];
    }
  } else {
    const double thr = msum + v[is.offset];
    switch (is.model) {
      case SqrtSModel::ExponentialPdf: {
        // As myGenerator::getRandCM: x1, x2 ~ Exp(tau) on (0,1], accept sqrt(s) >= threshold.
        const double tau = v[is.tau] + v[is.tau_slope] * (thr / 1000.0);
        bool ok = false;
        double x1 = 0, x2 = 0;
        for (int a = 0; a < is.max_attempts && !ok; ++a) {
          x1 = rng.exponential(tau);
          if (x1 > 1) continue;
          x2 = rng.exponential(tau);
          if (x2 > 1) continue;
          sqrtS = 2 * ebeam * std::sqrt(x1 * x2);
          ok = sqrtS >= thr;
        }
        if (!ok) throw std::runtime_error(ctx + ": sqrt(s) sampling failed above threshold " + num(thr));
        y = is.rapidity_from_pdf ? 0.5 * std::log(x1 / x2) : v[is.rapidity];
        break;
      }
      case SqrtSModel::ThresholdExponential:
        sqrtS = thr + rng.exponential(v[is.scale]);
        y = v[is.rapidity];
        break;
      case SqrtSModel::Fixed:
        sqrtS = v[is.value];
        if (sqrtS < thr) throw std::runtime_error(ctx + ": fixed sqrt(s) below threshold " + num(thr));
        y = v[is.rapidity];
        break;
      case SqrtSModel::Uniform: {
        const double lo = std::max(v[is.min], thr);
        if (v[is.max] < lo) throw std::runtime_error(ctx + ": sqrt(s) range below threshold " + num(thr));
        sqrtS = rng.uniform(lo, v[is.max]);
        y = v[is.rapidity];
        break;
      }
    }
  }
  double recoilPt = 0;
  if (is.recoil == RecoilModel::Exponential) recoilPt = v[is.recoil_min] + rng.exponential(v[is.recoil_mean] - v[is.recoil_min]);
  else if (is.recoil == RecoilModel::Param) recoilPt = v[is.recoil_pt];
  if (recoilPt < 0) throw std::runtime_error(ctx + ": negative recoil pT");
  const double phi = rng.uniform(0, 2 * M_PI);
  const P4 system = P4::fromMassRapidityPt(sqrtS, y, recoilPt, phi);
  P4 recoil;
  if (is.recoil != RecoilModel::None) {
    const double yr = v[is.recoil_rapidity];
    recoil = P4(-recoilPt * std::cos(phi), -recoilPt * std::sin(phi), recoilPt * std::sinh(yr), recoilPt * std::cosh(yr));
  }
  const P4 tot = system + recoil;
  const P4 in1(0, 0, 0.5 * (tot.e + tot.pz), 0.5 * (tot.e + tot.pz));
  const P4 in2(0, 0, -0.5 * (tot.e - tot.pz), 0.5 * (tot.e - tot.pz));
  ev.sqrt_s_hat = sqrtS;
  ev.x1 = ebeam > 0 ? in1.e / ebeam : 0;
  ev.x2 = ebeam > 0 ? in2.e / ebeam : 0;

  // 5) Production and decays.
  std::vector<int> roots;
  for (size_t i = 0; i < nodes.size(); ++i)
    if (nodes[i].parent < 0) roots.push_back(static_cast<int>(i));
  if (roots.size() == 1) {
    nodes[roots[0]].p4 = system;
  } else if (roots.size() == 2) {
    const double m[2] = {nodes[roots[0]].mass, nodes[roots[1]].mass};
    P4 out[2];
    decayTwoBody(system, m, comp.angular, rng, out);
    nodes[roots[0]].p4 = out[0];
    nodes[roots[1]].p4 = out[1];
  } else {
    const double m[3] = {nodes[roots[0]].mass, nodes[roots[1]].mass, nodes[roots[2]].mass};
    P4 out[3];
    decayThreeBody(system, m, ThreeBodyModel::PhaseSpace, comp.angular, rng, out);
    for (int k = 0; k < 3; ++k) nodes[roots[k]].p4 = out[k];
  }
  for (auto &n : nodes) {
    if (n.children.empty()) continue;
    const TaskDecayMode &mode = t.particles[n.species].modes[n.mode];
    // Use the sampled mass exactly (avoid round-off from p4.m()).
    P4 parent = P4::fromVectM(n.p4.vect(), n.mass);
    if (n.children.size() == 2) {
      const double m[2] = {nodes[n.children[0]].mass, nodes[n.children[1]].mass};
      P4 out[2];
      decayTwoBody(parent, m, mode.angular, rng, out);
      for (int k = 0; k < 2; ++k) nodes[n.children[k]].p4 = out[k];
    } else {
      const double m[3] = {nodes[n.children[0]].mass, nodes[n.children[1]].mass, nodes[n.children[2]].mass};
      P4 out[3];
      decayThreeBody(parent, m, mode.model, mode.angular, rng, out);
      for (int k = 0; k < 3; ++k) nodes[n.children[k]].p4 = out[k];
    }
  }

  // 6) Record: incoming partons, recoil, then the decay tree.
  auto parton = [&](const P4 &p, int idx) {
    TruthParticle tp;
    tp.index = idx;
    tp.name = "parton";
    tp.pdg = 21;
    tp.p4 = p;
    tp.status = TruthStatus::Incoming;
    tp.object_type = "invisible";
    return tp;
  };
  ev.particles.reserve(nodes.size() + 3);
  ev.particles.push_back(parton(in1, 0));
  ev.particles.push_back(parton(in2, 1));
  if (is.recoil != RecoilModel::None) {
    TruthParticle r;
    r.index = 2;
    r.name = "recoil";
    r.pdg = 21;
    r.p4 = recoil;
    r.status = is.emit_recoil ? TruthStatus::Final : TruthStatus::Recoil;
    r.object_type = is.emit_recoil ? is.recoil_object : "invisible";
    ev.particles.push_back(r);
  }
  const int offset = static_cast<int>(ev.particles.size());
  for (size_t i = 0; i < nodes.size(); ++i) {
    const Node &n = nodes[i];
    const TaskParticle &sp = t.particles[n.species];
    TruthParticle tp;
    tp.index = offset + static_cast<int>(i);
    tp.name = sp.name;
    tp.pdg = sp.pdg;
    tp.charge = sp.charge;
    tp.p4 = n.p4;
    tp.mass = n.mass;
    tp.parent = n.parent < 0 ? -1 : offset + n.parent;
    for (int c : n.children) tp.children.push_back(offset + c);
    tp.status = sp.stable ? TruthStatus::Final : TruthStatus::Decayed;
    tp.species = n.species;
    tp.branch = n.branch;
    tp.depth = n.depth;
    tp.decay_mode = n.mode;
    tp.object_type = sp.stable ? sp.object : "invisible";
    ev.particles.push_back(std::move(tp));
  }
  return ev;
}

GeneratedEvent generateEvent(const TaskConfig &task, const Environment &env, uint64_t event_index) {
  GeneratedEvent g;
  g.truth = generateTruthEvent(task, event_index);
  g.observed = DetectorModel(task).apply(g.truth, env);
  return g;
}

} // namespace kg
