// Truth generation: recursive ancestry, conservation, masses, multiplicities.
#include <map>
#include <set>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

namespace {

const char *kCascade = R"(
particles:
  X: {pdg: 1000001, mass: 1000, width: 10, decays: [{daughters: [Y, q], br: 1}]}
  Y: {pdg: 1000002, mass: 600, decays: [{daughters: [Z, q], br: 0.5}, {daughters: [N, q, q], br: 0.5}]}
  Z: {pdg: 23, mass: 91.19, width: 2.5, decays: [[l, l]]}
  N: {pdg: 1000022, mass: 100, visible: false}
  q: {pdg: 1}
  l: {pdg: 11}
production: [X, X]
initial_state:
  sqrt_s: {model: exponential_pdf}
  recoil: {model: exponential, min_pt: 10, mean_pt: 40}
)";

void checkEvent(const TruthEvent &ev) {
  // Energy-momentum conservation: incoming partons = final state + recoil.
  const P4 in = ev.incomingSum(), out = ev.finalSum();
  const double tol = 1e-8 * in.e;
  CHECK_CLOSE(in.px, out.px, tol);
  CHECK_CLOSE(in.py, out.py, tol);
  CHECK_CLOSE(in.pz, out.pz, tol);
  CHECK_CLOSE(in.e, out.e, tol);
  for (const auto &p : ev.particles) {
    CHECK(p.index == &p - &ev.particles[0]);
    if (p.species < 0) continue;
    // Generated mass is the four-vector mass.
    CHECK_CLOSE(p.p4.m(), p.mass, 1e-6 * std::max(1.0, p.p4.e));
    if (p.children.empty()) {
      CHECK(p.status == TruthStatus::Final);
      continue;
    }
    CHECK(p.status == TruthStatus::Decayed);
    P4 sum;
    double msum = 0;
    for (int c : p.children) {
      CHECK(ev.particles[c].parent == p.index);
      CHECK(ev.particles[c].depth == p.depth + 1);
      CHECK(ev.particles[c].branch == p.branch);
      sum += ev.particles[c].p4;
      msum += ev.particles[c].mass;
    }
    CHECK_CLOSE(sum.e, p.p4.e, 1e-8 * p.p4.e);
    CHECK_CLOSE(sum.pz, p.p4.pz, 1e-8 * p.p4.e);
    CHECK(msum <= p.mass * (1 + 1e-12)); // decay thresholds respected
  }
}

} // namespace

TEST(recursive_ancestry) {
  const TaskConfig t = sampleTask(kgtest::prior(kCascade), 1, 0);
  std::set<size_t> nfinal;
  for (int i = 0; i < 2000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    checkEvent(ev);
    const auto produced = ev.producedParticles();
    CHECK(produced.size() == 2);
    for (int f : ev.finalState()) {
      const auto anc = ev.ancestors(f);
      CHECK(!anc.empty());
      CHECK(anc.back() == produced[ev.particles[f].branch]); // ancestry ends at a produced particle
      CHECK(int(anc.size()) == ev.particles[f].depth);
      CHECK(ev.isDescendant(f, anc.back()));
      if (ev.particles[f].name == "l") { // l <- Z <- Y <- X
        CHECK(anc.size() == 3);
        CHECK(ev.particles[anc[0]].name == "Z" && ev.particles[anc[1]].name == "Y" && ev.particles[anc[2]].name == "X");
      }
    }
    nfinal.insert(ev.finalState().size());
  }
  // Y -> Z q (-> l l q) or Y -> N q q: per branch 4 final particles either way, + recoil (not final).
  CHECK(nfinal.size() == 1 && *nfinal.begin() == 8);
}

TEST(conservation_for_example_priors) {
  for (const char *f : {"fixed_process.yaml", "randomized_tasks.yaml", "detector_randomization.yaml", "collide_tuned_sm.yaml"}) {
    const PriorConfig p = kg::loadPrior(kgtest::examplePath(f));
    for (int it = 0; it < 20; ++it) {
      const TaskConfig t = sampleTask(p, 77, it);
      for (int i = 0; i < 200; ++i) checkEvent(generateTruthEvent(t, i));
    }
  }
}

TEST(branching_fraction_frequencies) {
  const TaskConfig t = sampleTask(kg::loadPrior(kgtest::examplePath("fixed_process.yaml")), 1, 0);
  const int wp = t.species("W+");
  std::map<int, int> modes;
  int n = 0;
  for (int i = 0; i < 20000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    for (const auto &p : ev.particles)
      if (p.species == wp) {
        modes[p.decay_mode]++;
        ++n;
      }
  }
  CHECK(n == 20000);
  CHECK_CLOSE(modes[2] / double(n), 0.784, 4 * std::sqrt(0.784 * 0.216 / n));
  CHECK_CLOSE(modes[0] / double(n), 0.108, 4 * std::sqrt(0.108 * 0.892 / n));
}

TEST(breit_wigner_masses_and_cutoff) {
  const PriorConfig p = kgtest::prior(R"(
bw_cutoff: 3
particles:
  Z: {pdg: 23, mass: 91.19, width: 2.5, decays: [[l, l]]}
  l: {pdg: 11}
production: [Z, Z]
)");
  const TaskConfig t = sampleTask(p, 1, 0);
  Stats m;
  int inside = 0;
  for (int i = 0; i < 20000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    for (int k : ev.producedParticles()) {
      const double x = ev.particles[k].mass;
      CHECK(x >= 91.19 - 7.5 - 1e-9 && x <= 91.19 + 7.5 + 1e-9);
      inside += std::fabs(x - 91.19) < 1.25;
      m.add(x);
    }
  }
  CHECK_CLOSE(m.mean(), 91.19, 0.1);
  CHECK(inside > 0.55 * m.n && inside < 0.75 * m.n); // ~0.5/(2/pi*atan(3*2)) = 0.63
}

TEST(variable_final_state_multiplicity) {
  const PriorConfig p = kgtest::prior(R"(
particles:
  X: {pdg: 1, mass: 300, object: jet, decays: [{daughters: [q, q], br: 0.5}, {daughters: [q, q, q], br: 0.5}]}
  q: {pdg: 1}
processes:
  - name: multi
    production: [X]
    multiplicity: {values: [1, 2, 3]}
)");
  std::set<size_t> produced, final;
  for (int it = 0; it < 30; ++it) {
    const TaskConfig t = sampleTask(p, 4, it);
    produced.insert(t.components[0].subprocesses[0].produced.size());
    for (int i = 0; i < 50; ++i) {
      const TruthEvent ev = generateTruthEvent(t, i);
      checkEvent(ev);
      final.insert(ev.finalState().size());
    }
  }
  CHECK(produced == (std::set<size_t>{1, 2, 3}));
  for (size_t n = 2; n <= 9; ++n) CHECK_MSG(final.count(n), "missing final-state multiplicity " << n);
}

TEST(single_pair_triple_production) {
  for (const char *prod : {"[Zp]", "[Zp, Zp]", "[Zp, Zp, Zp]"}) {
    const PriorConfig p = kgtest::prior(std::string(R"(
particles:
  Zp: {pdg: 32, mass: 500, width: 10, decays: [[l, l]]}
  l: {pdg: 13}
production: )") + prod + R"(
initial_state: {sqrt_s: {model: threshold_exponential, scale: 100, threshold_offset: 0}, rapidity: {dist: gaussian, mean: 0, sigma: 1, scope: event}}
)");
    const TaskConfig t = sampleTask(p, 2, 0);
    const size_t n = t.components[0].subprocesses[0].produced.size();
    for (int i = 0; i < 300; ++i) {
      const TruthEvent ev = generateTruthEvent(t, i);
      checkEvent(ev);
      double msum = 0;
      for (int k : ev.producedParticles()) msum += ev.particles[k].mass;
      if (n == 1) CHECK_CLOSE(ev.sqrt_s_hat, msum, 1e-9);
      else CHECK(ev.sqrt_s_hat >= msum);
      CHECK(ev.event_params.size() == 1 && ev.event_params[0].first == "initial_state.rapidity");
    }
  }
}

TEST(decay_mode_fixed_per_task) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("randomized_tasks.yaml"));
  for (int it = 0; it < 10; ++it) {
    const TaskConfig t = sampleTask(p, 3, it);
    const int y = t.species("Y");
    std::set<int> modes;
    for (int i = 0; i < 100; ++i)
      for (const auto &q : generateTruthEvent(t, i).particles)
        if (q.species == y) modes.insert(q.decay_mode);
    CHECK(modes.size() <= 1);
    if (!modes.empty()) CHECK(*modes.begin() == t.particles[y].fixed_mode);
  }
}

TEST(legacy_luminosity_threshold_and_recoil) {
  const PriorConfig p = kgtest::prior(kCascade);
  const TaskConfig t = sampleTask(p, 1, 0);
  Stats pt;
  for (int i = 0; i < 5000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    double msum = 0;
    for (int k : ev.producedParticles()) msum += ev.particles[k].mass;
    CHECK(ev.sqrt_s_hat >= msum);
    CHECK(ev.x1 > 0 && ev.x1 <= 1.01 && ev.x2 > 0 && ev.x2 <= 1.01);
    for (const auto &q : ev.particles)
      if (q.status == TruthStatus::Recoil) pt.add(q.p4.pt());
  }
  CHECK_CLOSE(pt.mean(), 40, 4 * pt.err());
}

TEST(long_run_with_widths_is_stable) {
  // Breit-Wigner tails must never push the hard system beyond the collider energy.
  const TaskConfig t = sampleTask(kg::loadPrior(kgtest::examplePath("fixed_process.yaml")), 1, 0);
  double mmax = 0;
  for (int i = 0; i < 100000; ++i)
    for (int k : generateTruthEvent(t, i).producedParticles()) mmax = std::max(mmax, generateTruthEvent(t, i).particles[k].mass);
  CHECK(mmax <= 172.5 + 15 * 1.42 + 1e-9);
}
