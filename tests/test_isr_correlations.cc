// Multi-jet ISR recoil and correlated detector response.

#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

namespace {

PriorConfig zmumu(const std::string &recoil, const std::string &detector = "") {
  return kgtest::prior(R"(
particles:
  Z: {pdg: 23, mass: 91.19, decays: [[mu+, mu-]]}
  mu+: {pdg: -13, charge: 1}
  mu-: {pdg: 13, charge: -1}
production: [Z]
initial_state:
  sqrt_s: {model: exponential_pdf}
  recoil: )" + recoil + "\n" + detector);
}

std::vector<int> recoilJets(const TruthEvent &ev) {
  std::vector<int> out;
  for (const auto &p : ev.particles)
    if (p.name == "recoil") out.push_back(p.index);
  return out;
}

double pearson(const std::vector<double> &a, const std::vector<double> &b) {
  Stats sa, sb;
  double sab = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    sa.add(a[i]);
    sb.add(b[i]);
    sab += a[i] * b[i];
  }
  return (sab / a.size() - sa.mean() * sb.mean()) / (sa.sd() * sb.sd());
}

// Pair of stable jets with a perfect detector apart from `jet`.
PriorConfig dijet(const std::string &jet, const std::string &corr) {
  return kgtest::prior(R"(
particles:
  q: {pdg: 1}
production: [q, q]
initial_state: {sqrt_s: {model: fixed, value: 500}, rapidity: 0}
detector:
  objects:
    jet: )" + jet + "\n  correlations:\n    jes: " + corr + "\n");
}

} // namespace

TEST(isr_fixed_jet_count_conserves_momentum) {
  const PriorConfig p = zmumu("{model: exponential, min_pt: 20, mean_pt: 40, n_jets: 3, emit: true, object: jet}");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  for (int i = 0; i < 2000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const auto jets = recoilJets(ev);
    CHECK(jets.size() == 3);
    P4 isr, sys;
    for (int j : jets) {
      CHECK(ev.particles[j].isFinal() && ev.particles[j].p4.pt() >= 20);
      CHECK_CLOSE(ev.particles[j].p4.m(), 0, 1e-6 * ev.particles[j].p4.e);
      isr += ev.particles[j].p4;
    }
    for (int k : ev.producedParticles()) sys += ev.particles[k].p4;
    CHECK_CLOSE(isr.px + sys.px, 0, 1e-8 * ev.sqrt_s_hat);
    CHECK_CLOSE(isr.py + sys.py, 0, 1e-8 * ev.sqrt_s_hat);
    const P4 in = ev.incomingSum(), out = ev.finalSum();
    CHECK_CLOSE(in.e, out.e, 1e-8 * in.e);
    CHECK_CLOSE(in.pz, out.pz, 1e-8 * in.e);
    const ObservedEvent ob = det.apply(ev);
    int nj = 0;
    for (const auto &o : ob.objects) nj += o.type == "jet";
    CHECK(nj == 3);
    CHECK_CLOSE(ob.met(), 0, 1e-6 * ev.sqrt_s_hat);
  }
}

TEST(isr_poisson_multiplicity_per_event) {
  const PriorConfig p = zmumu("{model: exponential, min_pt: 15, mean_pt: 30, n_jets: {dist: poisson, mean: 1.5, scope: event}}");
  const TaskConfig t = sampleTask(p, 1, 0);
  Stats n;
  int zero = 0;
  const int N = 20000;
  for (int i = 0; i < N; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const size_t k = recoilJets(ev).size();
    n.add(k);
    if (k == 0) {
      ++zero;
      P4 sys;
      for (int j : ev.producedParticles()) sys += ev.particles[j].p4;
      CHECK_CLOSE(sys.pt(), 0, 1e-9 * ev.sqrt_s_hat);
    }
  }
  CHECK_CLOSE(n.mean(), 1.5, 4 * n.err());
  CHECK_CLOSE(double(zero) / N, std::exp(-1.5), 4 * std::sqrt(0.22 * 0.78 / N));
}

TEST(isr_per_jet_parameter_draws) {
  const PriorConfig p = zmumu("{model: param, pt: {dist: exponential, scale: 30, scope: event}, n_jets: 2, "
                              "rapidity: {dist: uniform, min: -3, max: 3, scope: event}}");
  const TaskConfig t = sampleTask(p, 1, 0);
  Stats pt;
  for (int i = 0; i < 10000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const auto j = recoilJets(ev);
    CHECK(j.size() == 2);
    CHECK(ev.particles[j[0]].p4.pt() != ev.particles[j[1]].p4.pt());
    CHECK(ev.particles[j[0]].p4.eta() != ev.particles[j[1]].p4.eta());
    for (int k : j) pt.add(ev.particles[k].p4.pt());
  }
  CHECK_CLOSE(pt.mean(), 30, 4 * pt.err());
}

TEST(isr_single_jet_default_unchanged) {
  const TaskConfig a = sampleTask(zmumu("{model: exponential, min_pt: 10, mean_pt: 30}"), 1, 0);
  const TaskConfig b = sampleTask(zmumu("{model: exponential, min_pt: 10, mean_pt: 30, n_jets: 1}"), 1, 0);
  for (int i = 0; i < 50; ++i) CHECK(truthEventToYAML(generateTruthEvent(a, i)) == truthEventToYAML(generateTruthEvent(b, i)));
  // Event-scoped jet counts survive task serialisation.
  const TaskConfig c = sampleTask(zmumu("{model: exponential, min_pt: 10, mean_pt: 30, n_jets: {dist: poisson, mean: 2, scope: event}}"), 1, 0);
  const TaskConfig r = taskFromYAML(taskToYAML(c));
  for (int i = 0; i < 50; ++i) CHECK(truthEventToYAML(generateTruthEvent(c, i)) == truthEventToYAML(generateTruthEvent(r, i)));
}

TEST(isr_negative_count_rejected) {
  const PriorConfig p = zmumu("{model: exponential, n_jets: -1}", "task: {validation: error}\n");
  CHECK_THROWS_AS(sampleTask(p, 1, 0), InvalidTask, "n_jets");
}

TEST(correlated_scale_shift) {
  const PriorConfig p = dijet("{}", "{types: [jet], scale_sigma: 0.05}");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  Stats shift;
  for (int i = 0; i < 20000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev);
    CHECK(ob.objects.size() == 2 && ob.detector_params.size() == 1);
    CHECK(ob.detector_params[0].first == "detector.correlations.jes.shift");
    const double d = ob.detector_params[0].second;
    for (const auto &o : ob.objects) CHECK_CLOSE(o.pt() / ev.particles[o.truth_index].p4.pt(), 1 + d, 1e-9);
    shift.add(d);
  }
  CHECK_CLOSE(shift.mean(), 0, 4 * shift.err());
  CHECK_CLOSE(shift.sd(), 0.05, 0.001);
}

TEST(correlated_resolution_fluctuations) {
  for (double rho : {0.0, 0.5, 1.0}) {
    const PriorConfig p = dijet("{resolution: {constant: 0.1}}", "{types: [jet], resolution_correlation: " + std::to_string(rho) + "}");
    const TaskConfig t = sampleTask(p, 1, 0);
    const DetectorModel det(t);
    std::vector<double> r1, r2;
    Stats all;
    for (int i = 0; i < 20000; ++i) {
      const TruthEvent ev = generateTruthEvent(t, i);
      const ObservedEvent ob = det.apply(ev);
      // Order by truth index so the pair is not sorted by the smeared pT.
      const ObservedObject *a = &ob.objects[0], *b = &ob.objects[1];
      if (a->truth_index > b->truth_index) std::swap(a, b);
      r1.push_back(a->pt() / ev.particles[a->truth_index].p4.pt() - 1);
      r2.push_back(b->pt() / ev.particles[b->truth_index].p4.pt() - 1);
      all.add(r1.back());
    }
    CHECK_CLOSE(all.sd(), 0.1, 0.002); // marginal resolution unchanged
    CHECK_CLOSE(pearson(r1, r2), rho, 0.03);
  }
}

TEST(correlated_shift_propagates_to_unclustered_met) {
  const PriorConfig p = zmumu("{model: exponential, min_pt: 20, mean_pt: 50}",
                              "detector: {correlations: {jes: {types: [jet], scale_sigma: 0.1, unclustered: true}}}\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  for (int i = 0; i < 500; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev);
    const double d = ob.detector_params.at(0).second;
    const P4 &u = ev.particles[recoilJets(ev).at(0)].p4;
    // No invisibles: MET is the mis-measured part of the recoil.
    CHECK_CLOSE(ob.met_x, -d * u.px, 1e-8 * ev.sqrt_s_hat);
    CHECK_CLOSE(ob.met_y, -d * u.py, 1e-8 * ev.sqrt_s_hat);
  }
}

TEST(correlation_config_validated) {
  CHECK_THROWS_AS(dijet("{}", "{types: [jet], scale_sigma: 0.1}\n    other: {types: [jet]}"), std::invalid_argument,
                  "more than one correlation group");
  CHECK_THROWS_AS(dijet("{}", "{types: [], scale_sigma: 0.1}"), std::invalid_argument, "non-empty");
  CHECK_THROWS_AS(dijet("{}", "{types: [jet], scael_sigma: 0.1}"), std::invalid_argument, "unknown key");
  CHECK_THROWS_AS(sampleEnvironment(dijet("{}", "{types: [jet], resolution_correlation: 1.5}").detector, uint64_t(1), "e"),
                  std::invalid_argument, "resolution_correlation must be in [0, 1]");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}}\nproduction: [q]\ndetector: {correlation: {}}\n"), std::invalid_argument,
                  "unknown key 'correlation'");
}

TEST(poisson_distribution) {
  Rng rng(3);
  Stats s;
  const Distribution d = Distribution::poisson(4.2);
  for (int i = 0; i < 50000; ++i) {
    const double x = d.sample(rng);
    CHECK(x == std::floor(x) && x >= 0);
    s.add(x);
  }
  CHECK_CLOSE(s.mean(), 4.2, 4 * s.err());
  CHECK_CLOSE(s.sd() * s.sd(), 4.2, 0.15);
  CHECK_THROWS_AS(Distribution::poisson(-1).validate(), std::invalid_argument, "mean >= 0");
}
