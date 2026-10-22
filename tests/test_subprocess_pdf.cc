// Per-event production subprocesses (hard extra jets) and the PDF-like parton luminosity.
#include <map>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

namespace {

const char *kWJets = R"(
particles:
  W: {pdg: 24, mass: 80.4, width: 2.1, decays: [[l, v]]}
  l: {pdg: 13}
  v: {pdg: 14}
  q: {pdg: 1}
processes:
  - name: W
    subprocesses:
      - {production: [W], weight: 0.6}
      - {production: [W, q], weight: 0.3}
      - {production: [W, q, q], weight: 0.1, angular: isotropic}
initial_state:
  sqrt_s: {model: pdf, ebeam: 6500, a: 1.3, b: 5, threshold_offset: 20}
)";

} // namespace

TEST(subprocesses_chosen_per_event) {
  const TaskConfig t = sampleTask(kgtest::prior(kWJets), 1, 0);
  CHECK(t.components.size() == 1 && t.components[0].subprocesses.size() == 3);
  CHECK(t.topology() == "W: {W>(l,v) | W>(l,v) + q | W>(l,v) + q + q}");
  std::map<int, int> count;
  const int n = 20000;
  for (int i = 0; i < n; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const size_t nprod = ev.producedParticles().size();
    CHECK(nprod == size_t(ev.subprocess + 1));
    count[ev.subprocess]++;
    if (nprod > 1) {
      double msum = 0;
      for (int k : ev.producedParticles()) msum += ev.particles[k].mass;
      CHECK(ev.sqrt_s_hat >= msum + 20 - 1e-9);
    }
    const P4 in = ev.incomingSum(), out = ev.finalSum();
    CHECK_CLOSE(in.e, out.e, 1e-8 * in.e);
  }
  CHECK_CLOSE(count[0] / double(n), 0.6, 0.02);
  CHECK_CLOSE(count[2] / double(n), 0.1, 0.01);
  // Serialisation round trip keeps subprocesses and reproduces events.
  const TaskConfig r = taskFromYAML(taskToYAML(t));
  CHECK(r.components[0].subprocesses.size() == 3 && r.topology() == t.topology());
  for (int i = 0; i < 20; ++i) CHECK(truthEventToYAML(generateTruthEvent(t, i)) == truthEventToYAML(generateTruthEvent(r, i)));
  // Single-production processes still serialise in the compact form.
  CHECK(taskToYAML(sampleTask(kgtest::prior(kgtest::kTwoLeptons), 1, 0)).find("subprocesses") == std::string::npos);
}

TEST(subprocess_validation) {
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}}\nprocesses:\n  - {name: x, production: [q], subprocesses: [{production: [q]}]}\n"),
                  std::invalid_argument, "either 'subprocesses'");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}}\nprocesses:\n  - {name: x, subprocesses: []}\n"), std::invalid_argument,
                  "non-empty");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}}\nprocesses:\n  - {name: x, subprocesses: [{production: [q, q, q, q]}]}\n"),
                  std::invalid_argument, "at most three");
  const PriorConfig p = kgtest::prior("particles: {q: {pdg: 1}}\nprocesses:\n  - {name: x, subprocesses: "
                                      "[{production: [q, q], weight: 0}, {production: [q, q, q], weight: 0}]}\ntask: {validation: error}\n");
  CHECK_THROWS_AS(sampleTask(p, 1, 0), InvalidTask, "positive sum");
}

TEST(pdf_luminosity_distributions) {
  // Pair production of massless partons above a 100 GeV threshold: for b = 0, ln(tau) has density
  // ~ tau^(1-a) * (-ln tau) on [ln tau0, 0]; compare the mean with numerical integration.
  const double a = 1.3, thr = 100, eb = 6500;
  const PriorConfig p = kgtest::prior(R"(
particles: {q: {pdg: 1}}
production: [q, q]
initial_state: {sqrt_s: {model: pdf, ebeam: 6500, a: 1.3, b: 0, threshold_offset: 100}}
)");
  const TaskConfig t = sampleTask(p, 1, 0);
  Stats lt, yy;
  for (int i = 0; i < 40000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    CHECK(ev.sqrt_s_hat >= thr - 1e-9);
    lt.add(2 * std::log(ev.sqrt_s_hat / (2 * eb)));
    yy.add(0.5 * std::log(ev.x1 / ev.x2));
  }
  const double l0 = 2 * std::log(thr / (2 * eb));
  double num = 0, den = 0;
  for (int k = 0; k < 20000; ++k) {
    const double l = l0 + (0 - l0) * (k + 0.5) / 20000;
    const double w = std::exp((1 - a) * l) * (-l);
    num += w * l;
    den += w;
  }
  CHECK_CLOSE(lt.mean(), num / den, 4 * lt.err() + 0.01);
  CHECK_CLOSE(yy.mean(), 0, 4 * yy.err());
}

TEST(pdf_rapidity_and_heavy_states) {
  // Single resonance: |y| bounded by the kinematic limit, broader for lighter states.
  Stats y100, y2000;
  for (double m : {100.0, 2000.0}) {
    const PriorConfig p = kgtest::prior(std::string("particles:\n  R: {pdg: 32, mass: ") + std::to_string(m) +
                                        ", decays: [[l, l]]}\n  l: {pdg: 13}\nproduction: [R]\n"
                                        "initial_state: {sqrt_s: {model: pdf, ebeam: 6500, a: 1.3, b: 5}}\n");
    const TaskConfig t = sampleTask(p, 1, 0);
    for (int i = 0; i < 5000; ++i) {
      const TruthEvent ev = generateTruthEvent(t, i);
      const double y = 0.5 * std::log(ev.x1 / ev.x2);
      CHECK(std::fabs(y) <= -std::log(m / 13000.0) + 1e-6);
      (m < 1000 ? y100 : y2000).add(y * y);
    }
  }
  CHECK(y100.mean() > 2 * y2000.mean());
  // Very heavy pair production samples efficiently (exact inverse CDF, no rejection on sqrt(s)).
  const PriorConfig h = kgtest::prior(R"(
particles: {X: {pdg: 9, mass: 3000, object: jet}}
production: [X, X]
initial_state: {sqrt_s: {model: pdf, ebeam: 6500, a: 1.3, b: 5}, recoil: {model: exponential, min_pt: 10, mean_pt: 30}}
)");
  const TaskConfig th = sampleTask(h, 1, 0);
  for (int i = 0; i < 2000; ++i) CHECK(generateTruthEvent(th, i).sqrt_s_hat >= 6000);
  CHECK_THROWS_AS(sampleTask(kgtest::prior("particles: {X: {pdg: 9, mass: 7000, object: jet}}\nproduction: [X, X]\n"
                                           "initial_state: {sqrt_s: {model: pdf}}\ntask: {validation: error}\n"),
                             1, 0),
                  InvalidTask, "exceeds the collider energy");
}
