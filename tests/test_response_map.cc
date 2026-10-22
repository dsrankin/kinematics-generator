// Binned (|eta|, pT) detector response maps.
#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

namespace {

PriorConfig dijets(const std::string &jet) {
  return kgtest::prior(R"(
particles:
  q: {pdg: 1}
production: [q, q]
initial_state: {sqrt_s: {model: fixed, value: 400}, rapidity: {dist: uniform, min: -2, max: 2, scope: event}}
detector:
  objects:
    jet: )" + jet + "\n");
}

const char *kMap = R"({response_map: {eta_edges: [0, 1.5, 5], pt_edges: [0, 1000],
      scale: [[0.9], [1.1]], resolution: [[0.05], [0.2]], efficiency: [[1.0], [0.5]]}})";

} // namespace

TEST(response_map_bins_and_lookup) {
  ResponseMap m;
  m.eta_edges = {0, 1.5, 3};
  m.pt_edges = {10, 50, 100};
  CHECK(m.etaBin(0.2) == 0 && m.etaBin(1.6) == 1 && m.etaBin(9.0) == 1 && m.etaBin(1.5) == 1);
  CHECK(m.ptBin(5) == 0 && m.ptBin(60) == 1 && m.ptBin(1e6) == 1);
}

TEST(response_map_applied_per_eta_region) {
  const TaskConfig t = sampleTask(dijets(kMap), 1, 0);
  const DetectorModel det(t);
  Stats central, forward;
  int nCentral = 0, nForward = 0, seenForward = 0;
  for (int i = 0; i < 20000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    for (int f : ev.finalState())
      if (std::fabs(ev.particles[f].p4.eta()) >= 1.5) ++nForward;
    for (const auto &o : det.apply(ev).objects) {
      const P4 &tp = ev.particles[o.truth_index].p4;
      const double r = o.pt() / tp.pt();
      if (std::fabs(tp.eta()) < 1.5) {
        central.add(r);
        ++nCentral;
      } else {
        forward.add(r);
        ++seenForward;
      }
    }
  }
  CHECK_CLOSE(central.mean(), 0.9, 4 * central.err());
  CHECK_CLOSE(central.sd() / 0.9, 0.05, 0.002);
  CHECK_CLOSE(forward.mean(), 1.1, 4 * forward.err());
  CHECK_CLOSE(forward.sd() / 1.1, 0.2, 0.006);
  CHECK_CLOSE(double(seenForward) / nForward, 0.5, 0.02);
}

TEST(response_map_resolution_scale_and_serialisation) {
  const std::string jet = std::string(kMap).substr(0, std::string(kMap).size() - 1) + ", resolution_scale: 2.0}";
  const TaskConfig t = sampleTask(dijets(jet), 1, 0);
  Stats central;
  for (int i = 0; i < 10000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    for (const auto &o : DetectorModel(t).apply(ev).objects)
      if (std::fabs(ev.particles[o.truth_index].p4.eta()) < 1.5) central.add(o.pt() / ev.particles[o.truth_index].p4.pt());
  }
  CHECK_CLOSE(central.sd() / 0.9, 0.10, 0.004);
  // Environments with maps round-trip exactly.
  const Environment &e = t.environment;
  const Environment r = environmentFromYAML(environmentToYAML(e));
  CHECK(environmentToYAML(r) == environmentToYAML(e));
  const TruthEvent ev = generateTruthEvent(t, 5);
  CHECK(observedEventToYAML(DetectorModel(t).apply(ev, e)) == observedEventToYAML(DetectorModel(t).apply(ev, r)));
}

TEST(response_map_validation) {
  CHECK_THROWS_AS(dijets("{response_map: {eta_edges: [0, 2], pt_edges: [0, 10], scale: [[1.0, 1.0]], resolution: [[0.1]]}}"),
                  std::invalid_argument, "pT columns");
  CHECK_THROWS_AS(dijets("{response_map: {eta_edges: [2, 1], pt_edges: [0, 10], scale: [[1.0]], resolution: [[0.1]]}}"),
                  std::invalid_argument, "must increase");
  CHECK_THROWS_AS(dijets("{response_map: {eta_edges: [0, 2], pt_edges: [0, 10], scale: [[1.0]], resolution: [[0.1]], efficiency: [[1.5]]}}"),
                  std::invalid_argument, "out of range");
  CHECK_THROWS_AS(dijets("{response_map: {eta_edges: [0, 2], pt_edges: [0, 10], scale: [[1.0]]}}"), std::invalid_argument,
                  "needs eta_edges");
}
