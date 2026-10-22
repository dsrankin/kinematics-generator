// Auxiliary features: noise, proxies, spurious correlations, missing values, schemas.
#include <map>
#include <set>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

namespace {

PriorConfig withFeatures(const std::string &features, const std::string &extra = "") {
  return kgtest::prior(R"(
particles:
  X: {pdg: 32, mass: 400, decays: [[q, b]]}
  W: {pdg: 24, mass: 80.4, decays: [[q, q]]}
  q: {pdg: 1}
  b: {pdg: 5, mass: 4.7}
processes:
  - {name: sig, production: [X]}
  - {name: bkg, production: [W, W]}
)" + extra + "features:\n" + features);
}

} // namespace

TEST(independent_noise_features) {
  const PriorConfig p = withFeatures(R"(
  - {name: g, kind: gaussian, mean: 2, sigma: 3}
  - {name: u, kind: uniform, min: -1, max: 5, applies_to: [bjet]}
)", "task: {components: 2}\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  CHECK(t.feature_names->size() == 2);
  const DetectorModel det(t);
  Stats g, u;
  int nonb = 0;
  for (int i = 0; i < 5000; ++i) {
    const ObservedEvent ob = det.apply(generateTruthEvent(t, i));
    for (size_t k = 0; k < ob.objects.size(); ++k) {
      g.add(ob.feature(k, "g"));
      const double x = ob.feature(k, "u");
      if (ob.objects[k].type == "bjet") {
        CHECK(x >= -1 && x < 5);
        u.add(x);
      } else {
        CHECK(std::isnan(x)); // not applicable -> missing
        ++nonb;
      }
    }
  }
  CHECK(nonb > 0 && u.n > 1000);
  CHECK_CLOSE(g.mean(), 2, 4 * g.err());
  CHECK_CLOSE(g.sd(), 3, 0.06);
  CHECK_CLOSE(u.mean(), 2, 4 * u.err());
}

TEST(truth_and_observed_proxies) {
  const PriorConfig p = withFeatures(R"(
  - {name: true_pt, kind: truth_proxy, quantity: pt}
  - {name: parent_m, kind: truth_proxy, quantity: parent_mass, abs_sigma: 5}
  - {name: obs_eta, kind: observed_proxy, quantity: eta, rel_sigma: 0}
  - {name: fake_flag, kind: truth_proxy, quantity: is_fake}
)", "detector: {objects: {jet: {resolution: {constant: 0.2}}}, fakes: {rate: 0.5, types: {jet: 1}}}\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  Stats dm;
  int fakes = 0;
  for (int i = 0; i < 3000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev);
    for (size_t k = 0; k < ob.objects.size(); ++k) {
      const ObservedObject &o = ob.objects[k];
      CHECK_CLOSE(ob.feature(k, "obs_eta"), o.eta(), 1e-12);
      if (o.is_fake) {
        ++fakes;
        CHECK(std::isnan(ob.feature(k, "true_pt")));
        CHECK(ob.feature(k, "fake_flag") == 1);
        continue;
      }
      const TruthParticle &tp = ev.particles[o.truth_index];
      CHECK_CLOSE(ob.feature(k, "true_pt"), tp.p4.pt(), 1e-9);
      CHECK(ob.feature(k, "fake_flag") == 0);
      dm.add(ob.feature(k, "parent_m") - ev.particles[tp.parent].mass);
    }
  }
  CHECK(fakes > 100);
  CHECK_CLOSE(dm.mean(), 0, 4 * dm.err());
  CHECK_CLOSE(dm.sd(), 5, 0.15);
}

TEST(missing_values_and_task_dependent_schema) {
  const PriorConfig p = withFeatures(R"(
  - {name: always, kind: gaussian, missing_prob: 0.3}
  - {name: sometimes, kind: gaussian, presence_prob: 0.5}
)");
  int withSometimes = 0;
  for (int i = 0; i < 200; ++i) {
    const TaskConfig t = sampleTask(p, 2, i);
    CHECK((*t.feature_names)[0] == "always");
    withSometimes += t.feature_names->size() == 2;
  }
  CHECK(withSometimes > 70 && withSometimes < 130);
  const TaskConfig t = sampleTask(p, 2, 0);
  const DetectorModel det(t);
  int miss = 0, n = 0;
  for (int i = 0; i < 5000; ++i) {
    const ObservedEvent ob = det.apply(generateTruthEvent(t, i));
    for (size_t k = 0; k < ob.objects.size(); ++k, ++n) miss += std::isnan(ob.feature(k, "always"));
  }
  CHECK_CLOSE(double(miss) / n, 0.3, 4 * std::sqrt(0.21 / n));
}

TEST(spurious_label_correlations) {
  const PriorConfig p = withFeatures(R"(
  - {name: spur, kind: label_correlated, label: component, sigma: 0.5, shift: {dist: uniform, min: -3, max: 3}}
  - {name: tag, kind: label_correlated, label: truth_type, sigma: 0.1, shift: {dist: gaussian, mean: 0, sigma: 2}}
)", "task: {components: 2, fractions: equal}\n");
  const TaskConfig t = sampleTask(p, 3, 0);
  const auto &shift = t.features[0].shifts;
  CHECK(shift.size() == 2);
  const DetectorModel det(t);
  std::map<int, Stats> byComponent;
  std::map<std::string, Stats> byTruth;
  for (int i = 0; i < 4000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev);
    for (size_t k = 0; k < ob.objects.size(); ++k) {
      byComponent[ob.component].add(ob.feature(k, "spur"));
      byTruth[ev.particles[ob.objects[k].truth_index].name].add(ob.feature(k, "tag"));
    }
  }
  CHECK(byComponent.size() == 2);
  for (auto &kv : byComponent) CHECK_CLOSE(kv.second.mean(), shift.at(std::to_string(kv.first)), 4 * kv.second.err());
  for (auto &kv : byTruth) CHECK_CLOSE(kv.second.mean(), t.features[1].shifts.at(kv.first), 4 * kv.second.err());
  // The correlation strength itself is randomised between tasks.
  CHECK(sampleTask(p, 3, 1).features[0].shifts.at("0") != shift.at("0"));
}

TEST(environment_correlated_feature) {
  const PriorConfig p = withFeatures(R"(
  - {name: res, kind: environment, param: detector.jet.resolution.stochastic, coeff: 2, offset: 1}
)", "detector: {objects: {jet: {resolution: {stochastic: {dist: uniform, min: 0.5, max: 1.5}}}}}\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  for (int k = 0; k < 5; ++k) {
    const Environment env = sampleEnvironment(p.detector, 9, int64_t(k));
    const ObservedEvent ob = det.apply(generateTruthEvent(t, 0), env);
    CHECK(!ob.objects.empty());
    CHECK_CLOSE(ob.feature(0, "res"), 1 + 2 * env.value("detector.jet.resolution.stochastic"), 1e-12);
  }
}
