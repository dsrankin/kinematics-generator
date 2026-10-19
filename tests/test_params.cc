// Distributions, expressions, parameter sets and task-level sampling.
#include <set>

#include <yaml-cpp/yaml.h>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

TEST(distribution_moments) {
  Rng rng(11);
  const int n = 100000;
  auto stats = [&](const Distribution &d) {
    Stats s;
    for (int i = 0; i < n; ++i) s.add(d.sample(rng));
    return s;
  };
  Stats u = stats(Distribution::uniform(2, 6));
  CHECK_CLOSE(u.mean(), 4, 4 * u.err());
  CHECK_CLOSE(u.sd(), 4 / std::sqrt(12.0), 0.01);
  Stats g = stats(Distribution::gaussian(3, 2));
  CHECK_CLOSE(g.mean(), 3, 4 * g.err());
  CHECK_CLOSE(g.sd(), 2, 0.02);
  // log-uniform: log x uniform.
  Stats lu;
  const Distribution l = Distribution::logUniform(1, 100);
  for (int i = 0; i < n; ++i) {
    const double x = l.sample(rng);
    CHECK(x >= 1 && x <= 100);
    lu.add(std::log(x));
  }
  CHECK_CLOSE(lu.mean(), 0.5 * std::log(100.0), 4 * lu.err());
  const Distribution tg = Distribution::truncatedGaussian(0, 1, 0.5, 3);
  for (int i = 0; i < 10000; ++i) {
    const double x = tg.sample(rng);
    CHECK(x >= 0.5 && x <= 3);
  }
  // Far-tail truncation uses CDF inversion.
  const Distribution tail = Distribution::truncatedGaussian(0, 1, 5, 6);
  for (int i = 0; i < 1000; ++i) {
    const double x = tail.sample(rng);
    CHECK(x >= 5 && x <= 6);
  }
  Stats e = stats(Distribution::exponential(10, 5));
  CHECK_CLOSE(e.mean(), 15, 4 * e.err());
  // Breit-Wigner: median at the pole, half the events within +-Gamma/2.
  const Distribution bw = Distribution::breitWigner(91, 2.5);
  int inside = 0;
  for (int i = 0; i < n; ++i) inside += std::fabs(bw.sample(rng) - 91) < 1.25;
  CHECK_CLOSE(double(inside) / n, 0.5, 0.01);
  const Distribution bwt = Distribution::breitWigner(91, 2.5, 85, 95, true);
  for (int i = 0; i < 1000; ++i) {
    const double x = bwt.sample(rng);
    CHECK(x >= 85 && x <= 95);
  }
  // Categorical frequencies.
  const Distribution c = Distribution::categorical({1, 2, 3}, {1, 2, 7});
  int count3 = 0;
  for (int i = 0; i < n; ++i) count3 += c.sample(rng) == 3;
  CHECK_CLOSE(double(count3) / n, 0.7, 0.01);
  CHECK(Distribution::fixed(4.2).sample(rng) == 4.2);
}

TEST(distribution_validation) {
  CHECK_THROWS_AS(Distribution::uniform(3, 1).validate(), std::invalid_argument, "min <= max");
  CHECK_THROWS_AS(Distribution::logUniform(0, 1).validate(), std::invalid_argument, "loguniform");
  CHECK_THROWS_AS(Distribution::categorical({1, 2}, {1}).validate(), std::invalid_argument, "length");
  CHECK_THROWS_AS(Distribution::fromYAML(YAML::Load("{dist: weird}")), std::invalid_argument, "unknown distribution");
  CHECK_THROWS_AS(Distribution::fromYAML(YAML::Load("{dist: gaussian, mean: 1}")), std::invalid_argument, "sigma");
}

TEST(expressions) {
  Expression e("2*x + max(y, 3)^2 - -1 / 2");
  std::vector<std::string> vars = e.variables();
  CHECK(vars.size() == 2 && vars[0] == "x" && vars[1] == "y");
  e.bind([](const std::string &n) { return n == "x" ? 0 : 1; });
  CHECK_CLOSE(e.eval({1.5, 2.0}), 3 + 9 + 0.5, 1e-12);
  CHECK_CLOSE(e.eval({0.0, 4.0}), 16.5, 1e-12);
  Expression f("sqrt(X.mass) * pow(2, 3)");
  f.bind([](const std::string &) { return 0; });
  CHECK_CLOSE(f.eval({16.0}), 32, 1e-12);
  CHECK_THROWS_AS(Expression("1 +"), std::invalid_argument, "unexpected");
  CHECK_THROWS_AS(Expression("foo(1)"), std::invalid_argument, "unknown function");
}

TEST(param_dependencies_and_scope_propagation) {
  ParamSet ps;
  ps.parse(YAML::Load("{dist: uniform, min: 0, max: 1, scope: event}"), "a");
  ps.parse(YAML::Load("2 * a + c"), "b");
  ps.parse(YAML::Load("{dist: uniform, min: 10, max: 20}"), "c");
  ps.parse(YAML::Load("{dist: uniform, min: 0.1, max: 0.5, ratio_of: c}"), "d");
  ps.compile();
  CHECK(ps.isEvent(ps.at("a")) && ps.isEvent(ps.at("b")));
  CHECK(!ps.isEvent(ps.at("c")) && !ps.isEvent(ps.at("d")));
  Rng rng(3);
  std::vector<double> v = ps.sample(rng);
  const double c = v[ps.at("c")];
  CHECK_CLOSE(v[ps.at("b")], 2 * v[ps.at("a")] + c, 1e-12);
  CHECK(v[ps.at("d")] >= 0.1 * c && v[ps.at("d")] <= 0.5 * c);
  const double a0 = v[ps.at("a")];
  ps.resampleEvent(v, rng);
  CHECK(v[ps.at("a")] != a0);
  CHECK(v[ps.at("c")] == c); // task-level values untouched
  CHECK_CLOSE(v[ps.at("b")], 2 * v[ps.at("a")] + c, 1e-12);
  // Freezing keeps event-scoped definitions only.
  ParamSet fr = ps.frozen(v);
  CHECK(fr.def(fr.at("c")).kind == ParamDef::Kind::Fixed);
  CHECK(fr.def(fr.at("a")).kind == ParamDef::Kind::Dist);
}

TEST(param_cycles_and_unknown_references) {
  ParamSet ps;
  ps.parse(YAML::Load("b + 1"), "a");
  ps.parse(YAML::Load("a * 2"), "b");
  CHECK_THROWS_AS(ps.compile(), std::invalid_argument, "cycle");
  ParamSet q;
  q.parse(YAML::Load("nope + 1"), "a");
  CHECK_THROWS_AS(q.compile(), std::invalid_argument, "unknown parameter");
  ParamSet r;
  r.addFixed("x", 1);
  CHECK_THROWS_AS(r.addFixed("x", 2), std::invalid_argument, "twice");
}

TEST(mass_ratio_constraints_across_tasks) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("randomized_tasks.yaml"));
  for (int i = 0; i < 300; ++i) {
    const TaskConfig t = sampleTask(p, 5, i);
    const double mX = t.param("X.mass"), mY = t.param("Y.mass"), mN = t.param("N.mass");
    CHECK(mX >= 300 && mX <= 3000);
    CHECK(mY >= 0.2 * mX - 1e-9 && mY <= 0.8 * mX + 1e-9);
    CHECK(mN >= 0.05 * mY - 1e-9 && mN <= 0.9 * mY + 1e-9);
    CHECK_CLOSE(t.param("Y.width"), 0.01 * mY, 1e-9);
    const double g = t.param("X.width") / mX;
    CHECK(g >= 0.001 - 1e-12 && g <= 0.05 + 1e-12);
    CHECK_CLOSE(t.param("Zp.decay0.br") + t.param("Zp.decay1.br"), 1, 1e-12);
  }
}

TEST(task_level_topology_and_multiplicity_sampling) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("randomized_tasks.yaml"));
  std::set<std::string> processes;
  std::set<size_t> nproduced, ncomponents;
  std::set<int> ymodes;
  for (int i = 0; i < 400; ++i) {
    const TaskConfig t = sampleTask(p, 9, i);
    ncomponents.insert(t.components.size());
    double fsum = 0;
    for (const auto &c : t.components) {
      processes.insert(c.process);
      nproduced.insert(c.produced.size());
      fsum += c.fraction;
    }
    CHECK_CLOSE(fsum, 1, 1e-12);
    ymodes.insert(t.particles[t.species("Y")].fixed_mode);
  }
  CHECK(processes.size() == 4);
  CHECK(nproduced.count(1) && nproduced.count(2) && nproduced.count(3));
  CHECK(ncomponents.count(1) && ncomponents.count(2));
  CHECK(ymodes.size() == 3 && !ymodes.count(-1));
}

TEST(dirichlet_and_fixed_fractions) {
  const std::string base = std::string(kgtest::kTwoLeptons);
  auto mk = [&](const std::string &task) {
    std::string y = base;
    y.replace(y.find("production: [Z]"), 15, "processes:\n  - {name: a, production: [Z]}\n  - {name: b, production: [Z, Z]}");
    return kgtest::prior(y + task);
  };
  const PriorConfig fixed = mk("task: {components: 2, fractions: [3, 1]}\n");
  const TaskConfig t = sampleTask(fixed, 1, 0);
  CHECK(t.components.size() == 2);
  CHECK_CLOSE(t.components[0].fraction, 0.75, 1e-12);
  // Event-level component labels follow the fractions.
  int n0 = 0;
  const int n = 20000;
  for (int i = 0; i < n; ++i) n0 += generateTruthEvent(t, i).component == 0;
  CHECK_CLOSE(double(n0) / n, 0.75, 4 * std::sqrt(0.75 * 0.25 / n));
  CHECK_THROWS_AS(mk("task: {components: 3}\n"), std::invalid_argument, "exceeds");
}
