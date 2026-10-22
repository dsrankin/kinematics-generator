// Configuration and kinematic validation: invalid priors must fail loudly.
#include "helpers.h"
#include "test_framework.h"

using namespace kg;

TEST(cyclic_decay_graph_rejected) {
  CHECK_THROWS_AS(kgtest::prior(R"(
particles:
  A: {pdg: 1, mass: 100, object: jet, decays: [[B, q]]}
  B: {pdg: 2, mass: 50, object: jet, decays: [[C, q]]}
  C: {pdg: 3, mass: 20, object: jet, decays: [[A, q]]}
  q: {pdg: 1}
production: [A]
)"),
                  std::invalid_argument, "cyclic decay chain: A -> B -> C -> A");
}

TEST(unknown_daughter_and_particle_rejected) {
  CHECK_THROWS_AS(kgtest::prior("particles: {A: {mass: 10, decays: [[x, y]]}}\nproduction: [A]\n"), std::invalid_argument,
                  "unknown daughter 'x'");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}}\nproduction: [Z]\n"), std::invalid_argument, "unknown particle 'Z'");
}

TEST(branching_fractions_validated) {
  const std::string head = "particles:\n  q: {pdg: 1}\n  l: {pdg: 11}\n  A:\n    mass: 100\n    decays:\n";
  CHECK_THROWS_AS(kgtest::prior(head + "      - {daughters: [q, q], br: 0.5}\n      - {daughters: [l, l], br: 0.4}\nproduction: [A]\n"),
                  std::invalid_argument, "sum to 0.9");
  CHECK_THROWS_AS(kgtest::prior(head + "      - {daughters: [q, q], br: 1.5}\nproduction: [A]\n"), std::invalid_argument, "outside");
  CHECK_THROWS_AS(kgtest::prior(head + "      - [q, q]\n      - [l, l]\nproduction: [A]\n"), std::invalid_argument, "explicit 'br'");
  // Explicit normalisation is allowed.
  const PriorConfig p = kgtest::prior("particles:\n  q: {pdg: 1}\n  l: {pdg: 11}\n  A:\n    mass: 100\n    normalize_br: true\n"
                                      "    decays:\n      - {daughters: [q, q], br: 3}\n      - {daughters: [l, l], br: 1}\nproduction: [A]\n");
  const TaskConfig t = sampleTask(p, 1, 0);
  CHECK_CLOSE(t.particles[t.species("A")].modes[0].br, 0.75, 1e-12);
  // Sampled branching fractions that do not sum to one invalidate the task.
  const PriorConfig r = kgtest::prior(head + "      - {daughters: [q, q], br: {dist: uniform, min: 0, max: 1}}\n"
                                             "      - {daughters: [l, l], br: 0.5}\nproduction: [A]\ntask: {validation: error}\n");
  CHECK_THROWS_AS(sampleTask(r, 1, 0), InvalidTask, "sum to");
}

TEST(multiplicity_and_structure_rejected) {
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}, A: {mass: 100, decays: [[q, q, q, q]]}}\nproduction: [A]\n"),
                  std::invalid_argument, "two- and three-body");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}}\nproduction: [q, q, q, q]\n"), std::invalid_argument, "at most three");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1, stable: false}}\nproduction: [q]\n"), std::invalid_argument, "contradicts");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1, mas: 3}}\nproduction: [q]\n"), std::invalid_argument, "unknown key 'mas'");
  CHECK_THROWS_AS(kgtest::prior("particles: {X: {pdg: 9999}}\nproduction: [X]\n"), std::invalid_argument, "cannot infer object type");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}}\nproduction: [q]\ndetector: {objects: {jet: {efficency: 1}}}\n"),
                  std::invalid_argument, "unknown key 'efficency'");
  CHECK_THROWS_AS(kgtest::prior("particles: {q: {pdg: 1}}\nproduction: [q]\nprocesses: []\n"), std::invalid_argument, "either");
}

TEST(kinematically_forbidden_decay_rejected) {
  const std::string y = R"(
particles:
  X: {pdg: 1, mass: 100, object: jet, decays: [[Y, Y]]}
  Y: {pdg: 2, mass: 60, object: jet}
production: [X]
)";
  CHECK_THROWS_AS(sampleTask(kgtest::prior(y + "task: {validation: error}\n"), 1, 0), InvalidTask,
                  "kinematically forbidden: m(X)=100 < 120");
  // With rejection the prior has no valid task at all: still an error, not a silent repair.
  CHECK_THROWS_AS(sampleTask(kgtest::prior(y + "task: {max_attempts: 20}\n"), 1, 0), InvalidTask, "no valid task after 20 attempts");
}

TEST(rejection_sampling_respects_thresholds) {
  const PriorConfig p = kgtest::prior(R"(
particles:
  X: {pdg: 1, mass: {dist: uniform, min: 50, max: 200}, object: jet, decays: [[Y, Y]]}
  Y: {pdg: 2, mass: 60, object: jet}
production: [X]
)");
  int multi = 0;
  for (int i = 0; i < 200; ++i) {
    const TaskConfig t = sampleTask(p, 3, i);
    CHECK(t.param("X.mass") >= 120);
    multi += t.attempts > 1;
  }
  CHECK(multi > 20); // P(reject) = 70/150 per draw
}

TEST(off_shell_daughters_allowed_with_width) {
  // H -> Z Z* below the on-shell ZZ threshold is physical when the Z has a width.
  const PriorConfig p = kgtest::prior(R"(
particles:
  H: {pdg: 25, mass: 125, decays: [[Z, Z]]}
  Z: {pdg: 23, mass: 91.19, width: 2.5, decays: [[l, l]]}
  l: {pdg: 13}
production: [H]
initial_state: {sqrt_s: {model: fixed}, rapidity: 0}
)");
  const TaskConfig t = sampleTask(p, 1, 0);
  for (int i = 0; i < 500; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    double sum = 0;
    for (int c : ev.particles[ev.producedParticles()[0]].children) sum += ev.particles[c].mass;
    CHECK(sum <= 125 + 1e-9);
  }
  // Without a width the same decay is forbidden.
  const PriorConfig q = kgtest::prior(R"(
particles:
  H: {pdg: 25, mass: 125, decays: [[Z, Z]]}
  Z: {pdg: 23, mass: 91.19, decays: [[l, l]]}
  l: {pdg: 13}
production: [H]
task: {validation: error}
)");
  CHECK_THROWS_AS(sampleTask(q, 1, 0), InvalidTask, "kinematically forbidden");
}

TEST(breit_wigner_cutoff_validated) {
  const PriorConfig p = kgtest::prior(R"(
bw_cutoff: 5
particles:
  H: {pdg: 25, mass: 125, decays: [[Z, Z]]}
  Z: {pdg: 23, mass: 91.19, width: 2.5, decays: [[l, l]]}
  l: {pdg: 13}
production: [H]
task: {validation: error}
)");
  // Each Z needs m >= 91.19 - 5*2.5 = 78.69, so ZZ needs 157 > 125.
  CHECK_THROWS_AS(sampleTask(p, 1, 0), InvalidTask, "kinematically forbidden");
}

TEST(event_level_threshold_violation_throws) {
  const PriorConfig p = kgtest::prior(R"(
particles:
  X: {pdg: 1, mass: {dist: uniform, min: 50, max: 200, scope: event}, object: jet, decays: [[Y, Y]]}
  Y: {pdg: 2, mass: 60, object: jet}
production: [X]
task: {validation: error, max_attempts: 1}
)");
  // The task is only valid if its initial draw is; events then fail loudly when below threshold.
  bool sawFailure = false;
  for (int i = 0; i < 20 && !sawFailure; ++i) {
    TaskConfig t;
    try {
      t = sampleTask(p, 11, i);
    } catch (const InvalidTask &) {
      continue;
    }
    for (int e = 0; e < 50; ++e) {
      try {
        generateTruthEvent(t, e);
      } catch (const std::runtime_error &err) {
        CHECK(std::string(err.what()).find("below threshold") != std::string::npos);
        sawFailure = true;
        break;
      }
    }
  }
  CHECK(sawFailure);
}

TEST(detector_and_feature_config_validated) {
  const std::string base = "particles: {q: {pdg: 1}}\nproduction: [q]\n";
  CHECK_THROWS_AS(kgtest::prior(base + "detector: {objects: {jet: {misid: {jet: 0.1}}}}\n"), std::invalid_argument, "own type");
  CHECK_THROWS_AS(kgtest::prior(base + "features: [{name: f, kind: truth_proxy, quantity: spin}]\n"), std::invalid_argument,
                  "unknown quantity");
  CHECK_THROWS_AS(kgtest::prior(base + "features: [{name: f, kind: environment, param: detector.jet.nothing}]\n"),
                  std::invalid_argument, "unknown detector parameter");
  CHECK_THROWS_AS(kgtest::prior(base + "features: [{name: f, kind: gaussian, sigma: {dist: uniform, min: 0, max: 1, scope: event}}]\n"),
                  std::invalid_argument, "event-scoped");
  // Sampled detector values outside their physical range are errors, not clipped.
  const PriorConfig p = kgtest::prior(base + "detector: {objects: {jet: {efficiency: {dist: gaussian, mean: 1.0, sigma: 0.5}}}}\n"
                                             "task: {validation: error}\n");
  bool threw = false;
  for (int i = 0; i < 10 && !threw; ++i) {
    try {
      sampleEnvironment(p.detector, uint64_t(i), "e");
    } catch (const std::invalid_argument &e) {
      threw = std::string(e.what()).find("efficiency must be in [0, 1]") != std::string::npos;
    }
  }
  CHECK(threw);
}

TEST(misspelled_nested_keys_rejected) {
  // Typos inside parameter, detector, multiplicity and fraction maps are errors, not silent defaults.
  const std::string base = kgtest::kTwoLeptons;
  auto bad = [&](const std::string &extra, const char *msg) {
    CHECK_THROWS_AS(kgtest::prior(base + extra), std::invalid_argument, msg);
  };
  bad("parameters: {x: {dist: uniform, min: 0, max: 1, scop: event}}\n", "unknown key 'scop'");
  bad("parameters: {x: {expr: \"1 + 1\", scale: 2}}\n", "unknown key 'scale'");
  bad("detector: {objects: {muon: {resolution: {stochastc: 0.1}}}}\n", "unknown key 'stochastc'");
  bad("detector: {objects: {muon: {efficiency: {plateu: 0.9}}}}\n", "unknown key 'plateu'");
  bad("detector: {fakes: {rat: 0.1}}\n", "unknown key 'rat'");
  bad("detector: {met: {soft_sgima: 3}}\n", "unknown key 'soft_sgima'");
  CHECK_THROWS_AS(kgtest::prior(R"(
particles:
  X: {pdg: 13, mass: 0.1}
processes:
  - {name: m, production: [X], multiplicity: {values: [1, 2], weigths: [1, 1]}}
)"), std::invalid_argument, "unknown key 'weigths'");
  bad("task: {components: {dist: categorical, values: [1], scope: event}}\n", "unknown key 'scope'");
  bad("features: [{name: f, kind: label_correlated, label: component, shift: {dist: gaussian, mean: 0, sigma: 1, ratio_of: x}}]\n",
      "unknown key 'ratio_of'");
  // Correct spellings still parse.
  kgtest::prior(base + "parameters: {x: {dist: uniform, min: 0, max: 1, scope: event}, y: {expr: \"2 * x\", scope: event}}\n");
}
