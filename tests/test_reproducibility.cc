// Determinism, serialisation round trips and dataset splits.
#include <map>
#include <set>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;

TEST(fixed_seeds_reproduce_tasks_and_events) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("randomized_tasks.yaml"));
  const TaskConfig a = sampleTask(p, 123, 4), b = sampleTask(p, 123, 4);
  CHECK(taskToYAML(a) == taskToYAML(b));
  CHECK(taskToYAML(a) != taskToYAML(sampleTask(p, 123, 5)));
  CHECK(taskToYAML(a) != taskToYAML(sampleTask(p, 124, 4)));
  const DetectorModel da(a), db(b);
  for (int i = 0; i < 50; ++i) {
    CHECK(truthEventToYAML(generateTruthEvent(a, i)) == truthEventToYAML(generateTruthEvent(b, i)));
    CHECK(observedEventToYAML(da.apply(generateTruthEvent(a, i))) == observedEventToYAML(db.apply(generateTruthEvent(b, i))));
  }
  CHECK(truthEventToYAML(generateTruthEvent(a, 0)) != truthEventToYAML(generateTruthEvent(a, 1)));
  // Task samplers are sequential views of the same streams.
  TaskSampler s(p, 123);
  for (int i = 0; i < 5; ++i) CHECK(taskToYAML(s.next()) == taskToYAML(sampleTask(p, 123, i)));
}

TEST(event_order_independence) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("fixed_process.yaml"));
  TaskConfig t = sampleTask(p, 1, 0);
  std::vector<std::string> forward;
  for (int i = 0; i < 20; ++i) forward.push_back(truthEventToYAML(generateTruthEvent(t))); // uses next_event
  CHECK(t.next_event == 20);
  for (int i = 19; i >= 0; --i) CHECK(truthEventToYAML(generateTruthEvent(t, i)) == forward[i]);
}

TEST(task_and_event_randomness_are_distinct) {
  // Changing detector or feature configuration must not change truth events,
  // and truth is identical whichever environment observes it.
  const std::string physics = std::string(kgtest::kTwoLeptons);
  const TaskConfig a = sampleTask(kgtest::prior(physics), 1, 0);
  const TaskConfig b = sampleTask(kgtest::prior(physics + "detector: {objects: {muon: {efficiency: 0.5}}}\n"
                                                           "features: [{name: f, kind: gaussian}]\n"),
                                  1, 0);
  for (int i = 0; i < 20; ++i) CHECK(truthEventToYAML(generateTruthEvent(a, i)) == truthEventToYAML(generateTruthEvent(b, i)));
  // Detector randomness is independent between environments but fixed per (environment, event).
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("detector_randomization.yaml"));
  const TaskConfig t = sampleTask(p, 1, 0);
  const Environment e0 = sampleEnvironment(p.detector, 3, int64_t(0)), e1 = sampleEnvironment(p.detector, 3, int64_t(1));
  const TruthEvent ev = generateTruthEvent(t, 7);
  const DetectorModel det(t);
  CHECK(observedEventToYAML(det.apply(ev, e0)) == observedEventToYAML(det.apply(ev, e0)));
  CHECK(observedEventToYAML(det.apply(ev, e0)) != observedEventToYAML(det.apply(ev, e1)));
}

TEST(task_serialization_round_trip) {
  for (const char *f : {"fixed_process.yaml", "randomized_tasks.yaml", "detector_randomization.yaml"}) {
    const PriorConfig p = kg::loadPrior(kgtest::examplePath(f));
    for (int it = 0; it < 10; ++it) {
      const TaskConfig t = sampleTask(p, 31, it);
      const std::string y = taskToYAML(t);
      const TaskConfig r = taskFromYAML(y);
      CHECK(taskToYAML(r) == y);
      CHECK(r.topology() == t.topology());
      const DetectorModel dt(t), dr(r);
      for (int i = 0; i < 25; ++i) {
        const TruthEvent a = generateTruthEvent(t, i), b = generateTruthEvent(r, i);
        CHECK(truthEventToYAML(a) == truthEventToYAML(b));
        CHECK(observedEventToYAML(dt.apply(a)) == observedEventToYAML(dr.apply(b)));
      }
    }
  }
}

TEST(environment_and_event_serialization) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("detector_randomization.yaml"));
  const TaskConfig t = sampleTask(p, 1, 0);
  const Environment e = sampleEnvironment(p.detector, 8, int64_t(2));
  const Environment r = environmentFromYAML(environmentToYAML(e));
  CHECK(environmentToYAML(r) == environmentToYAML(e));
  CHECK(r.values == e.values);
  const DetectorModel det(t);
  for (int i = 0; i < 20; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev, e);
    CHECK(observedEventToYAML(ob) == observedEventToYAML(det.apply(ev, r)));
    CHECK(truthEventToYAML(truthEventFromYAML(truthEventToYAML(ev))) == truthEventToYAML(ev));
    CHECK(observedEventToYAML(observedEventFromYAML(observedEventToYAML(ob))) == observedEventToYAML(ob));
  }
}

TEST(dataset_splits) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("randomized_tasks.yaml"));
  SplitSpec byTask;
  std::map<std::string, int> counts;
  std::map<std::string, std::set<std::string>> topoSplits;
  SplitSpec byTopo;
  byTopo.by = SplitSpec::By::Topology;
  SplitSpec byMass;
  byMass.by = SplitSpec::By::Parameter;
  byMass.parameter = "X.mass";
  byMass.edges = {1500};
  byMass.labels = {"train", "test"};
  for (int i = 0; i < 1000; ++i) {
    const TaskConfig t = sampleTask(p, 50, i);
    const std::string s = assignSplit(t, byTask);
    CHECK(s == assignSplit(sampleTask(p, 50, i), byTask) || i > 50); // deterministic (spot-check first tasks)
    counts[s]++;
    topoSplits[t.topology()].insert(assignSplit(t, byTopo));
    CHECK(assignSplit(t, byMass) == (t.param("X.mass") < 1500 ? "train" : "test"));
  }
  CHECK_CLOSE(counts["train"] / 1000.0, 0.8, 0.05);
  CHECK_CLOSE(counts["test"] / 1000.0, 0.1, 0.04);
  for (const auto &kv : topoSplits) CHECK(kv.second.size() == 1); // a topology never straddles splits
  SplitSpec byEnv;
  byEnv.by = SplitSpec::By::Environment;
  const Environment e = sampleEnvironment(p.detector, 1, int64_t(3));
  CHECK(assignSplit(e, byEnv) == assignSplit(sampleEnvironment(p.detector, 1, int64_t(3)), byEnv));
}
