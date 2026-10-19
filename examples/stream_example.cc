// In-memory streaming example (no intermediate files):
// sample a randomised task, generate thousands of events from it, observe each
// event in several detector environments and inspect the variable-length
// objects together with truth and task metadata.
//
//   ./example_stream [prior.yaml] [n_events]
#include <cstdlib>
#include <iostream>
#include <map>

#include "kg/kinegen.h"

int main(int argc, char **argv) {
  const std::string path = argc > 1 ? argv[1] : KG_EXAMPLES_DIR "/randomized_tasks.yaml";
  const int nEvents = argc > 2 ? std::atoi(argv[2]) : 5000;

  kg::PriorConfig prior = kg::loadPrior(path);
  kg::TaskConfig task = kg::sampleTask(prior, /*master_seed=*/1, /*task_index=*/3);
  std::cout << "task " << task.task_id << ": " << task.topology() << "\n";
  for (const auto &kv : task.parameterValues())
    if (kv.first.find(".mass") != std::string::npos) std::cout << "  " << kv.first << " = " << kv.second << "\n";
  std::cout << "  features:";
  for (const auto &f : *task.feature_names) std::cout << " " << f;
  std::cout << "\n";

  // Three environments: the task's own and two domain-randomised ones.
  std::vector<kg::Environment> envs = {task.environment, kg::sampleEnvironment(prior.detector, 99, int64_t(0)),
                                       kg::sampleEnvironment(prior.detector, 99, int64_t(1))};
  kg::DetectorModel detector(task);

  std::map<std::string, double> nObj;
  std::map<size_t, int> multiplicity;
  for (int i = 0; i < nEvents; ++i) {
    kg::TruthEvent truth = kg::generateTruthEvent(task); // event i (task.next_event)
    multiplicity[truth.finalState().size()]++;
    for (const auto &env : envs) {
      kg::ObservedEvent reco = detector.apply(truth, env);
      nObj[env.env_id] += reco.objects.size();
      if (i == 0) {
        std::cout << "event 0 in " << env.env_id << ": " << reco.objects.size() << " objects, MET " << reco.met() << "\n";
        for (const auto &o : reco.objects) {
          std::cout << "    " << o.type << " pt=" << o.pt() << " eta=" << o.eta();
          if (o.truth_index >= 0) std::cout << " <- " << truth.particles[o.truth_index].name << " (ancestors:";
          if (o.truth_index >= 0)
            for (int a : truth.ancestors(o.truth_index)) std::cout << " " << truth.particles[a].name;
          std::cout << (o.truth_index >= 0 ? ")" : " fake") << "\n";
        }
      }
    }
  }
  for (const auto &kv : nObj) std::cout << kv.first << ": <n_objects> = " << kv.second / nEvents << "\n";
  std::cout << "truth final-state multiplicities:";
  for (const auto &kv : multiplicity) std::cout << " " << kv.first << ":" << kv.second;
  std::cout << "\n";
  return 0;
}
