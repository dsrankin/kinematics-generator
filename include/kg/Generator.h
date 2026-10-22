// Event generation: truth events from tasks, observed events from truth.
//
//   PriorConfig prior = kg::loadPrior("prior.yaml");
//   TaskConfig task = kg::sampleTask(prior, seed, task_index);
//   kg::DetectorModel detector(task);
//   for (uint64_t i = 0; i < 1000; ++i) {
//     TruthEvent truth = kg::generateTruthEvent(task, i);
//     ObservedEvent reco = detector.apply(truth, task.environment);
//   }
//
// Event i of a task is a pure function of (task, i); detector output is a pure
// function of (task, truth, environment). Events can therefore be generated in
// any order, in parallel, and re-generated on demand.
#ifndef KG_GENERATOR_H
#define KG_GENERATOR_H

#include "kg/Event.h"
#include "kg/Task.h"

namespace kg {

TruthEvent generateTruthEvent(const TaskConfig &task, uint64_t event_index);
// Uses and increments task.next_event.
TruthEvent generateTruthEvent(TaskConfig &task);

class DetectorModel {
 public:
  explicit DetectorModel(const TaskConfig &task) : task_(&task) {}
  ObservedEvent apply(const TruthEvent &truth, const Environment &env) const;
  ObservedEvent apply(const TruthEvent &truth) const { return apply(truth, task_->environment); }

 private:
  const TaskConfig *task_;
};

struct GeneratedEvent {
  TruthEvent truth;
  ObservedEvent observed;
};

GeneratedEvent generateEvent(const TaskConfig &task, const Environment &env, uint64_t event_index);

} // namespace kg

#endif
