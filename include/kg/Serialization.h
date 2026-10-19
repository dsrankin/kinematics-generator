// YAML serialisation of truth and observed events (tasks: see Task.h,
// environments: see Detector.h). Doubles are written with 17 significant
// digits, so round trips are exact.
#ifndef KG_SERIALIZATION_H
#define KG_SERIALIZATION_H

#include <string>

#include "kg/Event.h"

namespace kg {

std::string truthEventToYAML(const TruthEvent &ev);
TruthEvent truthEventFromYAML(const std::string &yaml);

std::string observedEventToYAML(const ObservedEvent &ev);
ObservedEvent observedEventFromYAML(const std::string &yaml);

} // namespace kg

#endif
