// Compatibility adapter for the original `process.card` format used by test.cc:
//
//   nevents <int>
//   outfile <root_file_name>
//   production <particleA> <particleB>
//   particle <name> <mass_min> <mass_max> <pid> <stable(0|1)> <visible(0|1)> <decay1> <decay2>
//
// The card is translated into an equivalent YAML prior that reproduces the
// legacy kinematics: exponential parton luminosity with tau = 0.05 + 0.02 *
// threshold/TeV and a 1 GeV threshold offset, soft-radiation recoil
// pT = 10 + Exp(5) GeV, isotropic decays, and the legacy rule that a decay
// whose daughters are too heavy becomes a three-body decay in which the
// heavier daughter is replaced by its own daughters (myGenerator's sequential
// three-body model). Mass ranges are re-sampled per event, as in test.cc,
// unless `task_level_masses` is set.
#ifndef KG_LEGACY_H
#define KG_LEGACY_H

#include <string>

#include "kg/Prior.h"

namespace kg {

struct LegacyOptions {
  bool task_level_masses = false; // true: a mass range becomes a task-level prior instead of per-event sampling
};

struct LegacyCard {
  int nevents = 1000000;
  std::string outfile = "test.root";
  std::string yaml; // equivalent prior
};

LegacyCard convertLegacyCard(const std::string &path, const LegacyOptions &opt = {});
LegacyCard convertLegacyCardString(const std::string &text, const LegacyOptions &opt = {});
PriorConfig priorFromLegacyCard(const std::string &path, const LegacyOptions &opt = {});

} // namespace kg

#endif
