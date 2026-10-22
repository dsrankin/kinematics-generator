#ifndef KG_TEST_HELPERS_H
#define KG_TEST_HELPERS_H

#include <cmath>
#include <string>
#include <vector>

#include "kg/kinegen.h"

namespace kgtest {

inline kg::PriorConfig prior(const std::string &yaml) { return kg::parsePriorString(yaml); }

inline std::string examplePath(const std::string &name) { return std::string(KG_SOURCE_DIR) + "/examples/" + name; }

struct Stats {
  double n = 0, s = 0, s2 = 0;
  void add(double x) {
    n += 1;
    s += x;
    s2 += x * x;
  }
  double mean() const { return s / n; }
  double sd() const { return std::sqrt(std::max(0.0, s2 / n - mean() * mean())); }
  double err() const { return sd() / std::sqrt(n); }
};

// A perfect detector (all defaults) appended to a physics prior.
inline const char *kTwoLeptons = R"(
particles:
  Z: {pdg: 23, mass: 91.19, width: 2.5, decays: [{daughters: [mu+, mu-], br: 1}]}
  mu+: {pdg: -13, charge: 1}
  mu-: {pdg: 13, charge: -1}
production: [Z]
initial_state:
  sqrt_s: {model: exponential_pdf}
  recoil: {model: exponential, min_pt: 10, mean_pt: 30}
)";

} // namespace kgtest

#endif
