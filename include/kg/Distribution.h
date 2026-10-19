// One-dimensional probability distributions used by priors.
#ifndef KG_DISTRIBUTION_H
#define KG_DISTRIBUTION_H

#include <limits>
#include <string>
#include <vector>

#include "kg/Random.h"

namespace YAML {
class Node;
class Emitter;
} // namespace YAML

namespace kg {

// Relativistic Breit-Wigner in m (density ~ 1/((m^2-M^2)^2 + M^2 G^2) in s=m^2),
// truncated exactly to [lo, hi] by inverse-CDF sampling.
double sampleRelativisticBW(Rng &rng, double mean, double width, double lo, double hi);

class Distribution {
 public:
  enum class Type { Fixed, Uniform, LogUniform, Gaussian, TruncatedGaussian, BreitWigner, Categorical, Exponential, Poisson };

  Type type = Type::Fixed;
  double value = 0;                                          // fixed
  double min = -std::numeric_limits<double>::infinity();     // uniform, loguniform, truncated_gaussian, breit_wigner
  double max = std::numeric_limits<double>::infinity();
  double mean = 0, sigma = 0;                                // gaussian, truncated_gaussian; mean also for breit_wigner
  double width = 0;                                          // breit_wigner (FWHM)
  bool relativistic = false;                                 // breit_wigner
  double scale = 1, offset = 0;                              // exponential: offset + Exp(scale)
                                                             // poisson: uses mean
  std::vector<double> values, weights;                       // categorical

  static Distribution fixed(double v);
  static Distribution uniform(double lo, double hi);
  static Distribution logUniform(double lo, double hi);
  static Distribution gaussian(double mu, double sd);
  static Distribution truncatedGaussian(double mu, double sd, double lo, double hi);
  static Distribution breitWigner(double mu, double gamma, double lo = -std::numeric_limits<double>::infinity(),
                                  double hi = std::numeric_limits<double>::infinity(), bool relativistic = false);
  static Distribution categorical(std::vector<double> vals, std::vector<double> w = {});
  static Distribution exponential(double scale, double offset = 0);
  static Distribution poisson(double mean);

  double sample(Rng &rng) const;

  // Throws std::invalid_argument with `context` in the message if parameters are inconsistent.
  void validate(const std::string &context = "") const;

  std::string typeName() const;

  // Parse from a YAML map containing `dist: <type>` and the type's fields.
  static Distribution fromYAML(const YAML::Node &node, const std::string &context = "");
  // Emit the `dist` key and fields into an already-open YAML map.
  void emitFields(YAML::Emitter &out) const;
};

} // namespace kg

#endif
