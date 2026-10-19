// Decay kinematics: angular models, two- and three-body decays, production.
#ifndef KG_KINEMATICS_H
#define KG_KINEMATICS_H

#include <string>
#include <vector>

#include "kg/FourVector.h"
#include "kg/Random.h"

namespace kg {

// Distribution of cos(theta*) of the first daughter in the parent rest frame,
// measured with respect to `axis`:
//   helicity: the parent's flight direction in the lab,
//   beam:     the z axis (after a pure boost into the rest frame).
// The azimuth around the axis is always uniform.
struct AngularModel {
  enum class Type { Isotropic, Polynomial };
  enum class Axis { Helicity, Beam };

  Type type = Type::Isotropic;
  Axis axis = Axis::Helicity;
  std::vector<double> coeffs; // f(c) = sum_i coeffs[i] * c^i, must be >= 0 on [-1, 1]

  static AngularModel isotropic(Axis a = Axis::Helicity) {
    AngularModel m;
    m.axis = a;
    return m;
  }
  static AngularModel polynomial(std::vector<double> c, Axis a = Axis::Helicity) {
    AngularModel m;
    m.type = Type::Polynomial;
    m.coeffs = std::move(c);
    m.axis = a;
    m.prepare();
    return m;
  }

  // Validate (non-negative, non-zero on [-1,1]) and compute the rejection bound.
  void prepare();
  double density(double c) const;
  double sampleCosTheta(Rng &rng) const;
  std::string describe() const;

  double fmax = 1.0;
};

enum class ThreeBodyModel {
  PhaseSpace,       // uniform Dalitz density (flat three-body phase space)
  LegacySequential  // 0 -> d0 + Y, Y -> d1 d2 with m_Y uniform (myGenerator::getRand3BodyDecay)
};

// |p| of either daughter in the rest frame of a parent of mass M. Throws below threshold.
double twoBodyMomentum(double M, double m1, double m2);

// Unit vector with cos(theta) ~ model about `axis`, uniform azimuth.
Vec3 sampleDirection(const AngularModel &model, const Vec3 &axis, Rng &rng);

// parent -> d[0] d[1]; masses m[0], m[1]. Output four-vectors are in the lab frame.
void decayTwoBody(const P4 &parent, const double m[2], const AngularModel &ang, Rng &rng, P4 out[2]);

// parent -> d[0] d[1] d[2]. The angular model governs the direction of d[0];
// the configuration is rolled uniformly about that direction.
void decayThreeBody(const P4 &parent, const double m[3], ThreeBodyModel model, const AngularModel &ang, Rng &rng,
                    P4 out[3]);

} // namespace kg

#endif
