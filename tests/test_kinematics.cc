// Two- and three-body decay kinematics and angular models.
#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

namespace {

P4 randomParent(Rng &rng, double m) {
  return P4::fromVectM({rng.uniform(-500, 500), rng.uniform(-500, 500), rng.uniform(-2000, 2000)}, m);
}

void checkConserved(const P4 &a, const P4 &b) {
  const double scale = std::max(1.0, a.e);
  CHECK_CLOSE(a.px, b.px, 1e-9 * scale);
  CHECK_CLOSE(a.py, b.py, 1e-9 * scale);
  CHECK_CLOSE(a.pz, b.pz, 1e-9 * scale);
  CHECK_CLOSE(a.e, b.e, 1e-9 * scale);
}

// cos(theta*) of `d` in the rest frame of `parent`, relative to the parent direction.
double helicityCos(const P4 &parent, P4 d) {
  const Vec3 b = parent.boostVector();
  d.boost(-b);
  return d.vect().unit().dot(parent.vect().unit());
}

} // namespace

TEST(two_body_decay_conserves_four_momentum) {
  Rng rng(1);
  for (int i = 0; i < 2000; ++i) {
    const double M = rng.uniform(100, 1000);
    const double m[2] = {rng.uniform(0, 0.45 * M), rng.uniform(0, 0.45 * M)};
    const P4 parent = randomParent(rng, M);
    P4 out[2];
    decayTwoBody(parent, m, AngularModel::isotropic(), rng, out);
    checkConserved(parent, out[0] + out[1]);
    CHECK_CLOSE(out[0].m(), m[0], 1e-6 * M);
    CHECK_CLOSE(out[1].m(), m[1], 1e-6 * M);
    // Rest-frame momentum equals the two-body formula.
    P4 d = out[0];
    d.boost(-parent.boostVector());
    CHECK_CLOSE(d.p(), twoBodyMomentum(M, m[0], m[1]), 1e-6 * M);
  }
}

TEST(two_body_below_threshold_throws) {
  Rng rng(1);
  const double m[2] = {60, 60};
  P4 out[2];
  CHECK_THROWS_AS(decayTwoBody(P4(0, 0, 0, 100), m, AngularModel::isotropic(), rng, out), std::runtime_error, "threshold");
}

TEST(three_body_phase_space_is_flat_dalitz) {
  Rng rng(2);
  const double M = 100;
  const double m[3] = {0, 0, 0};
  int below = 0;
  Stats s01;
  const int n = 40000;
  for (int i = 0; i < n; ++i) {
    const P4 parent = randomParent(rng, M);
    P4 out[3];
    decayThreeBody(parent, m, ThreeBodyModel::PhaseSpace, AngularModel::isotropic(), rng, out);
    checkConserved(parent, out[0] + out[1] + out[2]);
    const double s = (out[0] + out[1]).m2();
    s01.add(s / (M * M));
    below += s < 0.5 * M * M;
    // Dalitz identity for massless daughters.
    CHECK_CLOSE(s + (out[1] + out[2]).m2() + (out[0] + out[2]).m2(), M * M, 1e-6 * M * M);
  }
  // Flat Dalitz density for massless daughters: s01/M^2 has density 2(1-x).
  CHECK_CLOSE(s01.mean(), 1.0 / 3.0, 4 * s01.err());
  CHECK_CLOSE(double(below) / n, 0.75, 4 * std::sqrt(0.75 * 0.25 / n));
}

TEST(three_body_massive_and_legacy_models_conserve) {
  Rng rng(3);
  for (auto model : {ThreeBodyModel::PhaseSpace, ThreeBodyModel::LegacySequential}) {
    Stats mY;
    for (int i = 0; i < 5000; ++i) {
      const double M = 500;
      const double m[3] = {80.4, 4.7, 172.5};
      const P4 parent = randomParent(rng, M);
      P4 out[3];
      decayThreeBody(parent, m, model, AngularModel::isotropic(), rng, out);
      checkConserved(parent, out[0] + out[1] + out[2]);
      for (int k = 0; k < 3; ++k) CHECK_CLOSE(out[k].m(), m[k], 1e-5 * M);
      mY.add((out[1] + out[2]).m());
    }
    if (model == ThreeBodyModel::LegacySequential) // m(12) uniform in [m1+m2, M-m0]
      CHECK_CLOSE(mY.mean(), 0.5 * (4.7 + 172.5 + 500 - 80.4), 4 * mY.err());
  }
}

TEST(three_body_at_threshold) {
  Rng rng(4);
  const double m[3] = {10, 20, 30};
  P4 out[3];
  decayThreeBody(P4(0, 0, 0, 60), m, ThreeBodyModel::PhaseSpace, AngularModel::isotropic(), rng, out);
  for (int k = 0; k < 3; ++k) CHECK_CLOSE(out[k].p(), 0, 1e-6);
}

TEST(polynomial_angular_model_in_helicity_frame) {
  Rng rng(5);
  const AngularModel ang = AngularModel::polynomial({1, 2, 1}); // (1+c)^2, <c> = 1/2
  Stats c;
  for (int i = 0; i < 40000; ++i) {
    const P4 parent = randomParent(rng, 80.4);
    const double m[2] = {0, 0};
    P4 out[2];
    decayTwoBody(parent, m, ang, rng, out);
    c.add(helicityCos(parent, out[0]));
  }
  CHECK_CLOSE(c.mean(), 0.5, 4 * c.err());
}

TEST(beam_axis_angular_model) {
  Rng rng(6);
  const AngularModel ang = AngularModel::polynomial({1, 0, 1}, AngularModel::Axis::Beam); // 1 + c^2: <c^2> = 2/5
  Stats c2;
  for (int i = 0; i < 40000; ++i) {
    const double m[2] = {0, 0};
    P4 out[2];
    decayTwoBody(P4(0, 0, 0, 1000), m, ang, rng, out);
    const double c = out[0].pz / out[0].p();
    c2.add(c * c);
  }
  CHECK_CLOSE(c2.mean(), 0.4, 4 * c2.err());
}

TEST(three_body_angular_model_governs_first_daughter) {
  Rng rng(7);
  const AngularModel ang = AngularModel::polynomial({1, 1}, AngularModel::Axis::Beam); // <c> = 1/3
  Stats c;
  for (int i = 0; i < 40000; ++i) {
    const double m[3] = {0, 0, 0};
    P4 out[3];
    decayThreeBody(P4(0, 0, 0, 100), m, ThreeBodyModel::PhaseSpace, ang, rng, out);
    c.add(out[0].pz / out[0].p());
  }
  CHECK_CLOSE(c.mean(), 1.0 / 3.0, 4 * c.err());
}

TEST(negative_polynomial_rejected) {
  CHECK_THROWS_AS(AngularModel::polynomial({0, 1}), std::invalid_argument, "negative");
  CHECK_THROWS_AS(AngularModel::polynomial({0, 0}), std::invalid_argument, "vanishes");
}

TEST(four_vector_helpers) {
  const P4 p = P4::fromPtEtaPhiM(50, 1.2, -2.0, 10);
  CHECK_CLOSE(p.pt(), 50, 1e-9);
  CHECK_CLOSE(p.eta(), 1.2, 1e-9);
  CHECK_CLOSE(p.phi(), -2.0, 1e-9);
  CHECK_CLOSE(p.m(), 10, 1e-9);
  const P4 s = P4::fromMassRapidityPt(300, 0.7, 25, 1.0);
  CHECK_CLOSE(s.m(), 300, 1e-9);
  CHECK_CLOSE(s.rapidity(), 0.7, 1e-12);
  CHECK_CLOSE(s.pt(), 25, 1e-12);
}
