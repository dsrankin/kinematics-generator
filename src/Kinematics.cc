#include "kg/Kinematics.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace kg {

void AngularModel::prepare() {
  if (type == Type::Isotropic) {
    fmax = 1.0;
    return;
  }
  if (coeffs.empty()) throw std::invalid_argument("polynomial angular model needs at least one coefficient");
  double mx = 0, mn = 1e300;
  const int n = 4000;
  for (int i = 0; i <= n; ++i) {
    const double f = density(-1.0 + 2.0 * i / n);
    mx = std::max(mx, f);
    mn = std::min(mn, f);
  }
  if (mn < -1e-12) throw std::invalid_argument("polynomial angular model " + describe() + " is negative on [-1, 1]");
  if (!(mx > 0)) throw std::invalid_argument("polynomial angular model " + describe() + " vanishes on [-1, 1]");
  // Margin for maxima between grid points (|f'| bounded by the coefficients).
  double dsum = 0;
  for (size_t k = 1; k < coeffs.size(); ++k) dsum += k * std::fabs(coeffs[k]);
  fmax = mx + dsum * (2.0 / n);
}

double AngularModel::density(double c) const {
  if (type == Type::Isotropic) return 1.0;
  double f = 0, p = 1;
  for (double a : coeffs) {
    f += a * p;
    p *= c;
  }
  return f;
}

double AngularModel::sampleCosTheta(Rng &rng) const {
  if (type == Type::Isotropic) return rng.uniform(-1, 1);
  for (int it = 0; it < 1000000; ++it) {
    const double c = rng.uniform(-1, 1);
    if (rng.uniform() * fmax < density(c)) return c;
  }
  throw std::runtime_error("AngularModel: rejection sampling failed for " + describe());
}

std::string AngularModel::describe() const {
  std::ostringstream os;
  if (type == Type::Isotropic) {
    os << "isotropic";
  } else {
    os << "polynomial[";
    for (size_t i = 0; i < coeffs.size(); ++i) os << (i ? "," : "") << coeffs[i];
    os << "]";
  }
  os << (axis == Axis::Beam ? "@beam" : "@helicity");
  return os.str();
}

double twoBodyMomentum(double M, double m1, double m2) {
  const double tol = 1e-9 * std::max(1.0, M);
  if (M + tol < m1 + m2) {
    std::ostringstream os;
    os << "two-body decay below threshold: M=" << M << " < m1+m2=" << m1 + m2;
    throw std::runtime_error(os.str());
  }
  const double a = M * M - (m1 + m2) * (m1 + m2);
  const double b = M * M - (m1 - m2) * (m1 - m2);
  const double q2 = a * b;
  return q2 > 0 ? std::sqrt(q2) / (2 * M) : 0.0;
}

namespace {

// Orthonormal (e1, e2) perpendicular to n.
void perpBasis(const Vec3 &n, Vec3 &e1, Vec3 &e2) {
  const Vec3 t = std::fabs(n.z) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
  e1 = t.cross(n).unit();
  e2 = n.cross(e1);
}

Vec3 restAxis(const P4 &parent, AngularModel::Axis a) {
  if (a == AngularModel::Axis::Beam) return {0, 0, 1};
  const Vec3 p = parent.vect();
  return p.mag2() > 1e-24 * parent.e * parent.e ? p.unit() : Vec3(0, 0, 1);
}

// Map a configuration built with d0 along +z onto direction n, rolled by psi.
Vec3 orient(const Vec3 &v, const Vec3 &n, double psi) {
  Vec3 e1, e2;
  perpBasis(n, e1, e2);
  const double c = std::cos(psi), s = std::sin(psi);
  const Vec3 f1 = e1 * c + e2 * s;
  const Vec3 f2 = e2 * c - e1 * s;
  return f1 * v.x + f2 * v.y + n * v.z;
}

} // namespace

Vec3 sampleDirection(const AngularModel &model, const Vec3 &axis, Rng &rng) {
  const double c = model.sampleCosTheta(rng);
  const double s = std::sqrt(std::max(0.0, 1 - c * c));
  const double phi = rng.uniform(0, 2 * M_PI);
  Vec3 e1, e2;
  perpBasis(axis, e1, e2);
  return e1 * (s * std::cos(phi)) + e2 * (s * std::sin(phi)) + axis * c;
}

void decayTwoBody(const P4 &parent, const double m[2], const AngularModel &ang, Rng &rng, P4 out[2]) {
  const double M = parent.m();
  const double q = twoBodyMomentum(M, m[0], m[1]);
  const Vec3 n = sampleDirection(ang, restAxis(parent, ang.axis), rng);
  out[0] = P4::fromVectM(n * q, m[0]);
  out[1] = P4::fromVectM(n * (-q), m[1]);
  const Vec3 b = parent.boostVector();
  out[0].boost(b);
  out[1].boost(b);
}

void decayThreeBody(const P4 &parent, const double m[3], ThreeBodyModel model, const AngularModel &ang, Rng &rng,
                    P4 out[3]) {
  const double M = parent.m();
  const double tol = 1e-9 * std::max(1.0, M);
  if (M + tol < m[0] + m[1] + m[2]) {
    std::ostringstream os;
    os << "three-body decay below threshold: M=" << M << " < " << m[0] + m[1] + m[2];
    throw std::runtime_error(os.str());
  }
  Vec3 v[3]; // rest-frame momenta with d0 along +z
  if (model == ThreeBodyModel::PhaseSpace) {
    // Uniform in the Dalitz plane (s01, s12).
    const double lo01 = std::pow(m[0] + m[1], 2), hi01 = std::pow(std::max(M - m[2], m[0] + m[1]), 2);
    const double lo12 = std::pow(m[1] + m[2], 2), hi12 = std::pow(std::max(M - m[0], m[1] + m[2]), 2);
    double s01 = lo01, s12 = lo12;
    bool ok = hi01 <= lo01 || hi12 <= lo12; // degenerate: at threshold
    for (int it = 0; it < 100000 && !ok; ++it) {
      s01 = rng.uniform(lo01, hi01);
      s12 = rng.uniform(lo12, hi12);
      const double r = std::sqrt(s01);
      const double e1 = (s01 - m[0] * m[0] + m[1] * m[1]) / (2 * r); // d1 in (01) frame
      const double e2 = (M * M - s01 - m[2] * m[2]) / (2 * r);        // d2 in (01) frame
      const double p1 = std::sqrt(std::max(0.0, e1 * e1 - m[1] * m[1]));
      const double p2 = std::sqrt(std::max(0.0, e2 * e2 - m[2] * m[2]));
      const double smin = std::pow(e1 + e2, 2) - std::pow(p1 + p2, 2);
      const double smax = std::pow(e1 + e2, 2) - std::pow(p1 - p2, 2);
      ok = s12 >= smin && s12 <= smax;
    }
    if (!ok) throw std::runtime_error("three-body phase-space sampling failed");
    const double E0 = (M * M + m[0] * m[0] - s12) / (2 * M);
    const double E2 = (M * M + m[2] * m[2] - s01) / (2 * M);
    const double E1 = M - E0 - E2;
    const double q0 = std::sqrt(std::max(0.0, E0 * E0 - m[0] * m[0]));
    const double q1 = std::sqrt(std::max(0.0, E1 * E1 - m[1] * m[1]));
    const double q2 = std::sqrt(std::max(0.0, E2 * E2 - m[2] * m[2]));
    double c02 = (q0 > 0 && q2 > 0) ? (q1 * q1 - q0 * q0 - q2 * q2) / (2 * q0 * q2) : 0.0;
    c02 = std::min(1.0, std::max(-1.0, c02));
    const double s02 = std::sqrt(1 - c02 * c02);
    v[0] = {0, 0, q0};
    v[2] = {q2 * s02, 0, q2 * c02};
    v[1] = -(v[0] + v[2]);
  } else {
    const double mY = rng.uniform(m[1] + m[2], std::max(M - m[0], m[1] + m[2]));
    const double q = twoBodyMomentum(M, m[0], mY);
    v[0] = {0, 0, q};
    P4 pY = P4::fromVectM({0, 0, -q}, mY);
    const double qY = twoBodyMomentum(mY, m[1], m[2]);
    const Vec3 n = sampleDirection(AngularModel::isotropic(), {0, 0, 1}, rng);
    P4 a = P4::fromVectM(n * qY, m[1]), b = P4::fromVectM(n * (-qY), m[2]);
    a.boost(pY.boostVector());
    b.boost(pY.boostVector());
    v[1] = a.vect();
    v[2] = b.vect();
  }
  const Vec3 n = sampleDirection(ang, restAxis(parent, ang.axis), rng);
  const double psi = rng.uniform(0, 2 * M_PI);
  const Vec3 b = parent.boostVector();
  for (int i = 0; i < 3; ++i) {
    out[i] = P4::fromVectM(orient(v[i], n, psi), m[i]);
    out[i].boost(b);
  }
}

} // namespace kg
