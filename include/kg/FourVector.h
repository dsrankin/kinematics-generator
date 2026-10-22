// Minimal Lorentz-vector algebra (no ROOT dependency).
#ifndef KG_FOURVECTOR_H
#define KG_FOURVECTOR_H

#include <cmath>
#include <limits>

namespace kg {

struct Vec3 {
  double x = 0, y = 0, z = 0;
  Vec3() = default;
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  Vec3 operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  double dot(const Vec3 &o) const { return x * o.x + y * o.y + z * o.z; }
  Vec3 cross(const Vec3 &o) const { return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x}; }
  double mag2() const { return dot(*this); }
  double mag() const { return std::sqrt(mag2()); }
  Vec3 unit() const {
    const double m = mag();
    return m > 0 ? Vec3(x / m, y / m, z / m) : Vec3(0, 0, 1);
  }
};

struct P4 {
  double px = 0, py = 0, pz = 0, e = 0;

  P4() = default;
  P4(double px_, double py_, double pz_, double e_) : px(px_), py(py_), pz(pz_), e(e_) {}

  static P4 fromVectM(const Vec3 &p, double m) { return {p.x, p.y, p.z, std::sqrt(p.mag2() + m * m)}; }
  static P4 fromPtEtaPhiM(double pt, double eta, double phi, double m) {
    const double pz = pt * std::sinh(eta);
    const double px = pt * std::cos(phi), py = pt * std::sin(phi);
    return {px, py, pz, std::sqrt(px * px + py * py + pz * pz + m * m)};
  }
  // System of invariant mass m, rapidity y and transverse momentum pt at azimuth phi.
  static P4 fromMassRapidityPt(double m, double y, double pt, double phi) {
    const double mt = std::sqrt(m * m + pt * pt);
    return {pt * std::cos(phi), pt * std::sin(phi), mt * std::sinh(y), mt * std::cosh(y)};
  }

  P4 operator+(const P4 &o) const { return {px + o.px, py + o.py, pz + o.pz, e + o.e}; }
  P4 operator-(const P4 &o) const { return {px - o.px, py - o.py, pz - o.pz, e - o.e}; }
  P4 &operator+=(const P4 &o) {
    px += o.px; py += o.py; pz += o.pz; e += o.e;
    return *this;
  }
  P4 operator*(double s) const { return {px * s, py * s, pz * s, e * s}; }

  Vec3 vect() const { return {px, py, pz}; }
  double p2() const { return px * px + py * py + pz * pz; }
  double p() const { return std::sqrt(p2()); }
  double pt() const { return std::hypot(px, py); }
  double m2() const { return e * e - p2(); }
  double m() const {
    const double mm = m2();
    return mm >= 0 ? std::sqrt(mm) : -std::sqrt(-mm);
  }
  double phi() const { return (px == 0 && py == 0) ? 0.0 : std::atan2(py, px); }
  double eta() const {
    const double t = pt();
    if (t == 0) return pz == 0 ? 0.0 : (pz > 0 ? 1e10 : -1e10);
    return std::asinh(pz / t);
  }
  double rapidity() const { return 0.5 * std::log((e + pz) / (e - pz)); }
  Vec3 boostVector() const { return e != 0 ? Vec3(px / e, py / e, pz / e) : Vec3(); }

  // Active boost by velocity b (same convention as TLorentzVector::Boost).
  void boost(const Vec3 &b) {
    const double b2 = b.mag2();
    if (b2 <= 0) return;
    const double gamma = 1.0 / std::sqrt(1.0 - b2);
    const double bp = b.x * px + b.y * py + b.z * pz;
    const double gamma2 = (gamma - 1.0) / b2;
    px += gamma2 * bp * b.x + gamma * b.x * e;
    py += gamma2 * bp * b.y + gamma * b.y * e;
    pz += gamma2 * bp * b.z + gamma * b.z * e;
    e = gamma * (e + bp);
  }
};

inline double deltaPhi(double a, double b) {
  double d = a - b;
  while (d > M_PI) d -= 2 * M_PI;
  while (d <= -M_PI) d += 2 * M_PI;
  return d;
}

inline double wrapPhi(double phi) { return deltaPhi(phi, 0.0); }

} // namespace kg

#endif
