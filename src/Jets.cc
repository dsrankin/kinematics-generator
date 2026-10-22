#include "kg/Jets.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "kg/Task.h"

namespace kg {

namespace {

const double kCF = 4.0 / 3.0, kCA = 3.0;
const double kMZ = 91.1876;
const double kB0 = (33.0 - 2.0 * 5) / (12.0 * M_PI);
const double kInf = std::numeric_limits<double>::infinity();

// Unit vector at polar angle `alpha` around axis n (azimuth psi in the frame e1, e2).
Vec3 rotated(const Vec3 &n, const Vec3 &e1, const Vec3 &e2, double alpha, double psi) {
  const double s = std::sin(alpha);
  return n * std::cos(alpha) + e1 * (s * std::cos(psi)) + e2 * (s * std::sin(psi));
}

Vec3 cross(const Vec3 &a, const Vec3 &b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }

void frame(const Vec3 &n, Vec3 &e1, Vec3 &e2) {
  const Vec3 ref = std::fabs(n.x) < 0.6 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
  e1 = cross(n, ref);
  e1 = e1 * (1.0 / e1.mag());
  e2 = cross(n, e1);
}

P4 withEnergy(const Vec3 &dir, double e, double m) {
  const double p = e > m ? std::sqrt(e * e - m * m) : 0.0;
  return {dir.x * p, dir.y * p, dir.z * p, e};
}

double openingAngle(const P4 &a, const P4 &b) {
  const double pa = a.p(), pb = b.p();
  if (pa <= 0 || pb <= 0) return 0;
  return std::acos(std::max(-1.0, std::min(1.0, (a.px * b.px + a.py * b.py + a.pz * b.pz) / (pa * pb))));
}

// Quarks, gluons and their SUSY partners (PDG numbering).
bool coloured(int pdg) {
  const int a = std::abs(pdg), r = a % 1000000;
  return a <= 8 || a == 21 || (a > 1000000 && a < 3000000 && (r <= 6 || r == 21));
}

} // namespace

double alphaS(double alpha_s_mz, double q) {
  const double d = 1.0 + alpha_s_mz * kB0 * std::log(q * q / (kMZ * kMZ));
  return alpha_s_mz / std::max(d, 0.1);
}

bool showersParton(const TruthParticle &p) {
  const int a = std::abs(p.pdg);
  return p.visible() && ((a >= 1 && a <= 5) || a == 21);
}

void showerParton(const P4 &p0, int pdg, double mass, double theta_max, const ShowerSettings &s, Rng &rng,
                  std::vector<ShowerParticle> &out, int parent) {
  // Leading-log, angular-ordered shower in the emission angle theta and the energy fraction z
  // (<= 1/2) of the softer daughter:
  //   dP = alpha_s(z E theta) / pi * dtheta / theta * P(z) dz,  z E theta > kt_min,
  // with P(z) = C_F (1 + (1 - z)^2) / z for q -> q g and 2 C_A [z / (1 - z) + (1 - z) / z + z (1 - z)]
  // for g -> g g (both daughters counted). Emissions are generated with the veto algorithm from the
  // overestimate alpha_s(kt_min) / pi * c / z, c = 2 C_F or 2.25 C_A. Massive quarks do not radiate
  // inside their dead cone theta < m / E. Each emission splits the energy exactly; daughters are
  // placed at angles z theta and (1 - z) theta on opposite sides of the parent direction (transverse
  // momentum balance), so the pair mass is ~ z (1 - z) E^2 theta^2, and continue at smaller angles.
  struct Item {
    P4 p;
    int pdg;
    double mass, theta_max;
  };
  std::vector<Item> stack{{p0, pdg, mass, theta_max}};
  const double asMax = alphaS(s.alpha_s, s.kt_min);
  int emissions = 0;
  while (!stack.empty()) {
    const Item it = stack.back();
    stack.pop_back();
    const bool gluon = it.pdg == 21;
    const double E = it.p.e, pm = it.p.p();
    const double over = gluon ? 2.25 * kCA : 2.0 * kCF;
    const double zlo = s.kt_min / (E * it.theta_max);
    bool split = false;
    double theta = 0, z = 0;
    if (zlo < 0.5 && pm > 0 && emissions < s.max_emissions) {
      const double L = std::log(0.5 / zlo);
      const double rate = asMax / M_PI * over * L; // per unit log(theta_max / theta)
      double t = 0;
      while (true) {
        t += rng.exponential(1.0) / rate;
        theta = it.theta_max * std::exp(-t);
        if (0.5 * E * theta < s.kt_min) break;
        if (it.mass > 0 && theta < it.mass / E) break; // dead cone: all smaller angles are vetoed too
        z = zlo * std::exp(rng.uniform() * L);
        const double kt = z * E * theta;
        const double u = rng.uniform();
        if (kt < s.kt_min) continue;
        const double zP = gluon ? 2.0 * kCA * (z * z / (1.0 - z) + (1.0 - z) + z * z * (1.0 - z)) : kCF * (1.0 + (1.0 - z) * (1.0 - z));
        if (u < alphaS(s.alpha_s, kt) / asMax * zP / over) {
          split = true;
          break;
        }
      }
    }
    if (!split) {
      P4 f = it.p;
      if (s.hadron_mass > it.mass && pm > 0) {
        const double m = std::min(s.hadron_mass, 0.99 * E);
        const double q = std::sqrt(std::max(0.0, E * E - m * m)) / pm;
        f = P4(f.px * q, f.py * q, f.pz * q, E);
      }
      out.push_back({parent, it.pdg, f});
      continue;
    }
    ++emissions;
    const Vec3 n(it.p.px / pm, it.p.py / pm, it.p.pz / pm);
    Vec3 e1, e2;
    frame(n, e1, e2);
    const double psi = rng.uniform(-M_PI, M_PI);
    const Vec3 dHard = rotated(n, e1, e2, z * theta, psi + M_PI);
    const Vec3 dSoft = rotated(n, e1, e2, (1.0 - z) * theta, psi);
    stack.push_back({withEnergy(dSoft, z * E, 0.0), 21, 0.0, theta});
    stack.push_back({withEnergy(dHard, (1.0 - z) * E, it.mass), it.pdg, it.mass, theta});
  }
}

ShowerSettings showerSettings(const TaskConfig &t, const std::vector<double> &v) {
  ShowerSettings s;
  if (!t.shower.enabled) return s;
  s.alpha_s = v[t.shower.alpha_s];
  s.kt_min = v[t.shower.kt_min];
  s.theta_max = v[t.shower.theta_max];
  if (t.shower.hadron_mass >= 0) s.hadron_mass = v[t.shower.hadron_mass];
  if (t.shower.hadron_shift >= 0) s.hadron_shift = v[t.shower.hadron_shift];
  return s;
}

P4 hadronisedJet(const P4 &jet, int pdg, double R, const ShowerSettings &s) {
  if (s.hadron_shift <= 0) return jet;
  const double pt = jet.pt();
  if (pt <= 0) return jet;
  const double c = std::abs(pdg) == 21 ? kCA / kCF : 1.0;
  const double m0 = std::max(0.0, jet.m());
  return P4::fromPtEtaPhiM(pt, jet.eta(), jet.phi(), std::sqrt(m0 * m0 + s.hadron_shift * c * pt * R));
}

void applyPartonShower(const TaskConfig &t, const std::vector<double> &v, TruthEvent &ev) {
  if (!t.shower.enabled) return;
  const ShowerSettings s = showerSettings(t, v);
  Rng rng(mixSeed(ev.seed, kStreamShower));
  for (const TruthParticle &p : ev.particles) {
    if (!showersParton(p)) continue;
    // Colour coherence: the two partons of a colour-singlet decay (W, Z, H -> q q) radiate only inside
    // their opening angle; other partons radiate up to theta_max.
    double thetaMax = s.theta_max;
    if (p.parent >= 0) {
      const TruthParticle &par = ev.particles[p.parent];
      if (par.children.size() == 2 && !coloured(par.pdg)) {
        const TruthParticle &a = ev.particles[par.children[0]], &b = ev.particles[par.children[1]];
        if (showersParton(a) && showersParton(b)) thetaMax = std::min(thetaMax, openingAngle(a.p4, b.p4));
      }
    }
    showerParton(p.p4, p.pdg, std::max(0.0, p.mass), thetaMax, s, rng, ev.shower, p.index);
  }
}

std::vector<std::vector<int>> antiKtCluster(const std::vector<P4> &in, double R) {
  const size_t n = in.size();
  std::vector<std::vector<int>> jets;
  std::vector<P4> mom(in);
  std::vector<double> y(n), phi(n), k(n), nnd(n, kInf);
  std::vector<int> nn(n, -1);
  std::vector<bool> active(n, true);
  std::vector<std::vector<int>> members(n);
  auto kin = [&](size_t i) {
    const P4 &p = mom[i];
    const double pt2 = p.px * p.px + p.py * p.py;
    k[i] = pt2 > 0 ? 1.0 / pt2 : kInf;
    phi[i] = p.phi();
    const double ep = p.e + p.pz, em = p.e - p.pz;
    y[i] = (ep > 0 && em > 0) ? 0.5 * std::log(ep / em) : (p.pz >= 0 ? 1e3 : -1e3);
  };
  auto d2 = [&](size_t i, size_t j) {
    const double dy = y[i] - y[j], dp = deltaPhi(phi[i], phi[j]);
    return dy * dy + dp * dp;
  };
  auto findNN = [&](size_t i) {
    nn[i] = -1;
    nnd[i] = kInf;
    for (size_t j = 0; j < n; ++j)
      if (j != i && active[j]) {
        const double d = d2(i, j);
        if (d < nnd[i]) {
          nnd[i] = d;
          nn[i] = static_cast<int>(j);
        }
      }
  };
  for (size_t i = 0; i < n; ++i) {
    kin(i);
    members[i] = {static_cast<int>(i)};
  }
  for (size_t i = 0; i < n; ++i) findNN(i);
  const double R2 = R * R;
  size_t left = n;
  while (left > 0) {
    // Smallest of d_iB = 1/pT_i^2 and d_ij = min(1/pT_i^2, 1/pT_j^2) dR_ij^2 / R^2; for the pair
    // minimising d_ij, j is the geometric nearest neighbour of i (Cacciari-Salam lemma).
    double best = kInf;
    int bi = -1, bj = -1;
    for (size_t i = 0; i < n; ++i) {
      if (!active[i]) continue;
      if (k[i] < best || bi < 0) {
        best = k[i];
        bi = static_cast<int>(i);
        bj = -1;
      }
      if (nn[i] >= 0) {
        const double dij = std::min(k[i], k[nn[i]]) * nnd[i] / R2;
        if (dij < best) {
          best = dij;
          bi = static_cast<int>(i);
          bj = nn[i];
        }
      }
    }
    if (bj < 0) {
      jets.push_back(members[bi]);
      active[bi] = false;
      --left;
      for (size_t m = 0; m < n; ++m)
        if (active[m] && nn[m] == bi) findNN(m);
      continue;
    }
    mom[bi] += mom[bj];
    members[bi].insert(members[bi].end(), members[bj].begin(), members[bj].end());
    active[bj] = false;
    --left;
    kin(bi);
    findNN(bi);
    for (size_t m = 0; m < n; ++m) {
      if (!active[m] || static_cast<int>(m) == bi) continue;
      if (nn[m] == bi || nn[m] == bj) findNN(m);
      else {
        const double d = d2(m, bi);
        if (d < nnd[m]) {
          nnd[m] = d;
          nn[m] = bi;
        }
      }
    }
  }
  return jets;
}

} // namespace kg
