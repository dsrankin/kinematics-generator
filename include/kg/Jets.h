// Jets: a final-state parton shower and anti-kT clustering.
//
// The parton shower (opt-in, `parton_shower` in the prior) turns every visible final-state
// quark or gluon of a truth event into a set of collinear partons (TruthEvent::shower), using a
// leading-logarithmic, angular-ordered shower with running coupling, a transverse-momentum
// cutoff and a dead cone for massive quarks. Clustering (opt-in, `detector.jets`) groups
// reconstruction candidates into jets of a given radius with the anti-kT algorithm, so that jet
// masses arise from the shower (QCD jets) and from merging the partons of boosted heavy
// particles (W, Z, H, top).
#ifndef KG_JETS_H
#define KG_JETS_H

#include <vector>

#include "kg/Event.h"
#include "kg/FourVector.h"
#include "kg/Random.h"

namespace kg {

struct TaskConfig;

// Strong coupling at scale q (GeV), one loop with five flavours from alpha_s(mZ).
double alphaS(double alpha_s_mz, double q);

// True for truth particles the shower acts on: visible final-state quarks (|pdg| 1-5) and gluons.
bool showersParton(const TruthParticle &p);

struct ShowerSettings {
  double alpha_s = 0.118;  // alpha_s(mZ)
  double kt_min = 1.0;     // GeV: emissions below this transverse momentum are not resolved
  double theta_max = 1.0;  // largest emission angle (radians) for partons without a colour partner
  double hadron_mass = 0.0; // GeV: mass given to every final shower parton (energy conserved)
  // Non-perturbative jet-mass shift added to every jet built from shower partons:
  // m^2 -> m^2 + hadron_shift * (C_i / C_F) * pT * R, with C_i of the jet's leading parton.
  double hadron_shift = 0.0; // GeV
  int max_emissions = 1000; // per hard parton (safety limit)
};

// Shower one parton (pdg 21 = gluon, else a quark of the given mass) emitting at angles below
// theta_max. Appends the resulting partons (energy is conserved; momenta are rescaled for masses).
void showerParton(const P4 &p, int pdg, double mass, double theta_max, const ShowerSettings &s, Rng &rng,
                  std::vector<ShowerParticle> &out, int parent);

// Shower all eligible partons of a truth event (fills ev.shower); uses the task's parton_shower
// parameters `values` and the event's own shower random stream. No-op when the shower is disabled.
void applyPartonShower(const TaskConfig &t, const std::vector<double> &values, TruthEvent &ev);

// Settings of a task's shower (parameter values `values`); defaults when the shower is disabled.
ShowerSettings showerSettings(const TaskConfig &t, const std::vector<double> &values);

// Jet mass after the non-perturbative shift for a jet of radius R led by a parton `pdg`.
P4 hadronisedJet(const P4 &jet, int pdg, double R, const ShowerSettings &s);

// Underlying event: Poisson(density * A / pt_mean) massless particles over |eta| < eta_max
// (area A = 2 eta_max * 2 pi), pT ~ Exp(pt_mean), flat in eta and phi, shifted to zero total pT
// (they recoil against each other). Fills ev.underlying_event from the event's own random stream.
void applyUnderlyingEvent(const TaskConfig &t, const std::vector<double> &values, TruthEvent &ev);

// Anti-kT clustering with radius R in (rapidity, phi). Returns the input indices of each jet.
std::vector<std::vector<int>> antiKtCluster(const std::vector<P4> &particles, double R);

} // namespace kg

#endif
