// Truth and observed event records.
#ifndef KG_EVENT_H
#define KG_EVENT_H

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "kg/FourVector.h"

namespace kg {

enum class TruthStatus {
  Incoming = 0, // initial-state partons
  Recoil = 1,   // soft/ISR recoil against the hard system (unclustered unless emitted)
  Decayed = 2,  // intermediate particle
  Final = 3     // stable final-state particle
};

struct TruthParticle {
  int index = -1;
  std::string name;
  int pdg = 0;
  double charge = 0;
  P4 p4;
  double mass = 0;  // generated (possibly off-shell) mass
  int parent = -1;  // -1: produced in the hard process (or initial state)
  std::vector<int> children;
  TruthStatus status = TruthStatus::Final;
  int species = -1; // index into TaskConfig::particles (-1 for partons / recoil)
  int branch = -1;  // index of the produced particle this descends from
  int depth = 0;    // 0 for produced particles
  int decay_mode = -1;
  std::string object_type; // reconstructed type if visible, else "invisible"
  std::string reconstruct_as; // set when decayed through a composite mode: visible descendants form one object

  bool isFinal() const { return status == TruthStatus::Final; }
  bool visible() const { return isFinal() && object_type != "invisible"; }
};

struct TruthEvent {
  std::string task_id;
  uint64_t task_seed = 0;
  uint64_t event_index = 0;
  uint64_t seed = 0; // event seed (derived from task seed and event index)
  int component = 0; // index of the task component (process) this event belongs to
  std::string process;
  double sqrt_s_hat = 0;
  double x1 = 0, x2 = 0;
  std::vector<TruthParticle> particles;
  std::vector<std::pair<std::string, double>> event_params; // event-scoped parameter draws

  std::vector<int> finalState() const;
  std::vector<int> producedParticles() const;
  std::vector<int> ancestors(int i) const; // parent, grandparent, ...
  bool isDescendant(int i, int ancestor) const;
  P4 incomingSum() const;
  P4 finalSum() const; // final-state particles + recoil
  // Sum of invisible final-state transverse momentum.
  std::pair<double, double> trueMET() const;
  // Realised decay topology, e.g. "t>(W>(l,v),b) + t>(W>(q,q),b)".
  std::string topology() const;
};

struct ObservedObject {
  P4 p4;
  std::string type;
  int type_code = -1; // index into TaskConfig::object_types
  double charge = 0;
  std::vector<double> features; // aligned to ObservedEvent::feature_names; NaN = missing
  int truth_index = -1;          // TruthEvent particle index; -1 for fakes
  std::string true_type;         // truth object type (or "fake")
  bool is_fake = false;
  int n_constituents = 1;        // visible truth particles merged into this object (composites > 1 possible)
  int n_charged = 0;             // of which charged (e.g. tau prongs)

  double pt() const { return p4.pt(); }
  double eta() const { return p4.eta(); }
  double phi() const { return p4.phi(); }
  double m() const { return p4.m(); }
};

struct ObservedEvent {
  std::string task_id, env_id;
  uint64_t event_index = 0;
  uint64_t truth_seed = 0, detector_seed = 0;
  int component = 0;
  std::string process;
  std::vector<ObservedObject> objects;
  double met_x = 0, met_y = 0;
  double sum_et = 0;
  int n_lost = 0, n_fakes = 0, n_misid = 0;
  std::shared_ptr<const std::vector<std::string>> feature_names;
  std::vector<std::pair<std::string, double>> detector_params; // event-scoped detector draws

  double met() const { return std::hypot(met_x, met_y); }
  double metPhi() const { return (met_x == 0 && met_y == 0) ? 0.0 : std::atan2(met_y, met_x); }
  int featureIndex(const std::string &name) const;
  double feature(size_t object, const std::string &name) const;
};

} // namespace kg

#endif
