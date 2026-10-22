// Composite reconstruction (reconstruct_as): hadronic taus and merged objects.
#include <map>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;

namespace {

const char *kDitau = R"(
particles:
  Z: {pdg: 23, mass: 91.19, decays: [[tau+, tau-]]}
  tau+:
    pdg: -15
    charge: 1
    mass: 1.777
    decays:
      - {daughters: [e+, ve, vt], br: 0.2}
      - {daughters: [pi+, vt], br: 0.2, reconstruct_as: tau}
      - {daughters: [rho+, vt], br: 0.3, reconstruct_as: tau}
      - {daughters: [a1+, vt], br: 0.3, reconstruct_as: tau}
      - {daughters: [pi+, pi+, pi-, vt], br: 0.0, reconstruct_as: tau}
  tau-:
    pdg: 15
    charge: -1
    mass: 1.777
    decays:
      - {daughters: [e-, ve, vt], br: 0.2}
      - {daughters: [pi-, vt], br: 0.2, reconstruct_as: tau}
      - {daughters: [rho-, vt], br: 0.3, reconstruct_as: tau}
      - {daughters: [a1-, vt], br: 0.3, reconstruct_as: tau}
  rho+: {pdg: 213, charge: 1, mass: 0.775, width: 0.149, decays: [[pi+, pi0]]}
  rho-: {pdg: -213, charge: -1, mass: 0.775, width: 0.149, decays: [[pi-, pi0]]}
  a1-: {pdg: -20213, charge: -1, mass: 1.23, width: 0.4, decays: [[pi-, pi-, pi+]]}
  a1+: {pdg: 20213, charge: 1, mass: 1.23, width: 0.4, decays: [[pi+, pi+, pi-]]}
  pi0: {pdg: 111, mass: 0.135, decays: [[gamma, gamma]]}
  pi+: {pdg: 211, charge: 1, mass: 0.1396, object: jet}
  pi-: {pdg: -211, charge: -1, mass: 0.1396, object: jet}
  gamma: {pdg: 22}
  e+: {pdg: -11, charge: 1}
  e-: {pdg: 11, charge: -1}
  ve: {pdg: 12}
  vt: {pdg: 16}
production: [Z]
)";

} // namespace

TEST(composite_four_mode_validation) {
  // A four-body mode is still rejected even when composite.
  CHECK_THROWS_AS(kgtest::prior(kDitau), std::invalid_argument, "two- and three-body");
}

TEST(hadronic_taus_reconstructed_as_single_objects) {
  std::string y = kDitau;
  y.erase(y.find("      - {daughters: [pi+, pi+, pi-, vt], br: 0.0, reconstruct_as: tau}\n"), 71);
  const PriorConfig p = kgtest::prior(y);
  CHECK(p.objectTypes() == (std::vector<std::string>{"electron", "jet", "photon", "tau"}));
  const TaskConfig t = sampleTask(p, 1, 0);
  const DetectorModel det(t);
  std::map<int, int> prongs;
  int nTau = 0, nLep = 0;
  for (int i = 0; i < 4000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev);
    for (const auto &o : ob.objects) {
      CHECK(o.type != "jet" && o.type != "photon"); // pions and photons never appear on their own
      if (o.type == "tau") {
        ++nTau;
        const TruthParticle &tau = ev.particles[o.truth_index];
        CHECK(tau.name == "tau+" || tau.name == "tau-");
        CHECK(tau.reconstruct_as == "tau");
        // Perfect detector: the object is the visible part of the tau.
        P4 vis;
        int n = 0, q = 0;
        for (const auto &d : ev.particles)
          if (d.visible() && ev.isDescendant(d.index, tau.index)) {
            vis += d.p4;
            ++n;
            q += d.charge != 0;
          }
        CHECK_CLOSE(o.p4.pt(), vis.pt(), 1e-9 * vis.e);
        CHECK_CLOSE(o.p4.m(), vis.m(), 1e-6);
        CHECK(o.p4.m() < 1.777);
        CHECK(o.n_constituents == n && o.n_charged == q);
        CHECK_CLOSE(o.charge, tau.charge, 1e-12);
        prongs[o.n_charged]++;
      } else {
        CHECK(o.type == "electron" && o.n_constituents == 1 && o.n_charged == 1);
        ++nLep;
      }
    }
    // Neutrinos are the only missing momentum.
    const auto met = ev.trueMET();
    CHECK_CLOSE(ob.met_x, met.first, 1e-8 * ev.sqrt_s_hat);
  }
  CHECK(nTau > 4000 && nLep > 1000);
  CHECK(prongs[1] > 0 && prongs[3] > 0); // one- and three-prong taus
  // Composite types get their own detector response.
  const PriorConfig q = kgtest::prior(y + "detector: {objects: {tau: {efficiency: 0.5, misid: {jet: 0.2}}}}\n");
  const TaskConfig u = sampleTask(q, 1, 0);
  int tau = 0, asJet = 0, total = 0;
  for (int i = 0; i < 4000; ++i) {
    const TruthEvent ev = generateTruthEvent(u, i);
    for (const auto &pp : ev.particles) total += pp.reconstruct_as == "tau";
    for (const auto &o : DetectorModel(u).apply(ev).objects) {
      tau += o.type == "tau";
      asJet += o.type == "jet" && o.true_type == "tau";
    }
  }
  CHECK_CLOSE(double(tau + asJet) / total, 0.5, 0.03);
  CHECK_CLOSE(double(asJet) / (tau + asJet), 0.2, 0.03);
}

TEST(nested_composites_use_outermost) {
  const PriorConfig p = kgtest::prior(R"(
particles:
  W: {pdg: 24, mass: 80.4, decays: [{daughters: [q, X], br: 1, reconstruct_as: fatjet}]}
  X: {pdg: 9, mass: 20, decays: [{daughters: [q, q], br: 1, reconstruct_as: jet}]}
  q: {pdg: 1}
production: [W]
initial_state: {recoil: {model: exponential, min_pt: 300, mean_pt: 400}}
)");
  const TaskConfig t = sampleTask(p, 1, 0);
  for (int i = 0; i < 200; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = DetectorModel(t).apply(ev);
    CHECK(ob.objects.size() == 1);
    CHECK(ob.objects[0].type == "fatjet" && ob.objects[0].n_constituents == 3);
    CHECK_CLOSE(ob.objects[0].p4.m(), ev.particles[ob.objects[0].truth_index].mass, 1e-6);
  }
  // Serialisation keeps the composite information.
  const TaskConfig r = taskFromYAML(taskToYAML(t));
  CHECK(truthEventToYAML(generateTruthEvent(r, 3)) == truthEventToYAML(generateTruthEvent(t, 3)));
  const ObservedEvent ob = DetectorModel(t).apply(generateTruthEvent(t, 3));
  CHECK(observedEventToYAML(observedEventFromYAML(observedEventToYAML(ob))) == observedEventToYAML(ob));
}

TEST(composite_feature_quantities) {
  const PriorConfig p = kgtest::prior(R"(
particles:
  X: {pdg: 9, mass: 50, decays: [{daughters: [pi+, pi-, gamma], br: 1, reconstruct_as: tau}]}
  pi+: {pdg: 211, charge: 1, mass: 0.14, object: jet}
  pi-: {pdg: -211, charge: -1, mass: 0.14, object: jet}
  gamma: {pdg: 22}
production: [X, X]
features:
  - {name: prongs, kind: observed_proxy, quantity: n_charged}
  - {name: ncon, kind: observed_proxy, quantity: n_constituents}
)");
  const TaskConfig t = sampleTask(p, 1, 0);
  const ObservedEvent ob = DetectorModel(t).apply(generateTruthEvent(t, 0));
  CHECK(ob.objects.size() == 2);
  CHECK(ob.feature(0, "prongs") == 2 && ob.feature(0, "ncon") == 3);
}
