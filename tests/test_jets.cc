// Parton shower, anti-kT clustering and jet masses.
#include <algorithm>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;
using kgtest::Stats;

namespace {

P4 massless(double pt, double eta, double phi) { return P4::fromPtEtaPhiM(pt, eta, phi, 0.0); }

// Boosted heavy resonance -> W W (W -> q q) or a QCD dijet, with a perfect detector apart from
// the jet collections given in `jets`.
PriorConfig boosted(const std::string &process, const std::string &shower, const std::string &jets) {
  return kgtest::prior(R"(
particles:
  X: {pdg: 9000023, mass: 2500, width: 1, decays: [{daughters: [W, W], br: 1}]}
  W: {pdg: 24, mass: 80.4, width: 2.09, decays: [{daughters: [q, q], br: 1}]}
  q: {pdg: 1}
  g: {pdg: 21}
production: )" + process + R"(
initial_state: {sqrt_s: {model: fixed, value: 2500}, rapidity: {dist: uniform, min: -1, max: 1, scope: event}}
)" + shower + R"(
detector:
  objects:
    jet: {pt_min: 20}
    fatjet: {pt_min: 200}
)" + jets);
}

const char *kShower = "parton_shower: {}\n"; // calibrated defaults
const char *kJets = "  jets:\n    - {type: jet, radius: 0.4}\n    - {type: fatjet, radius: 0.8, in_met: false}\n";

const ObservedObject *leading(const ObservedEvent &ev, const std::string &type) {
  for (const auto &o : ev.objects)
    if (o.type == type) return &o; // objects are pT ordered
  return nullptr;
}

} // namespace

TEST(anti_kt_merges_within_radius) {
  const std::vector<P4> in = {massless(100, 0, 0), massless(50, 0.3, 0), massless(30, 1.0, 0), massless(5, 0.2, 0.25)};
  auto jets = antiKtCluster(in, 0.4);
  CHECK(jets.size() == 2);
  std::vector<int> hard;
  for (const auto &j : jets)
    if (std::find(j.begin(), j.end(), 0) != j.end()) hard = j;
  std::sort(hard.begin(), hard.end());
  CHECK((hard == std::vector<int>{0, 1, 3}));
  // Far apart, every particle is its own jet; with a huge radius everything merges.
  CHECK(antiKtCluster(in, 0.1).size() == 4);
  CHECK(antiKtCluster(in, 2.0).size() == 1);
  CHECK(antiKtCluster({}, 0.4).empty());
  // Every input ends up in exactly one jet.
  std::vector<P4> many;
  Rng rng(3);
  for (int i = 0; i < 300; ++i) many.push_back(massless(rng.exponential(20.0) + 1, rng.uniform(-3, 3), rng.uniform(-M_PI, M_PI)));
  std::vector<int> seen(many.size(), 0);
  for (const auto &j : antiKtCluster(many, 0.4))
    for (int k : j) ++seen[k];
  CHECK(std::all_of(seen.begin(), seen.end(), [](int s) { return s == 1; }));
}

TEST(shower_conserves_energy_and_gives_mass) {
  ShowerSettings s;
  Stats mq, mg, nq, ng;
  Rng rng(11);
  for (int i = 0; i < 4000; ++i) {
    for (int pdg : {1, 21}) {
      std::vector<ShowerParticle> out;
      const P4 p = massless(500, 0.5, 1.0);
      showerParton(p, pdg, 0.0, 0.4, s, rng, out, 7);
      P4 sum;
      for (const auto &q : out) {
        sum += q.p4;
        CHECK(q.parent == 7);
      }
      CHECK_CLOSE(sum.e, p.e, 1e-9 * p.e);
      CHECK(sum.m() > -1e-3 && sum.m() < 0.4 * p.e);
      (pdg == 21 ? mg : mq).add(sum.m() / sum.pt());
      (pdg == 21 ? ng : nq).add(out.size());
    }
  }
  // Gluons radiate more (C_A > C_F): heavier jets with more partons; m/pT of order 0.1 at R ~ 0.4.
  CHECK(mg.mean() > mq.mean() * 1.3);
  CHECK(ng.mean() > nq.mean());
  CHECK(mq.mean() > 0.03 && mg.mean() < 0.2);
  // A heavy quark does not radiate inside its dead cone.
  std::vector<ShowerParticle> out;
  showerParton(massless(10, 0, 0), 5, 4.7, 0.3, s, rng, out, 0);
  CHECK(out.size() == 1);
}

TEST(shower_is_opt_in_and_leaves_the_hard_event_unchanged) {
  const PriorConfig off = boosted("[X]", "", ""), noShower = boosted("[X]", "parton_shower: false\n", "");
  const PriorConfig on = boosted("[X]", "parton_shower: {hadron_shift: 0}\n", "");
  CHECK(!off.shower.enabled && !noShower.shower.enabled && on.shower.enabled);
  const TaskConfig a = sampleTask(off, 1, 0), b = sampleTask(noShower, 1, 0), c = sampleTask(on, 1, 0);
  const DetectorModel da(a), dc(c);
  int showered = 0;
  for (int i = 0; i < 50; ++i) {
    const TruthEvent ea = generateTruthEvent(a, i), ec = generateTruthEvent(c, i);
    CHECK(truthEventToYAML(ea) == truthEventToYAML(generateTruthEvent(b, i)));
    CHECK(ea.shower.empty());
    CHECK(ea.particles.size() == ec.particles.size());
    for (size_t k = 0; k < ea.particles.size(); ++k) CHECK(ea.particles[k].p4.e == ec.particles[k].p4.e);
    showered += !ec.shower.empty();
    // Without clustering the shower partons of each quark are summed back into one jet
    // (energy-conserving shower; no non-perturbative shift here).
    const ObservedEvent oa = da.apply(ea), oc = dc.apply(ec);
    CHECK(oa.objects.size() == oc.objects.size());
    for (const auto &o : oc.objects) CHECK_CLOSE(o.p4.e, ec.particles[o.truth_index].p4.e, 1e-6 * o.p4.e);
  }
  CHECK(showered == 50);
}

TEST(boosted_w_merges_into_one_large_radius_jet) {
  const TaskConfig t = sampleTask(boosted("[X]", kShower, kJets), 1, 0);
  const DetectorModel det(t);
  int n = 0, wTagged = 0, twoPartons = 0;
  Stats smallR;
  for (int i = 0; i < 2000; ++i) {
    const ObservedEvent ev = det.apply(generateTruthEvent(t, i));
    const ObservedObject *fj = leading(ev, "fatjet");
    if (!fj) continue;
    ++n;
    twoPartons += fj->n_partons == 2;
    wTagged += fj->m() > 65 && fj->m() < 100;
    const ObservedObject *j = leading(ev, "jet");
    if (j) smallR.add(j->m());
  }
  CHECK(n > 1900);
  CHECK(double(twoPartons) / n > 0.9);
  CHECK(double(wTagged) / n > 0.6);
  // W partons at pT ~ 1.2 TeV are ~0.13 apart, so R = 0.4 jets mostly contain both too.
  CHECK(smallR.mean() > 50);
}

TEST(qcd_jet_mass_grows_with_radius) {
  const TaskConfig t = sampleTask(boosted("[g, g]", kShower, kJets), 1, 0);
  const TaskConfig flat = sampleTask(boosted("[g, g]", "", kJets), 1, 0);
  const DetectorModel det(t), detFlat(flat);
  Stats small, large;
  for (int i = 0; i < 2000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const ObservedEvent ob = det.apply(ev);
    const ObservedObject *j = leading(ob, "jet"), *fj = leading(ob, "fatjet");
    if (j && fj) {
      small.add(j->m() / j->pt());
      large.add(fj->m() / fj->pt());
    }
    // Without a shower a gluon jet is massless.
    const ObservedObject *jf = leading(detFlat.apply(generateTruthEvent(flat, i)), "jet");
    if (jf) CHECK(jf->m() < 1e-6 * jf->pt());
  }
  CHECK(small.n > 1800);
  CHECK(small.mean() > 0.03 && small.mean() < 0.2);
  CHECK(large.mean() > 1.3 * small.mean());
}

TEST(extra_collection_does_not_change_jets_or_met) {
  const std::string smallOnly = "  jets:\n    - {type: jet, radius: 0.4}\n";
  const TaskConfig a = sampleTask(boosted("[X]", kShower, smallOnly), 1, 0), b = sampleTask(boosted("[X]", kShower, kJets), 1, 0);
  const DetectorModel da(a), db(b);
  for (int i = 0; i < 100; ++i) {
    const ObservedEvent oa = da.apply(generateTruthEvent(a, i)), ob = db.apply(generateTruthEvent(b, i));
    CHECK(oa.met_x == ob.met_x && oa.met_y == ob.met_y);
    std::vector<double> ja, jb;
    for (const auto &o : oa.objects) ja.push_back(o.p4.pt());
    for (const auto &o : ob.objects)
      if (o.type == "jet") jb.push_back(o.p4.pt());
    CHECK(ja == jb);
  }
}

TEST(jet_mass_response) {
  // Without a shower, a merged W -> q q jet has exactly the W mass; mass_scale and mass_resolution
  // then act on the jet mass only (the detector is otherwise perfect).
  const PriorConfig smeared = kgtest::prior(R"(
particles:
  X: {pdg: 9000023, mass: 2500, width: 1, decays: [{daughters: [W, W], br: 1}]}
  W: {pdg: 24, mass: 80.4, width: 2.09, decays: [{daughters: [q, q], br: 1}]}
  q: {pdg: 1, object: fatjet}
production: [X]
initial_state: {sqrt_s: {model: fixed, value: 2500}}
detector:
  objects:
    fatjet: {mass_scale: 1.1, mass_resolution: 0.1}
  jets: [{type: fatjet, radius: 0.8, inputs: [fatjet]}]
)");
  const TaskConfig t = sampleTask(smeared, 1, 0);
  const DetectorModel d(t);
  Stats r, pt;
  for (int i = 0; i < 4000; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    for (const auto &o : d.apply(ev).objects) {
      if (o.n_partons != 2) continue;
      const TruthParticle &w = ev.particles[ev.particles[o.truth_index].parent];
      r.add(o.m() / w.p4.m());
      pt.add(o.pt() / w.p4.pt());
    }
  }
  CHECK(r.n > 7000);
  CHECK_CLOSE(r.mean(), 1.1, 4 * r.err());
  CHECK_CLOSE(r.sd(), 0.11, 0.006);
  CHECK_CLOSE(pt.mean(), 1.0, 1e-9);
}

TEST(jet_configuration_round_trip_and_validation) {
  const PriorConfig p = kg::loadPrior(kgtest::examplePath("jet_substructure.yaml"));
  CHECK(p.shower.enabled);
  CHECK(p.detector.jets.size() == 2);
  const TaskConfig t = sampleTask(p, 5, 1);
  const TaskConfig r = taskFromYAML(taskToYAML(t));
  CHECK(r.shower.enabled);
  CHECK(taskToYAML(r) == taskToYAML(t));
  const DetectorModel dt(t), dr(r);
  for (int i = 0; i < 20; ++i) {
    const TruthEvent a = generateTruthEvent(t, i), b = generateTruthEvent(r, i);
    CHECK(!a.shower.empty());
    CHECK(truthEventToYAML(a) == truthEventToYAML(b));
    const TruthEvent back = truthEventFromYAML(truthEventToYAML(a));
    CHECK(back.shower.size() == a.shower.size());
    CHECK(observedEventToYAML(dt.apply(a)) == observedEventToYAML(dt.apply(back)));
    CHECK(observedEventToYAML(dt.apply(a)) == observedEventToYAML(dr.apply(b)));
    const ObservedEvent ob = dt.apply(a);
    CHECK(observedEventToYAML(observedEventFromYAML(observedEventToYAML(ob))) == observedEventToYAML(ob));
  }
  const std::string base = kgtest::kTwoLeptons;
  CHECK_THROWS_AS(kgtest::prior(base + "parton_shower: {alpha: 0.1}\n"), std::invalid_argument, "unknown key");
  CHECK_THROWS_AS(kgtest::prior(base + "detector: {jets: [{type: jet}, {type: jet2, inputs: [jet]}]}\n"), std::invalid_argument,
                  "already consumed");
  CHECK_THROWS_AS(kgtest::prior(base + "detector: {jets: [{radius: 0.4}]}\n"), std::invalid_argument, "needs a 'type'");
  CHECK_THROWS_AS(sampleTask(kgtest::prior(base + "parton_shower: {kt_min: 0}\n"), 1, 0), std::exception, "kt_min");
}
