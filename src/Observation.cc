#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "kg/Generator.h"

namespace kg {

int ObservedEvent::featureIndex(const std::string &name) const {
  if (!feature_names) return -1;
  for (size_t i = 0; i < feature_names->size(); ++i)
    if ((*feature_names)[i] == name) return static_cast<int>(i);
  return -1;
}

double ObservedEvent::feature(size_t object, const std::string &name) const {
  const int i = featureIndex(name);
  if (i < 0) return std::numeric_limits<double>::quiet_NaN();
  return objects.at(object).features.at(i);
}

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

double truthQuantity(const TruthEvent &ev, const TaskConfig &t, const TruthParticle &p, const std::string &q) {
  if (q == "pt") return p.p4.pt();
  if (q == "eta") return p.p4.eta();
  if (q == "phi") return p.p4.phi();
  if (q == "m") return p.mass;
  if (q == "e") return p.p4.e;
  if (q == "px") return p.p4.px;
  if (q == "py") return p.p4.py;
  if (q == "pz") return p.p4.pz;
  if (q == "charge") return p.charge;
  if (q == "pdg") return p.pdg;
  if (q == "is_fake") return 0;
  if (q == "depth") return p.depth;
  if (q == "branch") return p.branch;
  if (q == "parent_pdg") return p.parent >= 0 ? ev.particles[p.parent].pdg : 0;
  if (q == "parent_mass") return p.parent >= 0 ? ev.particles[p.parent].mass : kNaN;
  return kNaN;
}

double observedQuantity(const ObservedObject &o, const std::string &q) {
  if (q == "pt") return o.p4.pt();
  if (q == "eta") return o.p4.eta();
  if (q == "phi") return o.p4.phi();
  if (q == "m") return o.p4.m();
  if (q == "e") return o.p4.e;
  if (q == "px") return o.p4.px;
  if (q == "py") return o.p4.py;
  if (q == "pz") return o.p4.pz;
  if (q == "charge") return o.charge;
  if (q == "n_constituents") return o.n_constituents;
  if (q == "n_charged") return o.n_charged;
  return kNaN;
}

double noisy(double x, const FeatureSpec &f, Rng &rng) {
  return x + f.get("bias") + f.get("rel_sigma") * std::fabs(x) * rng.gauss() + f.get("abs_sigma") * rng.gauss();
}

void computeFeatures(const TaskConfig &t, const TruthEvent &ev, const Environment &env, const std::vector<double> &envValues,
                     ObservedObject &o, Rng &rng) {
  o.features.assign(t.features.size(), kNaN);
  const TruthParticle *tp = o.truth_index >= 0 ? &ev.particles[o.truth_index] : nullptr;
  for (size_t i = 0; i < t.features.size(); ++i) {
    const FeatureSpec &f = t.features[i];
    // Always consume the same random numbers so that features stay aligned across objects.
    const double g = rng.gauss();
    const double u = rng.uniform();
    const double miss = rng.uniform();
    if (!f.applies(o.type) || miss < f.get("missing_prob")) continue;
    double x = kNaN;
    switch (f.kind) {
      case FeatureKind::Gaussian: x = f.get("mean") + f.get("sigma") * g; break;
      case FeatureKind::Uniform: x = f.get("min") + (f.get("max") - f.get("min")) * u; break;
      case FeatureKind::TruthProxy:
        if (tp) x = noisy(truthQuantity(ev, t, *tp, f.quantity), f, rng);
        else if (f.quantity == "is_fake") x = noisy(1.0, f, rng);
        else if (f.fake_policy == "observed") x = noisy(observedQuantity(o, f.quantity), f, rng);
        break;
      case FeatureKind::ObservedProxy: x = noisy(observedQuantity(o, f.quantity), f, rng); break;
      case FeatureKind::LabelCorrelated: {
        std::string cat;
        if (f.label == "component") cat = std::to_string(ev.component);
        else if (f.label == "process") cat = ev.process;
        else if (f.label == "is_fake") cat = o.is_fake ? "1" : "0";
        else if (f.label == "truth_type") cat = tp ? tp->name : "fake";
        else if (f.label == "parent") cat = !tp ? "fake" : (tp->parent >= 0 ? ev.particles[tp->parent].name : "none");
        auto it = f.shifts.find(cat);
        x = f.get("base_mean") + (it == f.shifts.end() ? 0.0 : it->second) + f.get("sigma") * g;
        break;
      }
      case FeatureKind::Environment: {
        const int idx = env.detector->params.find(f.env_param);
        if (idx >= 0) x = f.get("offset") + f.get("coeff") * envValues[idx] + f.get("sigma") * g;
        break;
      }
    }
    o.features[i] = x;
  }
}

} // namespace

ObservedEvent DetectorModel::apply(const TruthEvent &truth, const Environment &env) const {
  const TaskConfig &t = *task_;
  ObservedEvent out;
  out.task_id = truth.task_id;
  out.env_id = env.env_id;
  out.event_index = truth.event_index;
  out.truth_seed = truth.seed;
  out.detector_seed = mixSeed(env.seed, truth.seed, kStreamDetector);
  out.component = truth.component;
  out.process = truth.process;
  out.feature_names = t.feature_names;
  Rng rng(out.detector_seed);
  Rng frng(mixSeed(out.detector_seed, kStreamFeature));

  // Event-scoped detector parameters.
  const DetectorState *state = &env.state;
  DetectorState eventState;
  std::vector<double> envValues = env.values;
  if (env.hasEventParams()) {
    env.detector->params.resampleEvent(envValues, rng);
    eventState = buildDetectorState(*env.detector, envValues, "environment '" + env.env_id + "' event " + std::to_string(truth.event_index));
    state = &eventState;
    for (size_t i = 0; i < envValues.size(); ++i)
      if (env.detector->params.isEvent(static_cast<int>(i)))
        out.detector_params.emplace_back(env.detector->params.def(static_cast<int>(i)).name, envValues[i]);
  }
  const DetectorState &ds = *state;

  // Correlated response: one scale shift and one common resolution draw per group.
  std::vector<double> groupShift(ds.groups.size()), groupCommon(ds.groups.size());
  for (size_t g = 0; g < ds.groups.size(); ++g) {
    groupShift[g] = ds.groups[g].scale_sigma * rng.gauss();
    groupCommon[g] = rng.gauss();
    out.detector_params.emplace_back("detector.correlations." + ds.groups[g].name + ".shift", groupShift[g]);
  }

  double sumx = 0, sumy = 0;         // observed visible objects
  double lostx = 0, losty = 0;       // truth pT of visible objects that were not reconstructed
  double recx = 0, recy = 0;         // unclustered recoil
  double sumEt = 0;

  // Reconstruction candidates: visible final-state particles, except that the
  // visible descendants of a particle decayed through a `reconstruct_as` mode are
  // merged into one composite object (e.g. a hadronic tau). The outermost merging
  // ancestor wins. Without composite modes this is exactly the list of visible
  // final-state particles in index order.
  struct Candidate {
    int index;
    P4 p4;
    double charge;
    std::string type;
    int n_constituents, n_charged;
  };
  const size_t np = truth.particles.size();
  std::vector<int> owner(np, -1);
  for (size_t i = 0; i < np; ++i) {
    const int par = truth.particles[i].parent;
    if (par >= 0) owner[i] = owner[par] >= 0 ? owner[par] : (truth.particles[par].reconstruct_as.empty() ? -1 : par);
  }
  std::vector<Candidate> cands;
  std::vector<int> compositeSlot(np, -1);
  for (size_t i = 0; i < np; ++i) {
    const TruthParticle &p = truth.particles[i];
    if (p.status == TruthStatus::Recoil) {
      recx += p.p4.px;
      recy += p.p4.py;
      continue;
    }
    if (owner[i] < 0 && !p.reconstruct_as.empty()) {
      compositeSlot[i] = static_cast<int>(cands.size());
      cands.push_back({p.index, P4(), 0.0, p.reconstruct_as, 0, 0});
      continue;
    }
    if (!p.visible()) continue;
    if (owner[i] >= 0) {
      Candidate &c = cands[compositeSlot[owner[i]]];
      c.p4 += p.p4;
      c.charge += p.charge;
      ++c.n_constituents;
      c.n_charged += p.charge != 0;
      continue;
    }
    cands.push_back({p.index, p.p4, p.charge, p.object_type, 1, p.charge != 0 ? 1 : 0});
  }

  for (const Candidate &c : cands) {
    if (c.n_constituents == 0) continue; // composite with only invisible descendants
    const Response &r = ds.response(c.type);
    // Draw all random numbers up front so that later choices do not shift streams.
    const double uEff = rng.uniform(), uMis = rng.uniform();
    const double gE = rng.gauss(), gEta = rng.gauss(), gPhi = rng.gauss();
    const double pt = c.p4.pt(), eta = c.p4.eta();
    double eff, mean, sigma;
    r.lookup(c.p4.e, pt, eta, eff, mean, sigma);
    bool keep = std::fabs(eta) <= r.eta_max && uEff < eff;
    std::string type = c.type;
    if (keep && !r.misid.empty()) {
      double acc = 0;
      for (const auto &m : r.misid) {
        acc += m.second;
        if (uMis < acc) {
          type = m.first;
          break;
        }
      }
    }
    P4 meas;
    if (keep) {
      double shift = 0, g = gE;
      auto grp = ds.group_of.find(c.type);
      if (grp != ds.group_of.end()) {
        const double rho = ds.groups[grp->second].rho;
        shift = groupShift[grp->second];
        g = std::sqrt(rho) * groupCommon[grp->second] + std::sqrt(1 - rho) * gE;
      }
      const double f = r.scale * mean * (1.0 + shift) * (1.0 + sigma * g);
      if (f <= 0) keep = false;
      else {
        meas = P4::fromPtEtaPhiM(pt * f, eta + r.sigma_eta * gEta, wrapPhi(c.p4.phi() + r.sigma_phi * gPhi), std::max(0.0, c.p4.m()) * f);
        if (meas.pt() < ds.response(type).pt_min) keep = false;
      }
    }
    if (!keep) {
      ++out.n_lost;
      lostx += c.p4.px;
      losty += c.p4.py;
      continue;
    }
    ObservedObject o;
    o.p4 = meas;
    o.type = type;
    o.type_code = t.typeCode(type);
    o.charge = c.charge;
    o.truth_index = c.index;
    o.true_type = c.type;
    o.n_constituents = c.n_constituents;
    o.n_charged = c.n_charged;
    if (type != c.type) ++out.n_misid;
    out.objects.push_back(std::move(o));
  }

  // Fakes.
  const int nf = ds.fake_rate > 0 ? rng.poisson(ds.fake_rate) : 0;
  for (int k = 0; k < nf; ++k) {
    ObservedObject o;
    o.type = ds.fake_types[rng.categorical(ds.fake_weights)];
    const double pt = ds.fake_pt_min + rng.exponential(ds.fake_pt_scale);
    const double eta = rng.uniform(-ds.fake_eta_max, ds.fake_eta_max);
    const double phi = rng.uniform(-M_PI, M_PI);
    o.p4 = P4::fromPtEtaPhiM(pt, eta, phi, 0.0);
    o.type_code = t.typeCode(o.type);
    o.is_fake = true;
    o.true_type = "fake";
    out.objects.push_back(std::move(o));
  }
  out.n_fakes = nf;

  for (const auto &o : out.objects) {
    sumx += o.p4.px;
    sumy += o.p4.py;
    sumEt += o.p4.pt();
  }
  // Missing transverse momentum: everything not measured as objects, with the
  // unclustered recoil measured with scale/resolution and a soft term.
  double rs = ds.recoil_scale * (1.0 + ds.recoil_resolution * rng.gauss());
  for (size_t g = 0; g < ds.groups.size(); ++g)
    if (ds.groups[g].unclustered) rs *= 1.0 + groupShift[g];
  sumEt += std::hypot(recx, recy);
  const double soft = std::sqrt(ds.soft_sigma * ds.soft_sigma + ds.sumet_coeff * ds.sumet_coeff * sumEt);
  const double keepLost = 1.0 - ds.lost_to_met;
  out.met_x = -(sumx + keepLost * lostx + rs * recx) + soft * rng.gauss();
  out.met_y = -(sumy + keepLost * losty + rs * recy) + soft * rng.gauss();
  out.sum_et = sumEt;

  if (ds.sort_by_pt)
    std::stable_sort(out.objects.begin(), out.objects.end(), [](const ObservedObject &a, const ObservedObject &b) { return a.p4.pt() > b.p4.pt(); });
  for (auto &o : out.objects) computeFeatures(t, truth, env, envValues, o, frng);
  return out;
}

} // namespace kg
