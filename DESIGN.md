# kinegen design

kinegen is a parton-level kinematic generator for synthetic pretraining and
meta-learning. It keeps the spirit of the original `myGenerator.h`/`test.cc`
(fast, two-body/three-body phase space, no matrix elements) and restructures it
around four stages:

```
PriorConfig ──sampleTask──▶ TaskConfig ──generateTruthEvent──▶ TruthEvent ──DetectorModel::apply──▶ ObservedEvent
     │                          │ (task.environment)                                   ▲
     └──sampleEnvironment──▶ Environment ──────────────────────────────────────────────┘
```

| Stage | Type | Randomness resolved | Header |
|---|---|---|---|
| Prior | `PriorConfig` | none: distributions over everything | `kg/Prior.h` |
| Task | `TaskConfig` | masses, widths, BRs, angular coefficients, processes/topology, multiplicity, mixture fractions, feature schema, default environment | `kg/Task.h` |
| Environment | `Environment` | detector parameters (environment scope) | `kg/Detector.h` |
| Truth event | `TruthEvent` | event-scoped physics parameters, process component, decay modes, off-shell masses, √ŝ, rapidity, recoil, angles | `kg/Generator.h`, `kg/Event.h` |
| Observed event | `ObservedEvent` | event-scoped detector parameters, efficiency, misID, smearing, fakes, MET, features | `kg/Generator.h`, `kg/Event.h` |

The core library (`libkinegen`) depends only on yaml-cpp. ROOT output is a
separate optional library (`libkinegen_root`), and Python bindings are an
optional pybind11 module.

## C++ API

```cpp
#include "kg/kinegen.h"

kg::PriorConfig prior = kg::loadPrior("examples/randomized_tasks.yaml");
kg::TaskConfig task = kg::sampleTask(prior, /*master_seed=*/1, /*task_index=*/0);
kg::DetectorModel detector(task);
kg::Environment other = kg::sampleEnvironment(prior.detector, /*master_seed=*/7, /*env_index=*/3);

for (uint64_t i = 0; i < 1000; ++i) {
  kg::TruthEvent truth = kg::generateTruthEvent(task, i);          // or generateTruthEvent(task): uses task.next_event
  kg::ObservedEvent reco = detector.apply(truth, task.environment); // or detector.apply(truth, other)
  for (const kg::ObservedObject &o : reco.objects) { o.p4; o.type; o.features; o.truth_index; }
}
```

* `sampleTask(prior)` uses `prior.seed` and index 0; `kg::TaskSampler` iterates indices.
* `kg::generateEvent(task, env, i)` returns truth and observed event together.
* `taskToYAML`/`taskFromYAML`/`saveTask`/`loadTask`, `environmentToYAML`/`environmentFromYAML`,
  `truthEventToYAML`/`observedEventToYAML` (and their inverses) serialise everything exactly
  (17 significant digits).
* `kg::assignSplit(task, SplitSpec)` assigns `train`/`val`/`test` (or any labels) by task, topology,
  process set, environment or parameter range.
* `kg::RootWriter` writes `events`, `tasks` and `environments` trees, or the legacy `test`/`info` trees.

## Randomness and reproducibility

All random numbers come from `kg::Rng` (xoshiro256** seeded through SplitMix64).
`std::*_distribution` is not used, so results do not depend on the standard library.
Seeds are derived, never shared:

| Stream | Seed |
|---|---|
| task *k* of master seed *S* | `task_seed = mix(S, "task", k)` |
| default environment of a task | `mix(task_seed, "env")` |
| environment *k* of master seed *S* | `mix(S, "env", k)` |
| truth event *i* | `mix(task_seed, "evn", i)` |
| detector response | `mix(env_seed, event_seed, "det")` |
| features | `mix(detector_seed, "fea")` |

Consequences:

* Event *i* of a task is a pure function of `(task, i)`: events can be generated in any order or in
  parallel (`generate_batch(..., threads=N)` in Python) and regenerated on demand.
* The truth record does not depend on the detector or feature configuration, and observing the same
  truth event in different environments uses independent detector randomness.
* Rejected task draws (see below) consume the task stream deterministically, so `(S, k)` still
  identifies one task. `TaskConfig::attempts` records how many draws were needed.
* A serialised task (`taskToYAML`) contains the frozen parameter set, the resolved particle table,
  components, initial-state model, feature schema, environment and seeds, and reproduces every event
  bit-for-bit after `taskFromYAML` (tested).

## Parameters

Every number in a prior is a `ParamDef` in a `ParamSet` (`kg/Param.h`):

```yaml
mass: 172.5                                       # fixed
mass: {dist: uniform, min: 100, max: 200}         # sampled per task (default scope)
mass: {dist: uniform, min: 0.1, max: 0.9, ratio_of: X}   # value = draw * X.mass
mass: "0.5 * X.mass + 10"                         # expression of other parameters
mass: {dist: gaussian, mean: 91, sigma: 1, scope: event} # re-drawn every event
```

Distributions: `fixed`, `uniform`, `loguniform`, `gaussian`, `truncated_gaussian`, `breit_wigner`
(optionally `relativistic`), `categorical` (`values`, `weights`), `exponential` (`scale`, `offset`),
`poisson` (`mean`).
Expressions support `+ - * / ^`, parentheses, `min max sqrt abs log exp pow`, and reference any
parameter by name. A bare particle name in `ratio_of` means its mass.

Inline parameters are named by convention, which is also how tasks and environments bind them
after deserialisation: `<particle>.mass`, `<particle>.width`, `<particle>.decay<i>.br`,
`<particle>.decay<i>.angular.c<k>`, `process.<name>.angular.c<k>`, `initial_state.<field>`,
`detector.<type>.<field>`, `detector.fakes.<field>`, `detector.met.<field>`,
`feature.<name>.<field>`. Named free parameters go in a top-level `parameters:` map.

`ParamSet::compile()` topologically sorts dependencies (cycles are errors) and propagates event
scope: anything depending on an event-scoped parameter is event-scoped. Branching fractions,
angular coefficients, process weights and feature parameters must be task-level.

## Physics definition

```yaml
particles:
  t:
    pdg: 6
    charge: 0.667
    mass: 172.5
    width: 1.42
    decays:                       # 2- or 3-body; recursive
      - {daughters: [W+, b], br: 1.0, angular: isotropic}
  W+:
    mass: 80.4
    width: 2.09
    decay_choice: event           # or: task (one mode fixed per task)
    decays:
      - {daughters: [mu+, vm], br: 0.108, angular: {type: polynomial, coeffs: [1, 2, 1], axis: helicity}}
      - {daughters: [q, q], br: 0.892}
  b: {pdg: 5, mass: 4.7}           # stable (no decays); object type from pdg -> bjet
  vm: {pdg: 14}                    # neutrinos -> invisible
  N: {mass: 100, visible: false}   # any invisible BSM state
```

* Names are arbitrary strings. Object types default from the PDG id (`jet`, `bjet`, `electron`,
  `muon`, `tau`, `photon`, `invisible`); unknown visible PDG ids require an explicit `object:`.
* Decays form an explicit graph: each mode lists daughters, a branching fraction, an angular model
  and, for three-body decays, `model: phase_space` (uniform Dalitz density, default) or
  `legacy_sequential` (myGenerator's `0 → 1 Y, Y → 2 3` with uniform `m_Y`).
* Angular models (`kg::AngularModel`) give the density of cos θ* of the first daughter in the parent
  rest frame about the `helicity` axis (parent flight direction, default for decays) or the `beam`
  axis (default for production): `isotropic` or `polynomial` with `coeffs` c₀…cₙ
  (f = Σ cᵢ cosⁱθ, validated non-negative). `W` reproduces myGenerator's `1 − cos θ` mode. For
  three-body decays the configuration is rolled uniformly about the first daughter's direction.
* Masses with `width > 0` are drawn per event from a relativistic Breit-Wigner, truncated exactly
  (inverse CDF) to `[threshold, M + bw_cutoff·Γ]` and, below, to `M − bw_cutoff·Γ`
  (`bw_cutoff` defaults to 15, MadGraph's convention; `0` keeps only physical limits). Daughter
  masses are drawn jointly so that their sum never exceeds the parent mass.

Processes:

```yaml
processes:
  - {name: pair, weight: 2, production: [X, X]}                      # 1-3 produced particles
  - {name: assoc, production: [X, {choice: [Zp, Z], weights: [1, 3]}]}
  - {name: multi, production: [X], multiplicity: {values: [1, 2, 3], weights: [5, 3, 2]}}
  - {name: res, production: [Zp], angular: {type: polynomial, coeffs: [1, 0, 1]}}
  - name: W                                   # per-event subprocesses: hard extra partons
    subprocesses:
      - {production: [W], weight: 0.6}
      - {production: [W, q], weight: 0.3, angular: {type: polynomial, coeffs: [1, 0, 2]}}
      - {production: [W, q, q], weight: 0.1}
production: [A, B]     # shorthand for a single process
```

`subprocesses` lists production alternatives of one process, chosen event by event with the given
(task-level parameter) weights; `TruthEvent::subprocess` records the choice. They model matched
multi-jet samples (V + 0/1/2 jets, tt̄ + 0/1 jet) with the extra partons produced in the hard
process (two- or three-body phase space above `threshold_offset`), unlike the soft independent recoil
jets. The process (and hence the pretraining family label) is the same for all alternatives.

Task assembly:

```yaml
task:
  components: {dist: categorical, values: [1, 2]}   # processes per task (event label = component)
  fractions: {dirichlet: 1.0}                        # or equal, or [0.7, 0.3] (a list covering every
                                                     # process applies in catalogue order)
  distinct_processes: true
  validation: reject                                 # or error
  max_attempts: 1000
```

Initial state (`initial_state.*` parameters):

```yaml
initial_state:
  sqrt_s: {model: exponential_pdf, ebeam: 6500, tau: 0.05, tau_slope: 0.02, threshold_offset: 0}
  #        pdf: effective parton densities f(x) ~ x^-a (1-x)^b, parameters a, b (see below)
  #        threshold_exponential: sqrt_s = threshold + Exp(scale)
  #        fixed: value ;  uniform: min, max
  rapidity: from_pdf              # or any parameter, e.g. {dist: gaussian, mean: 0, sigma: 1, scope: event}
  recoil: {model: exponential, min_pt: 10, mean_pt: 15, n_jets: 1, emit: false, object: jet, rapidity: <param>}
  #        model: param with pt: <param>;  or recoil: none (default)
  #        n_jets: e.g. {dist: poisson, mean: 1.5, scope: event} for multi-jet ISR
```

`pdf` is the more realistic luminosity model: both beams carry partons with density
f(x) = x^(−a) (1 − x)^b, so dL/dln τ dy ∝ τ^(1−a) (1 − x₁)^b (1 − x₂)^b with τ = x₁x₂ = ŝ/s and
y = ½ ln(x₁/x₂). ln τ is drawn exactly above the threshold from the tabulated marginal (inverse CDF,
efficient for any mass) and y from its conditional density, which also gives single resonances a
mass-dependent rapidity distribution. `exponential_pdf` is the legacy luminosity model: x₁, x₂ ~ Exp(τ) on (0, 1] with
τ = tau + tau_slope·threshold/TeV, accepted when √ŝ ≥ Σm + threshold_offset. For a single resonance
the rapidity is drawn from the same density conditional on x₁x₂ = m²/s.

The recoil (ISR) consists of `n_jets` massless jets (default 1; 0 is allowed and leaves the hard
system with no transverse momentum). Each jet draws its own pT (`exponential`: min_pt +
Exp(mean_pt − min_pt); `param`: a fresh draw of `pt` per jet when it is an event-scoped
distribution), its own rapidity (`rapidity`, default uniform in ±2.5 per jet) and a uniform
azimuth. The hard system takes the opposite of their summed transverse momentum. The jets are
recorded as truth particles named `recoil` (status `Recoil`, treated as unclustered energy by the
detector) or, with `emit: true`, as visible final-state objects of type `object`. Incoming partons
are reconstructed so that `Σ incoming = Σ final + recoil` exactly. With `n_jets: 1` the result is
identical to the single-recoil model. The ISR jet kinematics are uncorrelated with each other and
with the hard process; this is a recoil model, not a parton shower.

## Truth event

`TruthEvent::particles` is a flat list. Indices 0 and 1 are the incoming partons, followed by the recoil
jets (if any), then the decay tree in breadth-first order. Each `TruthParticle` has
name, PDG id, charge, four-vector, generated mass, `parent` (−1 for produced particles), `children`,
status (`Incoming`, `Recoil`, `Decayed`, `Final`), species index into the task's particle table,
`branch` (which produced particle it descends from), `depth`, chosen `decay_mode` and object type.
Helpers: `finalState()`, `producedParticles()`, `ancestors(i)`, `isDescendant(i, a)`, `trueMET()`,
`topology()` (e.g. `t>(W+>(mu+,vm),b) + tbar>(W->(q,q),bbar)`). Event-scoped parameter draws are
kept in `event_params`.

Generation of one event: draw event parameters → component → decay modes top-down → mass
thresholds bottom-up → masses top-down → √ŝ, rapidity, recoil → production (1, 2 or 3 bodies) →
decays in each parent's rest frame, boosted to the lab.

## Observation model

Per object type (keys under `detector.objects`, plus `default` for unlisted types; missing entries
mean a perfect detector):

| Field | Meaning |
|---|---|
| `efficiency: {plateau, pt50, width}` or a number | ε(pT) = plateau · ½[1 + erf((pT − pt50)/(√2·width))] |
| `eta_max`, `pt_min` | acceptance; `pt_min` of the *reconstructed* type is applied after smearing |
| `resolution: {stochastic, constant, noise, pt_term}` | σ/E = √(a²/E + b² + c²/E² + (d·pT)²) |
| `scale` | response scale (four-vector multiplied by scale·(1 + σ·N(0,1))) |
| `sigma_eta`, `sigma_phi` | Gaussian angular smearing |
| `misid: {type: p}` | probability to be reconstructed as another type (Σ ≤ 1) |
| `response_map: {eta_edges, pt_edges, scale, resolution, efficiency}` | binned response (below) |
| `resolution_scale` | multiplies the map resolution (default 1; a parameter, so it can be randomised) |

A binned response map replaces the parametric formulas with tables in (|η|, pT) of the true object,
for example measured from a full detector simulation:

```yaml
jet:
  efficiency: {dist: uniform, min: 0.97, max: 1.0}     # plateau factor applied on top of the map
  resolution_scale: {dist: uniform, min: 0.85, max: 1.15}
  response_map:
    eta_edges: [0, 1.3, 2.5, 4.0]                      # |eta| bin edges
    pt_edges: [15, 30, 60, 120, 1.0e5]                 # true-pT bin edges
    scale: [[0.97, 0.98, 0.99, 0.99], ...]             # median reco/true pT, [eta][pT]
    resolution: [[0.20, 0.14, 0.10, 0.08], ...]        # relative resolution sigma/pT
    efficiency: [[0.80, 0.97, 0.99, 1.00], ...]        # optional; else the parametric turn-on is used
```

With a map, an object in bin (i, j) is kept with probability `efficiency` × map efficiency[i][j]
(the parametric `pt50`/`width` turn-on when the map has no efficiency table) and its four-vector is
multiplied by `scale` · map scale[i][j] · (1 + `resolution_scale` · map resolution[i][j] · g), where g
includes the correlated-group fluctuation as for the parametric response. Values outside the edges
use the nearest bin. Map entries are fixed numbers (randomise through `resolution_scale`, `scale` and
the efficiency plateau); tasks serialise the full map.

Composite objects (`reconstruct_as` on a decay mode) model objects that a detector reconstructs from
several particles, most importantly hadronic taus:

```yaml
tau-:
  mass: 1.777
  decays:
    - {daughters: [e-, ve, vt], br: 0.178}                      # leptonic: the electron is its own object
    - {daughters: [pi-, vt], br: 0.108, reconstruct_as: tau}    # hadronic: one "tau" object
    - {daughters: [rho-, vt], br: 0.255, reconstruct_as: tau}   # rho- -> pi- pi0, pi0 -> gamma gamma
    - {daughters: [a1-, vt], br: 0.090, reconstruct_as: tau}    # three prongs
```

When a particle decays through such a mode, the four-momenta (and charges) of all its *visible*
final-state descendants are summed into one reconstruction candidate of the given type, which then
receives that type's efficiency, smearing, misidentification (e.g. `tau: {misid: {jet: 0.3}}`) and
correlation-group shift like any other object; invisible descendants (neutrinos) stay in the MET.
The same mechanism gives merged large-radius jets (`W -> q q` with `reconstruct_as: fatjet`). With
nested composites the outermost one wins. Observed objects carry `n_constituents` and `n_charged`
(prong count), available as `observed_proxy` feature quantities and as `n_constituents`/`n_charged`
arrays in `generate_batch`; truth particles record `reconstruct_as`. Without composite modes the
detector output is unchanged.

Correlated response (`detector.correlations`), e.g. a jet energy scale shared by light and b jets:

```yaml
correlations:
  jes: {types: [jet, bjet], scale_sigma: 0.02, resolution_correlation: 0.3, unclustered: true}
```

Per event and group a shift δ ~ N(0, scale_sigma) multiplies the response of every object of the
group's types, and the resolution fluctuation of each object is g = √ρ·g_common + √(1 − ρ)·g_object
with ρ = `resolution_correlation`, so the marginal resolution is unchanged while objects in the
same event fluctuate together. With `unclustered: true` the shift also applies to the unclustered
recoil in the MET (at most one group may do so). A type can belong to one group only. Groups are
matched on the true object type, as for all other response parameters. Each event records
`detector.correlations.<group>.shift` in `ObservedEvent::detector_params`; both parameters can be
randomised per environment or per event like any other. Because the MET is built from the smeared
objects, correlated shifts propagate into it consistently.

`fakes: {rate, types, pt_min, pt_scale, eta_max}` adds Poisson(rate) fake objects with
pT = pt_min + Exp(pt_scale), uniform η and φ. `met: {soft_sigma, sumet_coeff, recoil_scale,
recoil_resolution, lost_to_met}` reconstructs

  MET = −[Σ observed pT + (1 − lost_to_met)·Σ lost-visible pT + recoil_scale·(1 + recoil_resolution·N)·pT(recoil)] + N(0, √(soft_sigma² + sumet_coeff²·ΣET))

which equals the true invisible pT for a perfect detector (tested). Objects are sorted by pT unless
`sort_by_pt: false`. Detector values outside their physical range (e.g. a sampled efficiency above
one) are errors, not clipped; use `truncated_gaussian` priors.

`ObservedEvent` holds `objects` (four-vector, type, `type_code` into `TaskConfig::object_types`,
charge, `features`, `truth_index` or −1 for fakes, `true_type`, `is_fake`), `met_x/met_y`, `sum_et`,
counts of lost/fake/misidentified objects, ids, seeds, the component label, the feature schema and
any event-scoped detector draws.

## Auxiliary features

```yaml
features:
  - {name: noise, kind: gaussian, mean: 0, sigma: 1}
  - {name: u, kind: uniform, min: 0, max: 1, presence_prob: 0.5}          # task-dependent schema
  - {name: pt_proxy, kind: truth_proxy, quantity: pt, rel_sigma: 0.1, fake_policy: observed}
  - {name: eta_copy, kind: observed_proxy, quantity: eta, abs_sigma: 0.01}
  - {name: btag, kind: label_correlated, label: truth_type, shift: {dist: gaussian, mean: 0, sigma: 2}, applies_to: [jet, bjet]}
  - {name: spurious, kind: label_correlated, label: component, shift: {dist: uniform, min: -3, max: 3}}
  - {name: res_tag, kind: environment, param: detector.jet.resolution.stochastic, coeff: 1, sigma: 0.1}
```

Truth quantities: `pt eta phi m e px py pz charge pdg is_fake depth branch parent_pdg parent_mass`;
observed quantities additionally include `n_constituents` and `n_charged`.
Labels: `component`, `process`, `truth_type` (truth particle name or `fake`), `parent`, `is_fake`.
For label-correlated features a shift per category is drawn once per task, so both the presence
and the strength of a spurious correlation vary between tasks. `missing_prob` blanks values per
object (NaN), `applies_to` restricts object types (NaN elsewhere) and `presence_prob` drops a
feature from a task's schema. Truth ancestry is never exposed as a feature unless configured.

## Validation

Static (at `parsePrior`): unknown keys, unknown particles/daughters, more than three daughters or
produced particles, cyclic decay graphs (the cycle is printed), stable/decays contradictions,
uninferable object types, fixed BRs outside [0, 1] or not summing to one (unless
`normalize_br: true`), event-scoped BRs/angles/feature parameters, impossible component counts,
unknown feature quantities/labels/detector parameters, dependency cycles in parameters.

Per task (`validateTask`, `InvalidTask`): sampled BR sums, negative masses/widths, negative angular
polynomials, kinematically forbidden decays (narrow parents below the sum of daughter minima,
recursively, worst case over modes), thresholds beyond the Breit-Wigner cutoff, production
thresholds above the collider energy, invalid sampled detector values. With
`task.validation: reject` invalid draws are redrawn (conditioning the prior on validity); otherwise,
or after `max_attempts`, an error is raised. Nothing is silently clipped or rescaled.

Per event: event-scoped parameters that violate thresholds raise `std::runtime_error`.

## Output

ROOT (`kg::RootWriter`, `kg-generate --out`): `events` (ids, label, `obj_*` vectors, `obj_feat`
flattened with `n_feat` per object, MET, truth arrays `truth_*` unless `--no-truth`), `tasks`
(ids, seeds, topology, processes and fractions, feature names, object types, all parameter values,
full task YAML) and `environments` (id, seed, parameters, YAML). The legacy format writes the
original `test` and `info` trees.

Python (`kinegen.generate_batch`) returns padded NumPy arrays: `p4`, `ptetaphim` (n, M, 4),
`type` (−1 padding), `mask`, `features` (NaN = missing or padding), `truth_index`, `is_fake`, `met`,
`label`, `n_objects`, `event_index`, plus `feature_names` and `object_types`. Generation releases
the GIL and can use several threads.

## Splits

`SplitSpec{by, fractions, salt, parameter, edges, labels}`: `Task`, `Topology`, `Process` and
`Environment` modes hash a key (so e.g. a topology is never split across train and test);
`Parameter` bins a task-level parameter (e.g. `X.mass` below/above 1.5 TeV).

## Compatibility

`myGenerator.h`, `test.cc` and `process.card` are unchanged. `kg::convertLegacyCard` /
`kg-convert-card` translate a card into an equivalent YAML prior (legacy luminosity with the 1 GeV
offset, 10 + Exp(5) GeV recoil, isotropic decays, the "too-heavy daughters → three-body with the
heavier daughter replaced by its decay products" rule as `legacy_sequential`, per-event uniform
mass ranges). `kg-generate --legacy-card process.card` writes the `test`/`info` trees. Means and
RMS of all legacy branches agree with `test.cc` within statistical precision (checked on 2×10⁵
events) and a unit test compares the kinematics against `myGenerator.h` directly. Bitwise equality
with `test.cc` is not possible because ROOT's `TRandom3` is no longer used.

## Limitations

* No matrix elements, spin correlations beyond per-decay angular polynomials, parton showers,
  hadronisation or interference; widths enter only through Breit-Wigner line shapes.
* ISR jets are independent massless objects (no shower, no colour coherence); there is no
  underlying event or pile-up beyond the MET soft term and fakes.
* Detector response is per object apart from correlation groups and composites, without overlap
  removal or isolation.
* Parton momentum fractions are not re-checked after adding ISR jets, so very hard ISR can give
  x₁ or x₂ slightly above one.
