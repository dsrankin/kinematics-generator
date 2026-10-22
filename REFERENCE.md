# kinegen reference

Complete reference for the prior YAML schema (every key, its default and the parameter name it
creates), conventions, the Python API, the command-line tools and the output formats.
[DESIGN.md](DESIGN.md) explains the models behind the keys; [README.md](README.md) covers
installation and a quick start.

* [Conventions](#conventions)
* [Parameters](#parameters)
* [Prior schema](#prior-schema): [top level](#top-level), [particles](#particles),
  [processes](#processes), [task](#task), [initial_state](#initial_state),
  [parton_shower](#parton_shower), [underlying_event](#underlying_event),
  [detector](#detector), [features](#features)
* [Parameter names](#parameter-names)
* [Validation and errors](#validation-and-errors)
* [Python API](#python-api)
* [Command-line tools](#command-line-tools)
* [Output formats](#output-formats)
* [Splits](#splits)
* [Using the C++ library](#using-the-c-library)

## Conventions

* Energies, momenta and masses in GeV; angles in radians; φ in (−π, π].
* `ebeam` is the energy per beam: `ebeam: 6500` means √s = 13 TeV.
* Four-vectors are (px, py, pz, E); the beam is the z axis.
* A particle's charge, mass and width come only from its YAML entry (defaults 0), never from the
  PDG id: `l: {pdg: 13}` is a neutral, massless "muon". Only the reconstructed object type is
  inferred from the PDG id (see [particles](#particles)).
* "Task level" means drawn once per task (or once per environment for detector parameters);
  "event level" means redrawn for every event.

## Parameters

Almost every number in a prior is a parameter. A parameter value can be written as:

| Form | Meaning |
|---|---|
| `172.5` | fixed number |
| `"0.5 * X.mass + 10"` | expression (a string that is not a number) |
| `{expr: "...", scope: ...}` | expression, explicit form |
| `{value: 3.0, scope: ...}` | fixed number, explicit form |
| `{dist: <name>, ..., scope: ..., ratio_of: ...}` | distribution |

Distributions (`dist`) and their keys (required unless a default is given):

| `dist` | Keys |
|---|---|
| `fixed` | `value` |
| `uniform` | `min`, `max` |
| `loguniform` (alias `log_uniform`) | `min`, `max` (both > 0) |
| `gaussian` (alias `normal`) | `mean`, `sigma` |
| `truncated_gaussian` | `mean`, `sigma`, `min` (default −∞), `max` (default +∞) |
| `breit_wigner` | `mean`, `width` (full width at half maximum), `min`, `max` (default ±∞), `relativistic` (default false) |
| `categorical` | `values` (list), `weights` (list, default equal) |
| `exponential` | `scale` (mean), `offset` (default 0): offset + Exp(scale) |
| `poisson` | `mean` |

Common keys:

* `scope`: `task` (default; also accepted as `environment` or `env`) or `event`. For detector
  parameters `task` means "once per environment". Anything that depends on an event-scoped
  parameter is itself event-scoped.
* `ratio_of: <parameter>`: the value is the draw times that parameter (a bare particle name means
  its mass), e.g. `mass: {dist: uniform, min: 0.1, max: 0.9, ratio_of: X}`.

Expressions support `+ - * / ^`, unary minus, parentheses and the functions `min(a, b)`,
`max(a, b)`, `sqrt`, `abs`, `log`, `exp`, `pow(a, b)`. Identifiers must match
`[A-Za-z_][A-Za-z0-9_.]*`, so parameters of particles whose names contain other characters
(`W+`, `tau-`) cannot appear in expressions; use `ratio_of`, or name such particles `Wp`,
`taum`. A bare particle name means its mass.

There are three separate parameter namespaces, which cannot reference each other: physics
(top-level `parameters`, particles, processes, `initial_state`, `parton_shower`,
`underlying_event`), `detector`, and `features`.

Unknown keys in parameter and distribution maps are errors.

## Prior schema

### Top level

| Key | Default | Meaning |
|---|---|---|
| `name` | `prior` | name recorded in tasks |
| `seed` | 0 | default master seed (used when no seed is given) |
| `parameters` | none | map name → parameter; free parameters usable in expressions |
| `particles` | required | map name → particle |
| `processes` | — | list of processes (give this or `production`) |
| `production` | — | shorthand for one process named `default` |
| `task` | see [task](#task) | how tasks are assembled from processes |
| `initial_state` | see [initial_state](#initial_state) | luminosity, rapidity and recoil (ISR) model |
| `bw_cutoff` | 15 | Breit-Wigner masses limited to M ± bw_cutoff·Γ (0: physical limits only) |
| `parton_shower` | none | final-state parton shower |
| `underlying_event` | none | soft underlying-event particles |
| `detector` | perfect detector | detector prior |
| `features` | none | auxiliary per-object features |

### Particles

```yaml
particles:
  <name>:
    pdg: 6                 # default 0
    charge: 0.667          # default 0 (not taken from pdg)
    mass: 172.5            # parameter, default 0  -> "<name>.mass"
    width: 1.42            # parameter, default 0  -> "<name>.width"
    stable: false          # default: true without decays, false with decays (a contradiction is an error)
    visible: true          # default true; false forces object "invisible"
    object: jet            # reconstructed type; default from pdg (table below)
    decay_choice: event    # event (default): mode drawn per event; task: one mode fixed per task
    normalize_br: false    # allow fixed BRs that do not sum to one
    decays:
      - daughters: [W, b]  # 2 or 3 daughters (particle names)
        br: 1.0            # parameter (task level) -> "<name>.decay<i>.br"; required with several modes, default 1 for one
        angular: isotropic # see below; -> "<name>.decay<i>.angular.c<k>"
        model: phase_space # 3-body only: phase_space (default) or legacy_sequential
        reconstruct_as: tau  # optional composite: visible descendants form one object of this type
    # shorthand for a single mode with br 1:  decays: [[W, b]]
```

Particle names are free strings (subject to the expression caveat above). Default object types
from |pdg|: 1–4 and 21 → `jet`, 5 → `bjet`, 11 → `electron`, 13 → `muon`, 15 → `tau`,
22 → `photon`, 12/14/16 → `invisible`; any other stable visible particle needs `object:`.
Unstable particles have no object type.

Angular models (`angular:` on decays and productions):

| Form | Meaning |
|---|---|
| `isotropic` (alias `flat`), or omitted | flat in cos θ* |
| `{type: polynomial, coeffs: [c0, c1, ...], axis: helicity\|beam}` | f(cos θ*) = Σ cᵢ cosⁱθ*, coefficients are task-level parameters, must be non-negative on [−1, 1] |
| `W` (alias `w_leptonic`) | f = 1 − cos θ* about the **beam** axis (myGenerator's W mode) |

θ* is the angle of the first daughter in the parent rest frame. The default axis is `helicity`
(parent flight direction) for decays and `beam` for production.

### Processes

```yaml
processes:
  - name: pair                     # default process<i>
    weight: 2                      # plain number (not a parameter), default 1: task-level sampling weight
    production: [X, X]             # 1-3 slots; a slot is a particle name or {choice: [...], weights: [...]}
    multiplicity: {values: [1, 2, 3], weights: [5, 3, 2]}   # or an integer; replicates the single slot
    angular: {type: polynomial, coeffs: [1, 0, 1]}          # production angle, default isotropic, beam axis
  - name: W
    subprocesses:                  # instead of production/multiplicity/angular
      - {production: [W], weight: 0.6}           # weight: parameter, default 1 -> "process.W.sub0.weight"
      - {production: [W, q], weight: 0.3, angular: {type: polynomial, coeffs: [1, 0, 2]}}
production: [A, B]                 # top-level shorthand (process "default"); not together with processes
```

When things are drawn:

* per task: which processes form the task's components, production `choice` slots, `multiplicity`,
  task-level parameters, decay modes when `decay_choice: task`;
* per event: the component (by mixture fraction), the subprocess (by weight), decay modes
  (default), event-scoped parameters, masses within widths, kinematics.

### task

```yaml
task:
  components: 1                      # number or distribution (rounded); default 1
  fractions: equal                   # equal (default) | [f1, f2, ...] | {dirichlet: alpha}
  distinct_processes: true           # draw component processes without replacement (default true)
  validation: reject                 # reject (default): redraw invalid tasks | error: raise at once
  max_attempts: 1000                 # redraws before giving up (ignored with validation: error)
```

Components are drawn from `processes` by `weight`. A fixed `fractions` list needs a fixed
`components` of the same length; if that length equals the number of processes, component *i* is
process *i* in catalogue order, otherwise processes are drawn by weight and the fractions are
assigned in draw order. The event label (`component`) is the component index.

### initial_state

```yaml
initial_state:
  sqrt_s:
    model: exponential_pdf   # exponential_pdf (default) | pdf | threshold_exponential | fixed | uniform
    ebeam: 6500              # GeV per beam
    tau: 0.05                # exponential_pdf: x ~ Exp(tau + tau_slope * threshold / TeV)
    tau_slope: 0.02
    threshold_offset: 0      # GeV above the produced masses (2-3 body production only)
    scale: 100               # threshold_exponential: sqrt_s = threshold + Exp(scale)
    value: 1000              # fixed
    min: 0                   # uniform
    max: 13000               # uniform
    a: 1.3                   # pdf: f(x) ~ x^-a (1 - x)^b  -> "initial_state.pdf_a"
    b: 5.0                   #                              -> "initial_state.pdf_b"
  rapidity: from_pdf         # default for exponential_pdf and pdf; otherwise a parameter (default fixed 0)
  recoil:                    # default: none (also `recoil: none`)
    model: exponential       # exponential (default) | param | none
    min_pt: 10               # exponential: pT = min_pt + Exp(mean_pt - min_pt) per jet
    mean_pt: 15
    pt: <parameter>          # param model: pT per jet (redrawn per jet if event-scoped)
    n_jets: 1                # number of recoil jets (parameter; e.g. Poisson, event scope)
    rapidity: <parameter>    # per jet; default uniform(-2.5, 2.5), event scope
    emit: false              # true: recoil jets are visible final-state objects
    object: jet              # their object type
  max_attempts: 10000000     # rejection-sampling attempts per event
```

All numbers are parameters named `initial_state.<field>` except `a`, `b` (→ `pdf_a`, `pdf_b`) and
the recoil fields (→ `recoil_min_pt`, `recoil_mean_pt`, `recoil_pt`, `recoil_n_jets`,
`recoil_rapidity`).

For a **single produced particle**, √ŝ is that particle's (Breit-Wigner) mass whatever the
`sqrt_s` model; `fixed`, `uniform` and `threshold_exponential` settings and `threshold_offset` are
then not used, and only the rapidity follows the luminosity model (`from_pdf`) or the `rapidity`
parameter.

### parton_shower

```yaml
parton_shower:           # absent, false or none: no shower; {} or true: all defaults
  alpha_s: 0.13          # effective alpha_s(mZ), one-loop running
  kt_min: 1.0            # GeV, shower cutoff
  hadron_mass: 0.5       # GeV, mass of each final shower parton
  hadron_shift: 2.0      # GeV, jet m^2 += hadron_shift * (C_i / C_F) * pT * R
  theta_max: 1.0         # rad, largest emission angle (partons without a colour partner)
```

Parameters `parton_shower.<field>`. Showered: visible final-state particles with |pdg| 1–5 or 21
(including emitted recoil jets). See DESIGN.md for the model.

### underlying_event

```yaml
underlying_event:        # absent, false or none: none; {} or true: all defaults
  density: 4.0           # GeV of pT per unit eta-phi area
  pt_mean: 1.5           # GeV per particle
  eta_max: 5.0
```

Parameters `underlying_event.<field>`. The particles only enter jet clustering (see
`detector.jets`); without a clustering collection with `jet` inputs they have no effect.

### detector

```yaml
detector:
  objects:
    <type>:                # one entry per reconstructed type, plus `default` for unlisted types
      efficiency: {plateau: 0.95, pt50: 20, width: 4}  # or a number (plateau); flat form: efficiency + eff_pt50 + eff_width
      pt_min: 0            # applied to the reconstructed type after smearing
      eta_max: .inf        # acceptance in |eta|
      resolution: {stochastic: 0, constant: 0, noise: 0, pt_term: 0}
      scale: 1
      sigma_eta: 0
      sigma_phi: 0
      misid: {<other type>: p, ...}      # probabilities, sum <= 1
      response_map: {eta_edges, pt_edges, scale, resolution, efficiency}   # optional binned response
      resolution_scale: 1  # factor on the map resolution
      mass_scale: 1        # optional mass response (any of these three enables it)
      mass_resolution: 0   # (1 with a mass_response map: a factor on its resolution)
      mass_response: {pt_edges, scale, resolution}   # optional binned mass response in true pT
  jets:                    # optional anti-kT collections
    - type: jet            # output object type (required)
      inputs: [jet]        # candidate types clustered (default [jet])
      radius: 0.4          # parameter, in (0, 3)  -> "detector.jets.<type>.radius"
      in_met: true         # true: replaces its inputs and enters MET; false: extra view outside MET
  fakes:
    rate: 0                # Poisson mean per event
    types: {jet: 1}        # type -> relative weight (parameters); required when rate > 0
    pt_min: 20             # pT = pt_min + Exp(pt_scale)
    pt_scale: 10
    eta_max: 2.5           # uniform in |eta| < eta_max, uniform phi
  correlations:
    <group>: {types: [jet, bjet], scale_sigma: 0, resolution_correlation: 0, unclustered: false}
  met: {soft_sigma: 0, sumet_coeff: 0, recoil_scale: 1, recoil_resolution: 0, lost_to_met: 1}
  sort_by_pt: true
```

Notes:

* Missing object entries mean a perfect detector for that type (efficiency 1, no smearing, no
  acceptance cut).
* An efficiency without `width` (or `width: 0`) is a step function at `pt50`.
* `misid` and fakes create objects of other types; every type that can appear is listed in
  `TaskConfig::object_types` (sorted), whose index is `type_code`.
* Fakes are massless, bypass efficiency, smearing, misidentification and `pt_min`, and enter MET
  and ΣET.
* With a parton shower, b quarks (default type `bjet`) are clustered only if `bjet` is among
  `inputs` (or the prior sets `b: {object: jet}`); composite (`reconstruct_as`) objects are never
  clustered.
* ΣET (`sum_et`) is the scalar pT sum of the observed objects that enter MET (fakes included,
  `in_met: false` collections excluded) plus the magnitude of the unclustered recoil pT.
  `lost_to_met: 1` (default) means objects that are not reconstructed are not subtracted, so their
  momentum shows up as MET; `0` subtracts them (as if measured in the soft term).

### features

```yaml
features:
  - name: <unique name>
    kind: gaussian | uniform | truth_proxy | observed_proxy | label_correlated | environment
    applies_to: [jet, ...]   # object types (default all); NaN elsewhere
    missing_prob: 0          # per-object probability of NaN
    presence_prob: 1         # probability that the feature is in a task's schema
    # kind-specific:
    mean: 0, sigma: 1                         # gaussian
    min: 0, max: 1                            # uniform
    quantity: pt, rel_sigma: 0, abs_sigma: 0, bias: 0   # truth_proxy / observed_proxy
    fake_policy: missing                      # truth_proxy: missing (NaN for fakes) | observed
    label: component, base_mean: 0, sigma: 1, shift: <distribution or number>   # label_correlated
    param: detector.jet.resolution.stochastic, coeff: 1, offset: 0, sigma: 0     # environment
```

* Proxies: value = x + bias + rel_sigma·|x|·N(0,1) + abs_sigma·N(0,1). Truth quantities: `pt eta
  phi m e px py pz charge pdg is_fake depth branch parent_pdg parent_mass`. Observed quantities:
  `pt eta phi m e px py pz charge n_constituents n_charged n_partons`. With `fake_policy: missing`
  fakes get NaN, except `quantity: is_fake`, which gives 1.
* Labels: `component`, `process`, `truth_type` (truth particle name or `fake`), `parent`,
  `is_fake`. A shift per category is drawn once per task; value = base_mean + shift + sigma·N(0,1).
* Environment features: offset + coeff · (detector parameter value, per event for event-scoped
  parameters) + sigma·N(0,1).
* Feature parameters are named `feature.<name>.<field>` and must be task-level.

## Parameter names

Inline parameters get these names (used in expressions, `ratio_of`, splits, environment features,
`TaskConfig::param` and the serialised tasks):

| Where | Name |
|---|---|
| top-level `parameters` | the key itself |
| particle | `<particle>.mass`, `<particle>.width`, `<particle>.decay<i>.br`, `<particle>.decay<i>.angular.c<k>` |
| process | `process.<name>.angular.c<k>`; with subprocesses `process.<name>.sub<i>.weight`, `process.<name>.sub<i>.angular.c<k>` |
| initial state | `initial_state.ebeam`, `.tau`, `.tau_slope`, `.threshold_offset`, `.scale`, `.value`, `.min`, `.max`, `.pdf_a`, `.pdf_b`, `.rapidity`, `.recoil_min_pt`, `.recoil_mean_pt`, `.recoil_pt`, `.recoil_n_jets`, `.recoil_rapidity` |
| shower, underlying event | `parton_shower.<field>`, `underlying_event.<field>` |
| detector objects | `detector.<type>.efficiency` (plateau), `.eff_pt50`, `.eff_width`, `.pt_min`, `.eta_max`, `.resolution.<term>`, `.scale`, `.sigma_eta`, `.sigma_phi`, `.resolution_scale`, `.misid.<type>`, `.mass_scale`, `.mass_resolution` |
| detector | `detector.jets.<type>.radius`, `detector.fakes.<field>`, `detector.fakes.types.<type>`, `detector.met.<field>`, `detector.correlations.<group>.scale_sigma`, `.resolution_correlation` |
| features | `feature.<name>.<field>` |

Process `weight` is a plain number, not a parameter.

## Validation and errors

Unknown keys are errors at every level of the prior, including parameter and distribution maps.
Static checks at parse time, per-task checks (`InvalidTask`, redrawn with `validation: reject`) and
per-event errors are listed in DESIGN.md. Nothing is clipped silently. C++ throws
`std::invalid_argument` (bad configuration), `kg::InvalidTask` and `std::runtime_error`
(generation failures); Python maps them to `ValueError`, `kinegen.InvalidTask` (a
`RuntimeError`) and `RuntimeError`. One failing event aborts a whole `generate_batch` call.

## Python API

```python
import kinegen as kg
kg.__version__
```

| Function | Returns |
|---|---|
| `load_prior(path)`, `parse_prior(yaml)` | `PriorConfig` (`name`, `seed`, `particles`, `processes`, `object_types`) |
| `prior_from_legacy_card(path, task_level_masses=False)`, `convert_legacy_card(path)` | `PriorConfig`, YAML string |
| `sample_task(prior, seed=None, index=0)` | `TaskConfig` (seed None: the prior's seed) |
| `task_from_yaml(yaml)`, `load_task(path)` | `TaskConfig` |
| `sample_environment(prior, seed, index)` | `Environment` |
| `generate_truth(task, index)` | `TruthEvent` |
| `generate_event(task, index, env=None)` | `(TruthEvent, ObservedEvent)` |
| `generate_batch(task, n, start=0, env=None, max_objects=None, threads=1, include_truth=False)` | dict of arrays (below) |
| `assign_split(task, by="task", fractions=None, salt=0, parameter="", edges=[], labels=[], env=None)` | split label |
| `DetectorModel(task).apply(truth, env=task env)` | `ObservedEvent` |

Classes (read-only attributes):

* `TaskConfig`: `task_id`, `seed`, `master_seed`, `task_index`, `attempts`, `prior_name`,
  `object_types`, `feature_names`, `environment`, `topology`, `parameters` (dict),
  `components` (list of `{process, produced, subprocesses, fraction}`), `param(name)`,
  `to_yaml()`, `from_yaml(yaml)`, `save(path)`.
* `Environment`: `env_id`, `seed`, `parameters`, `value(name)`, `to_yaml()`, `from_yaml(yaml)`.
* `TruthEvent`: `task_id`, `event_index`, `seed`, `component`, `subprocess`, `process`,
  `sqrt_s_hat`, `particles`, `shower`, `underlying_event`, `event_params`, `final_state()`,
  `ancestors(i)`, `topology()`, `to_yaml()`. `TruthParticle`: `index`, `name`, `pdg`, `charge`,
  `mass`, `parent`, `children`, `branch`, `depth`, `decay_mode`, `object_type`,
  `reconstruct_as`, `status`, `p4`. `ShowerParticle`: `parent`, `pdg`, `p4`.
* `ObservedEvent`: `task_id`, `env_id`, `event_index`, `component`, `process`, `objects`,
  `met_x`, `met_y`, `met`, `sum_et`, `n_lost`, `n_fakes`, `n_misid`, `detector_params`,
  `feature_names`, `to_yaml()`. `ObservedObject`: `type`, `type_code`, `charge`, `features`,
  `truth_index`, `true_type`, `is_fake`, `n_constituents`, `n_charged`, `n_partons`, `p4`, `pt`,
  `eta`, `phi`, `m`.

`generate_batch` generates events `[start, start + n)` with `threads` threads (the GIL is
released) and returns, with n events and M object slots (M = largest event, or `max_objects`;
events with more objects keep the first M, i.e. the highest-pT ones when `sort_by_pt` is on):

| Key | Shape, dtype | Content |
|---|---|---|
| `p4` | (n, M, 4) float32 | px, py, pz, E |
| `ptetaphim` | (n, M, 4) float32 | pT, η, φ, m |
| `type` | (n, M) int32 | index into `object_types`; −1 padding |
| `mask` | (n, M) bool | real object |
| `features` | (n, M, F) float32 | NaN = missing or padding |
| `truth_index` | (n, M) int32 | truth particle index; −1 for fakes and padding |
| `is_fake` | (n, M) bool | |
| `charge` | (n, M) float32 | |
| `n_constituents`, `n_charged`, `n_partons` | (n, M) int32 | see DESIGN.md |
| `met` | (n, 2) float32 | met_x, met_y |
| `label` | (n,) int32 | component index |
| `n_objects` | (n,) int32 | object count before truncation |
| `event_index` | (n,) uint64 | |
| `truth_met` | (n, 2) float32 | true invisible pT (targets; not model inputs) |
| `sqrt_s_hat` | (n,) float32 | |
| `max_produced_mass` | (n,) float32 | largest mass among produced particles |
| `n_visible_truth`, `n_invisible_truth` | (n,) int32 | final-state truth counts |
| `feature_names`, `object_types` | lists | |
| `truth` | list of `TruthEvent` | only with `include_truth=True` |

## Command-line tools

`kg-generate (--prior FILE | --legacy-card FILE) [options]` (built with the C++ build, not by pip):

| Option | Default | Meaning |
|---|---|---|
| `--tasks N` | 1 | tasks to sample |
| `--first-task I` | 0 | index of the first task |
| `--events M` | 1000 (legacy card: its `nevents`) | events per task |
| `--seed S` | the prior's `seed` | master seed |
| `--envs K` | 0 | observe every event in K environments sampled from the master seed instead of the task's own environment (each truth event is written K times) |
| `--out FILE.root` | none (legacy card: its `outfile`) | ROOT output; requires a build with ROOT |
| `--yaml-out FILE` | none | YAML stream of tasks and observed events (no truth); works without ROOT |
| `--legacy-format` / `--new-format` | legacy for `--legacy-card` | `test`/`info` trees or the new trees |
| `--no-truth` | | omit truth branches |
| `--save-tasks DIR` | | write each task as `DIR/<task_id>.yaml` |
| `--task-level-masses` | | legacy card: mass ranges become task-level priors |
| `--quiet` | | no progress output |

With neither `--out` nor `--yaml-out` nothing is written (useful for timing). `kg-generate` is
single-threaded. `kg-convert-card process.card [out.yaml] [--task-level-masses]` converts a legacy
card (prints to stdout without `out.yaml`).

## Output formats

ROOT (`kg::RootWriter`, `kg-generate --out`):

* `events`: `task_id`, `task_index`, `env_id`, `event_index`, `truth_seed`, `detector_seed`,
  `component` (the event label), `process`, `n_obj`, `n_feat`, `obj_px/py/pz/e/pt/eta/phi/m`,
  `obj_charge`, `obj_type` (integer index into `tasks.object_types`), `obj_truth`, `obj_fake`,
  `obj_feat` (flattened, `n_feat` per object), `met`, `met_phi`, `met_x`, `met_y`, `sum_et`,
  `n_lost`, `n_fakes`, `n_misid`, `true_met_x`, `true_met_y`, `sqrt_s_hat`, and unless
  `--no-truth` `truth_pdg/parent/status/species/branch/depth/px/py/pz/e/mass`. Not written:
  `n_constituents`, `n_charged`, `n_partons`, `true_type`, `subprocess`, shower and
  underlying-event particles (use Python or the YAML stream for these).
* `tasks`: `task_id`, `task_index`, `prior_name`, `seed`, `master_seed`, `topology`,
  `processes`, `fractions`, `feature_names`, `object_types`, `param_names`, `param_values`,
  `yaml` (full task).
* `environments`: `env_id`, `seed`, `param_names`, `param_values`, `yaml`.
* Legacy format: the original `test` (`P1*`, `P2*`, `T1*`, `T2*`, `MET*`) and `info` (`truthM*`,
  `truthId`) trees.

YAML (`--yaml-out`, `taskToYAML`, `truthEventToYAML`, `observedEventToYAML`): documents separated
by `--- # task` / `--- # event`; every value at 17 significant digits, so reading back reproduces
events exactly.

## Splits

`assignSplit(task, SplitSpec)` / `kinegen.assign_split`:

* `by`: `task` (key = task id and seed), `topology` (`TaskConfig::topology()`, which lists only
  decay modes with BR > 0), `process` (the set of component processes), `environment` (the
  environment id and seed; the task's own environment unless `env` is given), `parameter`.
* `fractions`: label → weight, default `{train: 0.8, val: 0.1, test: 0.1}`; the key is hashed with
  `salt` (change the salt for a different but equally deterministic assignment).
* `parameter` mode: `labels[i]` for `edges[i−1] ≤ value < edges[i]` of the task-level parameter
  `parameter`; needs `len(labels) == len(edges) + 1` and sorted edges.

## Using the C++ library

The CMake project builds the library target `kinegen` (static unless `BUILD_SHARED_LIBS` is on;
plus `kinegen_root` with ROOT); it has no
install rules, so use it from another CMake project with `add_subdirectory(kinematics-generator)`
and `target_link_libraries(myapp PRIVATE kinegen)` (headers under `include/kg/`, umbrella header
`kg/kinegen.h`). CMake options: `KG_BUILD_PYTHON` (OFF), `KG_BUILD_TESTS`, `KG_BUILD_EXAMPLES`,
`KG_BUILD_TOOLS` (ON), `KG_WITH_ROOT` (`AUTO`: use ROOT if found; `ON`/`OFF`).
