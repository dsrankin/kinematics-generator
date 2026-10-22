## Kinematic Generator

A lightweight collider event generator for synthetic pretraining and
meta-learning of collider-physics models (hard process with explicit decay
trees, optional final-state parton shower, underlying event and jet clustering,
parametrised or measured detector response):

**Prior → Task → Truth Event → Observed Event**

* A **prior** (YAML) defines distributions over particles, masses, widths,
  branching fractions, angular distributions, production topologies,
  multiplicities and per-event production alternatives (hard extra partons),
  the parton luminosity and ISR recoil, the parton shower and underlying event,
  detector environments and auxiliary features.
* A **task** is one sampled physics configuration; any number of independent
  events can be generated from it, in any order, reproducibly.
* A **truth event** keeps the full decay tree (four-vectors, PDG ids, parents,
  children, ancestry).
* An **observed event** is a variable-length list of reconstructed objects
  (four-vector, type, arbitrary named features, optional truth association)
  plus MET, produced by a configurable detector model (parametric or binned
  measured response, efficiencies, fakes, misidentification, correlated
  energy-scale shifts, MET). Composite reconstruction turns e.g. hadronic tau
  decays into single tau objects; anti-kT clustering of shower partons gives
  jets with masses (QCD jets and merged boosted W/Z/H/top).

Everything runs in memory (C++ or Python); ROOT and YAML output are optional.
[DESIGN.md](DESIGN.md) explains the data flow and the physics and detector
models; [REFERENCE.md](REFERENCE.md) is the complete reference (every YAML key
with its default, parameter names, the Python API, command-line options and
output formats).

### Install the Python module (any machine)

```
pip install "git+https://github.com/dsrankin/kinematics-generator@foundation-model-generator"
# or, from a checkout:  pip install ./kinematics-generator
```

This compiles the extension for the Python running `pip` (needs a C++17 compiler, git and network
access for yaml-cpp/pybind11; CMake is fetched from PyPI if the system one is missing or too old).
It installs only the Python module `kinegen`: the command-line tools, ROOT output and the C++
library come from the CMake build below.

### Build

Requirements: a C++17 compiler and CMake ≥ 3.16. yaml-cpp and pybind11 are used
from the system if found and fetched otherwise. ROOT is optional: it is needed
for ROOT output (`kg-generate --out`, including `--legacy-card` runs, which
default to the card's ROOT file) and the legacy `test.cc`; without ROOT use
`--yaml-out`. Build options and linking the library from your own CMake project
are in [REFERENCE.md](REFERENCE.md#using-the-c-library).

```
cmake -S . -B build -DKG_BUILD_PYTHON=ON      # -DKG_WITH_ROOT=OFF to skip ROOT
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Single-threaded throughput in memory (`kg-generate --events 100000` without output, Apple M2,
Release build):

| Example | Events/s |
|---|---|
| `fixed_process.yaml` | 2.5×10⁵ |
| `detector_randomization.yaml` | 2.2×10⁵ |
| `randomized_tasks.yaml` | 1.5×10⁵ |
| `collide_tuned_sm.yaml` (shower, underlying event, jet clustering) | 3×10⁴ |
| `jet_substructure.yaml` (shower, underlying event, two jet collections) | 2×10⁴ |

The parton shower, underlying event and clustering cost roughly a factor 10. `kg-generate` is
single-threaded; Python's `generate_batch(..., threads=N)` scales with cores (about 1×10⁵ events/s
on 8 threads for a COLLIDE-like prior with shower). ROOT output with full truth runs at about
5×10⁴ events per second.

### Quick start

C++ (no intermediate files):

```cpp
#include "kg/kinegen.h"

kg::PriorConfig prior = kg::loadPrior("examples/randomized_tasks.yaml");
kg::TaskConfig task = kg::sampleTask(prior, /*seed=*/1, /*task_index=*/0);
kg::DetectorModel detector(task);
for (uint64_t i = 0; i < 1000; ++i) {
  kg::TruthEvent truth = kg::generateTruthEvent(task, i);
  kg::ObservedEvent reco = detector.apply(truth, task.environment);
}
```

See `examples/stream_example.cc` (built as `build/example_stream`), which also
observes each event through several randomised detector environments.

Python:

```python
import kinegen  # PYTHONPATH=build

prior = kinegen.load_prior("examples/randomized_tasks.yaml")
task = kinegen.sample_task(prior, seed=1, index=0)
env = kinegen.sample_environment(prior, seed=2, index=0)
batch = kinegen.generate_batch(task, 4096, start=0, env=env, threads=8)
batch["ptetaphim"], batch["features"], batch["mask"], batch["label"]   # padded numpy arrays
```

See `examples/train_feed.py` for an episodic (support/query) training feed with
topology-based train/val/test splits.

Command line (`--out` needs a ROOT build; `--yaml-out` works everywhere; options in
[REFERENCE.md](REFERENCE.md#command-line-tools)):

```
build/kg-generate --prior examples/fixed_process.yaml --events 100000 --out ttbar.root
build/kg-generate --prior examples/randomized_tasks.yaml --tasks 100 --events 2000 --seed 1 --out tasks.root
build/kg-generate --prior examples/detector_randomization.yaml --envs 10 --events 10000 --out dr.root
build/kg-generate --prior examples/randomized_tasks.yaml --tasks 3 --save-tasks tasks/   # task YAML files
build/kg-generate --prior examples/jet_substructure.yaml --events 1000 --yaml-out jets.yaml # no ROOT needed
```

### Examples

| File | Shows |
|---|---|
| `examples/fixed_process.yaml` | fixed-process generation: t t̄ with widths, W branching fractions, helicity angular distributions, Poisson multi-jet ISR, a realistic detector with a correlated jet energy scale |
| `examples/randomized_tasks.yaml` | randomised tasks: log-uniform masses, daughter/parent mass ratios, widths as fractions of masses, per-task decay-mode choice, single/pair/1–3 multiplicity production, two-process mixtures, randomised features |
| `examples/detector_randomization.yaml` | detector-domain randomisation: one physics process, environment-level and event-level detector parameters, randomised correlated jet energy scale |
| `examples/collide_tuned_sm.yaml` | realistic SM defaults: W, Z, tt̄, dibosons, QCD, γ+jet, V+γ with hard extra partons, PDF luminosity, ISR, parton shower, underlying event, anti-kT jets and a detector tuned to COLLIDE-1M (a Delphes simulation with pile-up; kinegen itself has no pile-up beyond the MET soft term). Each task is one SM process (`task: {components: 1}`); see the file header for what was and was not tuned |
| `examples/jet_substructure.yaml` | jet masses: final-state parton shower, underlying event, anti-kT R = 0.4 jets plus an R = 0.8 large-radius collection, boosted W/Z/H/top from a heavy resonance and QCD dijets, jet-mass response |
| `examples/migrated_process.yaml` | `process.card` converted with `kg-convert-card` |
| `examples/stream_example.cc`, `examples/train_feed.py` | in-memory streaming from C++ and Python |

### Migrating a process card

The original card format is still supported:

```
production A B
particle B 165.0 180.0 9998 0 0 w b
particle w 80.4 80.4 23 0 0 l v
...
```

Convert it once and edit the YAML from there:

```
build/kg-convert-card process.card examples/migrated_process.yaml    # --task-level-masses to sample ranges per task
```

which yields, e.g.,

```yaml
production: [A, B]
initial_state:
  sqrt_s: {model: exponential_pdf, ebeam: 6500, tau: 0.05, tau_slope: 0.02, threshold_offset: 1}
  rapidity: from_pdf
  recoil: {model: exponential, min_pt: 10, mean_pt: 15}
particles:
  A:
    pdg: 9999
    mass: 5.5
    decays:
      - {daughters: [q, b, q], model: legacy_sequential}   # legacy "too heavy" rule made explicit
  B:
    pdg: 9998
    mass: {dist: uniform, min: 165, max: 180, scope: event}
    decays:
      - {daughters: [w, b]}
  v: {pdg: 12, mass: 0, visible: false}
```

or run a card directly, producing the original `test`/`info` ROOT trees
(`P1*`, `P2*`, `T1*`, `T2*`, `MET*`, `truthM*`):

```
build/kg-generate --legacy-card process.card                 # nevents/outfile from the card
build/kg-generate --legacy-card process.card --new-format --out new.root
```

The converted prior reproduces the legacy kinematics: on 2×10⁵ events the means
and RMS of every `test` branch agree with `test.cc` within statistical errors.
Single-character names are no longer required; particles can be called `W+`,
`chi_10`, etc. (names used inside expressions must be identifiers, `[A-Za-z_][A-Za-z0-9_.]*`, so
prefer `Wp` over `W+` for particles whose mass you want to reference in an expression). The old mass-grid script `crank.sh` can be replaced by a prior
with mass distributions (or a `ratio_of` constraint) and `--tasks N`.

### Legacy macro

`myGenerator.h`, `test.cc` and `process.card` are kept unchanged:

```
root -l test.cc
root -l -q 'test.cc("process.card")'
```

Card syntax:
```
nevents <int>
outfile <root_file_name>
production <particleA> <particleB>
particle <name> <mass_min> <mass_max> <pid> <stable(0|1)> <visible(0|1)> <decay1> <decay2>
```

If `mass_min == mass_max`, the mass is fixed; otherwise it is sampled per event
from a uniform prior in `[mass_min, mass_max)`. For stable particles, `decay1`
and `decay2` are ignored, but placeholders are required.

### Tests

`ctest` runs the C++ suite (decay kinematics and conservation, recursive
ancestry, Dalitz flatness, angular models, parameter distributions and mass-ratio
constraints, graph and threshold validation, variable multiplicities, detector
smearing/efficiency/fakes/misidentification/MET, multi-jet ISR, correlated
smearing, composite (hadronic-tau) reconstruction, features, deterministic
reproducibility, task/environment/event serialisation round trips, splits, the
legacy adapter, parton shower, jet clustering and masses, underlying event, response maps,
validation of misspelled keys), the ROOT suite (writer, legacy trees, comparison with
`myGenerator.h`; only when built with ROOT) and the Python suite (only with
`-DKG_BUILD_PYTHON=ON` and NumPy installed).
