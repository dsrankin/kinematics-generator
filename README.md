## Kinematic Generator

A lightweight parton-level collider event generator for synthetic pretraining
and meta-learning of collider-physics models:

**Prior → Task → Truth Event → Observed Event**

* A **prior** (YAML) defines distributions over particles, masses, widths,
  branching fractions, angular distributions, production topologies and
  multiplicities, initial-state models, detector environments and auxiliary
  features.
* A **task** is one sampled physics configuration; any number of independent
  events can be generated from it, in any order, reproducibly.
* A **truth event** keeps the full decay tree (four-vectors, PDG ids, parents,
  children, ancestry).
* An **observed event** is a variable-length list of reconstructed objects
  (four-vector, type, arbitrary named features, optional truth association)
  plus MET, produced by a configurable detector model (smearing, efficiencies,
  fakes, misidentification, correlated energy-scale shifts, MET). Composite
  reconstruction turns e.g. hadronic tau decays into single tau objects.

Everything runs in memory (C++ or Python); ROOT and YAML output are optional.
See [DESIGN.md](DESIGN.md) for the API, data flow and full configuration
schema.

### Install the Python module (any machine)

```
pip install "git+https://github.com/dsrankin/kinematics-generator@foundation-model-generator"
# or, from a checkout:  pip install ./kinematics-generator
```

This compiles the extension for the Python running `pip` (needs a C++17 compiler, git and network
access for yaml-cpp/pybind11; CMake is fetched from PyPI if the system one is missing or too old).

### Build

Requirements: a C++17 compiler and CMake ≥ 3.16. yaml-cpp and pybind11 are used
from the system if found and fetched otherwise. ROOT is optional (needed only
for ROOT output and the legacy `test.cc`).

```
cmake -S . -B build -DKG_BUILD_PYTHON=ON      # -DKG_WITH_ROOT=OFF to skip ROOT
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Single-threaded throughput is roughly 2.5–4×10⁵ observed events per second
in memory (Apple M-series, examples below); ROOT output with full truth runs at
about 5×10⁴ events per second.

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

Command line:

```
build/kg-generate --prior examples/fixed_process.yaml --events 100000 --out ttbar.root
build/kg-generate --prior examples/randomized_tasks.yaml --tasks 100 --events 2000 --seed 1 --out tasks.root
build/kg-generate --prior examples/detector_randomization.yaml --envs 10 --events 10000 --out dr.root
build/kg-generate --prior examples/randomized_tasks.yaml --tasks 3 --save-tasks tasks/   # task YAML files
```

### Examples

| File | Shows |
|---|---|
| `examples/fixed_process.yaml` | fixed-process generation: t t̄ with widths, W branching fractions, helicity angular distributions, Poisson multi-jet ISR, a realistic detector with a correlated jet energy scale |
| `examples/randomized_tasks.yaml` | randomised tasks: log-uniform masses, daughter/parent mass ratios, widths as fractions of masses, per-task decay-mode choice, single/pair/1–3 multiplicity production, two-process mixtures, randomised features |
| `examples/detector_randomization.yaml` | detector-domain randomisation: one physics process, environment-level and event-level detector parameters, randomised correlated jet energy scale |
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
`chi_10`, etc. The old mass-grid script `crank.sh` can be replaced by a prior
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
legacy adapter), the ROOT suite (writer, legacy trees, comparison with
`myGenerator.h`) and the Python suite.
