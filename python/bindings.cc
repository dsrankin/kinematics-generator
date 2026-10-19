// pybind11 bindings: `import kinegen`.
//
//   prior = kinegen.load_prior("examples/randomized_tasks.yaml")
//   task = kinegen.sample_task(prior, seed=1, index=0)
//   batch = kinegen.generate_batch(task, 4096, start=0, threads=8)   # dict of numpy arrays
#include <algorithm>
#include <cmath>
#include <limits>
#include <thread>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "kg/kinegen.h"

namespace py = pybind11;
using namespace kg;

namespace {

py::tuple p4tuple(const P4 &p) { return py::make_tuple(p.px, p.py, p.pz, p.e); }

const char *statusName(TruthStatus s) {
  switch (s) {
    case TruthStatus::Incoming: return "incoming";
    case TruthStatus::Recoil: return "recoil";
    case TruthStatus::Decayed: return "decayed";
    case TruthStatus::Final: return "final";
  }
  return "final";
}

// Generate events [start, start+n) in parallel; each event is independent and deterministic.
std::vector<GeneratedEvent> generateRange(const TaskConfig &task, const Environment &env, uint64_t start, size_t n, int threads) {
  std::vector<GeneratedEvent> out(n);
  const DetectorModel det(task);
  threads = std::max(1, std::min<int>(threads, static_cast<int>(n)));
  std::vector<std::thread> pool;
  std::vector<std::exception_ptr> errors(threads);
  for (int t = 0; t < threads; ++t)
    pool.emplace_back([&, t]() {
      try {
        for (size_t i = t; i < n; i += threads) {
          out[i].truth = generateTruthEvent(task, start + i);
          out[i].observed = det.apply(out[i].truth, env);
        }
      } catch (...) {
        errors[t] = std::current_exception();
      }
    });
  for (auto &th : pool) th.join();
  for (auto &e : errors)
    if (e) std::rethrow_exception(e);
  return out;
}

py::dict toArrays(const TaskConfig &task, const std::vector<GeneratedEvent> &evs, py::object max_objects, bool include_truth) {
  const size_t n = evs.size();
  size_t M = 0;
  for (const auto &e : evs) M = std::max(M, e.observed.objects.size());
  if (!max_objects.is_none()) M = max_objects.cast<size_t>();
  const size_t F = task.feature_names->size();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  using shape = std::vector<py::ssize_t>;
  const auto sn = static_cast<py::ssize_t>(n), sm = static_cast<py::ssize_t>(M), sf = static_cast<py::ssize_t>(F);
  py::array_t<float> p4(shape{sn, sm, 4}), kin(shape{sn, sm, 4}), feat(shape{sn, sm, sf}), met(shape{sn, 2});
  py::array_t<int32_t> type(shape{sn, sm}), truth(shape{sn, sm}), label(shape{sn}), nobj(shape{sn});
  py::array_t<bool> mask(shape{sn, sm}), fake(shape{sn, sm});
  py::array_t<uint64_t> index(shape{sn});
  auto P = p4.mutable_unchecked<3>();
  auto K = kin.mutable_unchecked<3>();
  auto Fe = feat.mutable_unchecked<3>();
  auto Me = met.mutable_unchecked<2>();
  auto T = type.mutable_unchecked<2>();
  auto Tr = truth.mutable_unchecked<2>();
  auto Ma = mask.mutable_unchecked<2>();
  auto Fk = fake.mutable_unchecked<2>();
  auto L = label.mutable_unchecked<1>();
  auto N = nobj.mutable_unchecked<1>();
  auto I = index.mutable_unchecked<1>();
  for (py::ssize_t i = 0; i < sn; ++i) {
    const ObservedEvent &o = evs[i].observed;
    L(i) = o.component;
    N(i) = static_cast<int32_t>(o.objects.size());
    I(i) = o.event_index;
    Me(i, 0) = static_cast<float>(o.met_x);
    Me(i, 1) = static_cast<float>(o.met_y);
    for (py::ssize_t j = 0; j < sm; ++j) {
      const bool have = static_cast<size_t>(j) < o.objects.size();
      Ma(i, j) = have;
      if (!have) {
        for (int k = 0; k < 4; ++k) P(i, j, k) = K(i, j, k) = 0;
        for (py::ssize_t f = 0; f < sf; ++f) Fe(i, j, f) = static_cast<float>(nan);
        T(i, j) = -1;
        Tr(i, j) = -1;
        Fk(i, j) = false;
        continue;
      }
      const ObservedObject &ob = o.objects[j];
      P(i, j, 0) = static_cast<float>(ob.p4.px);
      P(i, j, 1) = static_cast<float>(ob.p4.py);
      P(i, j, 2) = static_cast<float>(ob.p4.pz);
      P(i, j, 3) = static_cast<float>(ob.p4.e);
      K(i, j, 0) = static_cast<float>(ob.pt());
      K(i, j, 1) = static_cast<float>(ob.eta());
      K(i, j, 2) = static_cast<float>(ob.phi());
      K(i, j, 3) = static_cast<float>(ob.m());
      for (py::ssize_t f = 0; f < sf; ++f) Fe(i, j, f) = static_cast<float>(ob.features[f]);
      T(i, j) = ob.type_code;
      Tr(i, j) = ob.truth_index;
      Fk(i, j) = ob.is_fake;
    }
  }
  py::dict d;
  d["p4"] = p4;
  d["ptetaphim"] = kin;
  d["type"] = type;
  d["mask"] = mask;
  d["features"] = feat;
  d["truth_index"] = truth;
  d["is_fake"] = fake;
  d["met"] = met;
  d["label"] = label;
  d["n_objects"] = nobj; // before any max_objects truncation
  d["event_index"] = index;
  d["feature_names"] = *task.feature_names;
  d["object_types"] = task.object_types;
  if (include_truth) {
    py::list tr;
    for (const auto &e : evs) tr.append(py::cast(e.truth));
    d["truth"] = tr;
  }
  return d;
}

SplitSpec makeSplit(const std::string &by, py::object fractions, uint64_t salt, const std::string &parameter,
                    const std::vector<double> &edges, const std::vector<std::string> &labels) {
  SplitSpec s;
  s.by = splitByFromString(by);
  if (!fractions.is_none()) {
    s.fractions.clear();
    for (auto kv : fractions.cast<py::dict>()) s.fractions.emplace_back(kv.first.cast<std::string>(), kv.second.cast<double>());
  }
  s.salt = salt;
  s.parameter = parameter;
  s.edges = edges;
  s.labels = labels;
  return s;
}

} // namespace

PYBIND11_MODULE(kinegen, m) {
  m.doc() = "Prior -> Task -> Truth -> Observed collider event generator";

  py::register_exception<InvalidTask>(m, "InvalidTask", PyExc_RuntimeError);

  py::class_<PriorConfig>(m, "PriorConfig")
      .def_readonly("name", &PriorConfig::name)
      .def_readonly("seed", &PriorConfig::seed)
      .def_property_readonly("object_types", &PriorConfig::objectTypes)
      .def_property_readonly("processes", [](const PriorConfig &p) {
        std::vector<std::string> n;
        for (const auto &pr : p.processes) n.push_back(pr.name);
        return n;
      })
      .def_property_readonly("particles", [](const PriorConfig &p) {
        std::vector<std::string> n;
        for (const auto &x : p.particles) n.push_back(x.name);
        return n;
      });

  py::class_<Environment>(m, "Environment")
      .def_readonly("env_id", &Environment::env_id)
      .def_readonly("seed", &Environment::seed)
      .def_property_readonly("parameters", &Environment::parameterValues)
      .def("value", &Environment::value)
      .def("to_yaml", &environmentToYAML)
      .def_static("from_yaml", &environmentFromYAML);

  py::class_<TaskConfig>(m, "TaskConfig")
      .def_readonly("task_id", &TaskConfig::task_id)
      .def_readonly("seed", &TaskConfig::seed)
      .def_readonly("master_seed", &TaskConfig::master_seed)
      .def_readonly("task_index", &TaskConfig::task_index)
      .def_readonly("attempts", &TaskConfig::attempts)
      .def_readonly("prior_name", &TaskConfig::prior_name)
      .def_readonly("object_types", &TaskConfig::object_types)
      .def_readonly("environment", &TaskConfig::environment)
      .def_property_readonly("topology", &TaskConfig::topology)
      .def_property_readonly("feature_names", [](const TaskConfig &t) { return *t.feature_names; })
      .def_property_readonly("parameters", &TaskConfig::parameterValues)
      .def_property_readonly("components", [](const TaskConfig &t) {
        py::list out;
        for (const auto &c : t.components) {
          py::dict d;
          std::vector<std::string> prod;
          for (int s : c.produced) prod.push_back(t.particles[s].name);
          d["process"] = c.process;
          d["produced"] = prod;
          d["fraction"] = c.fraction;
          out.append(d);
        }
        return out;
      })
      .def("param", &TaskConfig::param)
      .def("to_yaml", &taskToYAML)
      .def_static("from_yaml", &taskFromYAML)
      .def("save", [](const TaskConfig &t, const std::string &p) { saveTask(t, p); });

  py::class_<TruthParticle>(m, "TruthParticle")
      .def_readonly("index", &TruthParticle::index)
      .def_readonly("name", &TruthParticle::name)
      .def_readonly("pdg", &TruthParticle::pdg)
      .def_readonly("charge", &TruthParticle::charge)
      .def_readonly("mass", &TruthParticle::mass)
      .def_readonly("parent", &TruthParticle::parent)
      .def_readonly("children", &TruthParticle::children)
      .def_readonly("branch", &TruthParticle::branch)
      .def_readonly("depth", &TruthParticle::depth)
      .def_readonly("decay_mode", &TruthParticle::decay_mode)
      .def_readonly("object_type", &TruthParticle::object_type)
      .def_property_readonly("status", [](const TruthParticle &p) { return statusName(p.status); })
      .def_property_readonly("p4", [](const TruthParticle &p) { return p4tuple(p.p4); })
      .def("__repr__", [](const TruthParticle &p) { return "<TruthParticle " + std::to_string(p.index) + " " + p.name + ">"; });

  py::class_<TruthEvent>(m, "TruthEvent")
      .def_readonly("task_id", &TruthEvent::task_id)
      .def_readonly("event_index", &TruthEvent::event_index)
      .def_readonly("seed", &TruthEvent::seed)
      .def_readonly("component", &TruthEvent::component)
      .def_readonly("process", &TruthEvent::process)
      .def_readonly("sqrt_s_hat", &TruthEvent::sqrt_s_hat)
      .def_readonly("particles", &TruthEvent::particles)
      .def_readonly("event_params", &TruthEvent::event_params)
      .def("final_state", &TruthEvent::finalState)
      .def("ancestors", &TruthEvent::ancestors)
      .def("topology", &TruthEvent::topology)
      .def("to_yaml", &truthEventToYAML);

  py::class_<ObservedObject>(m, "ObservedObject")
      .def_readonly("type", &ObservedObject::type)
      .def_readonly("type_code", &ObservedObject::type_code)
      .def_readonly("charge", &ObservedObject::charge)
      .def_readonly("features", &ObservedObject::features)
      .def_readonly("truth_index", &ObservedObject::truth_index)
      .def_readonly("true_type", &ObservedObject::true_type)
      .def_readonly("is_fake", &ObservedObject::is_fake)
      .def_property_readonly("p4", [](const ObservedObject &o) { return p4tuple(o.p4); })
      .def_property_readonly("pt", &ObservedObject::pt)
      .def_property_readonly("eta", &ObservedObject::eta)
      .def_property_readonly("phi", &ObservedObject::phi)
      .def_property_readonly("m", &ObservedObject::m);

  py::class_<ObservedEvent>(m, "ObservedEvent")
      .def_readonly("task_id", &ObservedEvent::task_id)
      .def_readonly("env_id", &ObservedEvent::env_id)
      .def_readonly("event_index", &ObservedEvent::event_index)
      .def_readonly("component", &ObservedEvent::component)
      .def_readonly("process", &ObservedEvent::process)
      .def_readonly("objects", &ObservedEvent::objects)
      .def_readonly("met_x", &ObservedEvent::met_x)
      .def_readonly("met_y", &ObservedEvent::met_y)
      .def_readonly("n_lost", &ObservedEvent::n_lost)
      .def_readonly("n_fakes", &ObservedEvent::n_fakes)
      .def_readonly("detector_params", &ObservedEvent::detector_params)
      .def_property_readonly("met", &ObservedEvent::met)
      .def_property_readonly("feature_names", [](const ObservedEvent &e) {
        return e.feature_names ? *e.feature_names : std::vector<std::string>{};
      })
      .def("to_yaml", &observedEventToYAML);

  py::class_<DetectorModel>(m, "DetectorModel")
      .def(py::init<const TaskConfig &>(), py::keep_alive<1, 2>())
      .def("apply", [](const DetectorModel &d, const TruthEvent &t, const Environment &e) { return d.apply(t, e); })
      .def("apply", [](const DetectorModel &d, const TruthEvent &t) { return d.apply(t); });

  m.def("load_prior", &loadPrior, py::arg("path"));
  m.def("parse_prior", &parsePriorString, py::arg("yaml"));
  m.def("prior_from_legacy_card", [](const std::string &path, bool task_level) {
    LegacyOptions o;
    o.task_level_masses = task_level;
    return priorFromLegacyCard(path, o);
  }, py::arg("path"), py::arg("task_level_masses") = false);
  m.def("convert_legacy_card", [](const std::string &path) { return convertLegacyCard(path).yaml; }, py::arg("path"));
  m.def("sample_task", [](const PriorConfig &p, py::object seed, int64_t index) {
    return sampleTask(p, seed.is_none() ? p.seed : seed.cast<uint64_t>(), index);
  }, py::arg("prior"), py::arg("seed") = py::none(), py::arg("index") = 0);
  m.def("task_from_yaml", &taskFromYAML);
  m.def("load_task", &loadTask);
  m.def("sample_environment", [](const PriorConfig &p, uint64_t seed, int64_t index) {
    return sampleEnvironment(p.detector, seed, index);
  }, py::arg("prior"), py::arg("seed"), py::arg("index"));
  m.def("generate_truth", [](const TaskConfig &t, uint64_t i) { return generateTruthEvent(t, i); }, py::arg("task"), py::arg("index"));
  m.def("generate_event", [](const TaskConfig &t, uint64_t i, const Environment *env) {
    const GeneratedEvent g = generateEvent(t, env ? *env : t.environment, i);
    return py::make_tuple(g.truth, g.observed);
  }, py::arg("task"), py::arg("index"), py::arg("env") = nullptr);
  m.def("generate_batch", [](const TaskConfig &t, size_t n, uint64_t start, const Environment *env, py::object max_objects,
                             int threads, bool include_truth) {
    std::vector<GeneratedEvent> evs;
    {
      py::gil_scoped_release release;
      evs = generateRange(t, env ? *env : t.environment, start, n, threads);
    }
    return toArrays(t, evs, max_objects, include_truth);
  }, py::arg("task"), py::arg("n"), py::arg("start") = 0, py::arg("env") = nullptr, py::arg("max_objects") = py::none(),
        py::arg("threads") = 1, py::arg("include_truth") = false,
        "Generate events [start, start+n) and return padded numpy arrays.\n"
        "max_objects=None pads to the largest event; otherwise events are truncated to the\n"
        "max_objects highest-pT objects (n_objects keeps the untruncated count).");
  m.def("assign_split", [](const TaskConfig &t, const std::string &by, py::object fractions, uint64_t salt,
                           const std::string &parameter, const std::vector<double> &edges, const std::vector<std::string> &labels,
                           const Environment *env) {
    return assignSplit(t, makeSplit(by, fractions, salt, parameter, edges, labels), env);
  }, py::arg("task"), py::arg("by") = "task", py::arg("fractions") = py::none(), py::arg("salt") = 0, py::arg("parameter") = "",
        py::arg("edges") = std::vector<double>{}, py::arg("labels") = std::vector<std::string>{}, py::arg("env") = nullptr);
}
