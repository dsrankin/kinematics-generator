// Command-line driver: sample tasks from a prior, generate events, apply
// detector environments and write ROOT (or YAML) output.
//
//   kg-generate --prior examples/randomized_tasks.yaml --tasks 10 --events 1000 --seed 1 --out out.root
//   kg-generate --prior examples/detector_randomization.yaml --envs 5 --out dr.root
//   kg-generate --legacy-card process.card            # old test.cc output format
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "kg/kinegen.h"
#ifdef KG_HAVE_ROOT
#include "kg/RootWriter.h"
#endif

namespace {

void usage() {
  std::cerr << "usage: kg-generate (--prior FILE | --legacy-card FILE) [options]\n"
               "  --tasks N           number of tasks to sample (default 1)\n"
               "  --first-task I      index of the first task (default 0)\n"
               "  --events M          events per task (default 1000; legacy card: its nevents)\n"
               "  --seed S            master seed (default: prior seed)\n"
               "  --envs K            observe every event in K sampled environments (default 0: task environment)\n"
               "  --out FILE.root     ROOT output (legacy card: its outfile)\n"
               "  --yaml-out FILE     YAML stream of tasks and observed events\n"
               "  --legacy-format     write test.cc trees (test/info); default for --legacy-card\n"
               "  --new-format        with --legacy-card, write the new event format\n"
               "  --no-truth          omit truth branches\n"
               "  --save-tasks DIR    write each task as DIR/<task_id>.yaml\n"
               "  --task-level-masses legacy card: mass ranges become task-level priors\n"
               "  --quiet\n";
}

} // namespace

int main(int argc, char **argv) {
  std::string priorPath, cardPath, out, yamlOut, saveDir;
  long long nTasks = 1, firstTask = 0, nEvents = -1, nEnvs = 0;
  bool haveSeed = false, legacyFormat = false, newFormat = false, noTruth = false, quiet = false;
  uint64_t seed = 0;
  kg::LegacyOptions lopt;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        usage();
        std::exit(2);
      }
      return argv[++i];
    };
    if (a == "--prior") priorPath = next();
    else if (a == "--legacy-card") cardPath = next();
    else if (a == "--tasks") nTasks = std::stoll(next());
    else if (a == "--first-task") firstTask = std::stoll(next());
    else if (a == "--events") nEvents = std::stoll(next());
    else if (a == "--seed") { seed = std::stoull(next()); haveSeed = true; }
    else if (a == "--envs") nEnvs = std::stoll(next());
    else if (a == "--out") out = next();
    else if (a == "--yaml-out") yamlOut = next();
    else if (a == "--legacy-format") legacyFormat = true;
    else if (a == "--new-format") newFormat = true;
    else if (a == "--no-truth") noTruth = true;
    else if (a == "--save-tasks") saveDir = next();
    else if (a == "--task-level-masses") lopt.task_level_masses = true;
    else if (a == "--quiet") quiet = true;
    else { usage(); return 2; }
  }
  if (priorPath.empty() == cardPath.empty()) { usage(); return 2; }
  try {
    kg::PriorConfig prior;
    if (!cardPath.empty()) {
      const kg::LegacyCard card = kg::convertLegacyCard(cardPath, lopt);
      prior = kg::parsePriorString(card.yaml);
      if (nEvents < 0) nEvents = card.nevents;
      if (out.empty() && yamlOut.empty()) out = card.outfile;
      if (!newFormat) legacyFormat = true;
    } else {
      prior = kg::loadPrior(priorPath);
    }
    if (nEvents < 0) nEvents = 1000;
    if (!haveSeed) seed = prior.seed;

#ifdef KG_HAVE_ROOT
    std::unique_ptr<kg::RootWriter> writer;
    if (!out.empty()) {
      kg::RootWriterOptions ro;
      ro.legacy_format = legacyFormat;
      ro.write_truth = !noTruth;
      writer.reset(new kg::RootWriter(out, ro));
    }
#else
    if (!out.empty()) throw std::runtime_error("built without ROOT: use --yaml-out");
#endif
    std::ofstream yaml;
    if (!yamlOut.empty()) yaml.open(yamlOut);

    if (!saveDir.empty()) std::filesystem::create_directories(saveDir);
    std::vector<kg::Environment> envs;
    for (long long k = 0; k < nEnvs; ++k) envs.push_back(kg::sampleEnvironment(prior.detector, seed, static_cast<int64_t>(k)));

    const auto t0 = std::chrono::steady_clock::now();
    long long nWritten = 0;
    for (long long it = firstTask; it < firstTask + nTasks; ++it) {
      const kg::TaskConfig task = kg::sampleTask(prior, seed, it);
      if (!quiet) std::cout << task.task_id << " (" << task.attempts << " draw(s)): " << task.topology() << "\n";
      if (!saveDir.empty()) kg::saveTask(task, saveDir + "/" + task.task_id + ".yaml");
      if (yaml) yaml << "--- # task\n" << kg::taskToYAML(task);
      const kg::DetectorModel detector(task);
      std::vector<const kg::Environment *> use;
      if (envs.empty()) use.push_back(&task.environment);
      for (const auto &e : envs) use.push_back(&e);
#ifdef KG_HAVE_ROOT
      if (writer) {
        writer->writeTask(task);
        for (const auto *e : use) writer->writeEnvironment(*e);
      }
#endif
      for (long long i = 0; i < nEvents; ++i) {
        const kg::TruthEvent truth = kg::generateTruthEvent(task, static_cast<uint64_t>(i));
        for (const auto *e : use) {
          const kg::ObservedEvent reco = detector.apply(truth, *e);
#ifdef KG_HAVE_ROOT
          if (writer) writer->fill(task, truth, reco);
#endif
          if (yaml) yaml << "--- # event\n" << kg::observedEventToYAML(reco);
          ++nWritten;
        }
        if (!quiet && nEvents >= 100000 && i % 100000 == 0) std::cout << "  " << i << "/" << nEvents << "\n";
      }
    }
#ifdef KG_HAVE_ROOT
    if (writer) writer->close();
#endif
    const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!quiet)
      std::cout << "wrote " << nWritten << " observed events in " << dt << " s (" << (dt > 0 ? nWritten / dt : 0) << " /s)"
                << (out.empty() ? "" : " to " + out) << "\n";
  } catch (const std::exception &e) {
    std::cerr << "kg-generate: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
