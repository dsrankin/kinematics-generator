// Convert a legacy process card (test.cc format) into a YAML prior.
//
//   kg-convert-card process.card [out.yaml] [--task-level-masses]
#include <cstring>
#include <fstream>
#include <iostream>

#include "kg/Legacy.h"

int main(int argc, char **argv) {
  std::string in, out;
  kg::LegacyOptions opt;
  for (int i = 1; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--task-level-masses")) opt.task_level_masses = true;
    else if (in.empty()) in = argv[i];
    else out = argv[i];
  }
  if (in.empty()) {
    std::cerr << "usage: kg-convert-card process.card [out.yaml] [--task-level-masses]\n";
    return 2;
  }
  try {
    const kg::LegacyCard card = kg::convertLegacyCard(in, opt);
    kg::parsePriorString(card.yaml); // validate the result
    std::string header = "# nevents " + std::to_string(card.nevents) + " and outfile " + card.outfile +
                         " are run options: pass them to kg-generate (--events/--out)\n";
    if (out.empty()) std::cout << header << card.yaml;
    else {
      std::ofstream f(out);
      f << header << card.yaml;
      std::cerr << "wrote " << out << "\n";
    }
  } catch (const std::exception &e) {
    std::cerr << "kg-convert-card: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
