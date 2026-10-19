#include <cstring>

#include "test_framework.h"

int main(int argc, char **argv) {
  int nfail = 0, nrun = 0;
  for (const auto &c : kgtest::registry()) {
    if (argc > 1 && c.name.find(argv[1]) == std::string::npos) continue;
    ++nrun;
    try {
      c.fn();
      std::cout << "[ PASS ] " << c.name << "\n";
    } catch (const std::exception &e) {
      ++nfail;
      std::cout << "[ FAIL ] " << c.name << "\n         " << e.what() << "\n";
    }
  }
  std::cout << nrun - nfail << "/" << nrun << " tests passed\n";
  return nfail ? 1 : 0;
}
