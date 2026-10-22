// Minimal self-registering test harness.
#ifndef KG_TEST_FRAMEWORK_H
#define KG_TEST_FRAMEWORK_H

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace kgtest {

struct Case {
  std::string name;
  std::function<void()> fn;
};

inline std::vector<Case> &registry() {
  static std::vector<Case> r;
  return r;
}

struct Registrar {
  Registrar(const char *name, std::function<void()> fn) { registry().push_back({name, std::move(fn)}); }
};

struct Failure : std::runtime_error {
  using std::runtime_error::runtime_error;
};

inline void fail(const std::string &msg, const char *file, int line) {
  std::ostringstream os;
  os << file << ":" << line << ": " << msg;
  throw Failure(os.str());
}

} // namespace kgtest

#define KG_CAT2(a, b) a##b
#define KG_CAT(a, b) KG_CAT2(a, b)
#define TEST(name)                                                       \
  static void name();                                                    \
  static kgtest::Registrar KG_CAT(reg_, name)(#name, name);              \
  static void name()

#define CHECK(cond)                                                      \
  do {                                                                   \
    if (!(cond)) kgtest::fail("CHECK(" #cond ") failed", __FILE__, __LINE__); \
  } while (0)

#define CHECK_MSG(cond, msg)                                             \
  do {                                                                   \
    if (!(cond)) {                                                       \
      std::ostringstream os_;                                            \
      os_ << "CHECK(" #cond ") failed: " << msg;                         \
      kgtest::fail(os_.str(), __FILE__, __LINE__);                       \
    }                                                                    \
  } while (0)

#define CHECK_CLOSE(a, b, tol)                                           \
  do {                                                                   \
    const double a_ = (a), b_ = (b), t_ = (tol);                         \
    if (!(std::fabs(a_ - b_) <= t_)) {                                   \
      std::ostringstream os_;                                            \
      os_ << "CHECK_CLOSE(" #a ", " #b ") failed: " << a_ << " vs " << b_ << " (tol " << t_ << ")"; \
      kgtest::fail(os_.str(), __FILE__, __LINE__);                       \
    }                                                                    \
  } while (0)

#define CHECK_THROWS_AS(expr, type, substr)                              \
  do {                                                                   \
    bool thrown_ = false;                                                \
    try {                                                                \
      (void)(expr);                                                      \
    } catch (const type &e_) {                                           \
      thrown_ = true;                                                    \
      if (std::string(e_.what()).find(substr) == std::string::npos)      \
        kgtest::fail(std::string("exception message '") + e_.what() + "' lacks '" + substr + "'", __FILE__, __LINE__); \
    }                                                                    \
    if (!thrown_) kgtest::fail("expected " #type " from " #expr, __FILE__, __LINE__); \
  } while (0)

#endif
