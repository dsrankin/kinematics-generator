// Portable, seedable random numbers.
//
// All sampling in kinegen goes through kg::Rng so that results are bit-for-bit
// reproducible across platforms and standard libraries (std::*_distribution is
// implementation defined and therefore not used).
#ifndef KG_RANDOM_H
#define KG_RANDOM_H

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace kg {

// Independent random streams derived from one seed.
enum Stream : uint64_t {
  kStreamTask = 0x7461736bULL,      // "task"
  kStreamEnvironment = 0x656e76ULL, // "env"
  kStreamEvent = 0x65766eULL,       // "evn"
  kStreamDetector = 0x646574ULL,    // "det"
  kStreamFeature = 0x666561ULL,     // "fea"
  kStreamShower = 0x736872ULL,      // "shr": parton shower (truth)
  kStreamOverlay = 0x6f766cULL,     // "ovl": extra jet collections that do not enter MET
  kStreamMass = 0x6d6173ULL,        // "mas": jet-mass smearing
};

inline uint64_t splitmix64(uint64_t &x) {
  uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

// Combine seeds into a new, statistically independent seed.
inline uint64_t mixSeed(uint64_t a, uint64_t b) {
  uint64_t s = a;
  uint64_t h = splitmix64(s);
  s = h ^ (b * 0xd1b54a32d192ed03ULL + 0x8bb84b93962eacc9ULL);
  return splitmix64(s);
}
inline uint64_t mixSeed(uint64_t a, uint64_t b, uint64_t c) { return mixSeed(mixSeed(a, b), c); }

// FNV-1a, used to hash identifiers (e.g. for split assignment).
inline uint64_t hashString(const std::string &s) {
  uint64_t h = 0xcbf29ce484222325ULL;
  for (unsigned char c : s) {
    h ^= c;
    h *= 0x100000001b3ULL;
  }
  return h;
}

// xoshiro256** generator.
class Rng {
 public:
  explicit Rng(uint64_t seed = 0) { reseed(seed); }

  void reseed(uint64_t seed) {
    uint64_t x = seed;
    for (auto &s : s_) s = splitmix64(x);
  }

  uint64_t next() {
    const uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

  // Uniform in [0, 1).
  double uniform() { return (next() >> 11) * 0x1.0p-53; }
  double uniform(double a, double b) { return a + (b - a) * uniform(); }

  // Standard normal (Box-Muller, no cached second value).
  double gauss() {
    const double u1 = 1.0 - uniform(); // (0, 1]
    const double u2 = uniform();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
  }
  double gauss(double mean, double sigma) { return mean + sigma * gauss(); }

  double exponential(double mean) { return -mean * std::log(1.0 - uniform()); }

  bool bernoulli(double p) { return uniform() < p; }

  int poisson(double mean) {
    if (mean <= 0) return 0;
    if (mean > 500) {
      const double x = std::floor(gauss(mean, std::sqrt(mean)) + 0.5);
      return x < 0 ? 0 : static_cast<int>(x);
    }
    // Count unit-rate exponential arrivals within `mean`.
    int n = 0;
    double t = exponential(1.0);
    while (t < mean) {
      ++n;
      t += exponential(1.0);
    }
    return n;
  }

  // Gamma(shape, 1) via Marsaglia-Tsang.
  double gamma(double shape) {
    if (shape <= 0) throw std::invalid_argument("Rng::gamma: shape must be > 0");
    if (shape < 1.0) {
      const double u = 1.0 - uniform();
      return gamma(shape + 1.0) * std::pow(u, 1.0 / shape);
    }
    const double d = shape - 1.0 / 3.0;
    const double c = 1.0 / std::sqrt(9.0 * d);
    while (true) {
      double x, v;
      do {
        x = gauss();
        v = 1.0 + c * x;
      } while (v <= 0);
      v = v * v * v;
      const double u = 1.0 - uniform();
      if (std::log(u) < 0.5 * x * x + d - d * v + d * std::log(v)) return d * v;
    }
  }

  // Index drawn with probability proportional to weights (must be >= 0, sum > 0).
  size_t categorical(const std::vector<double> &w) {
    double tot = 0;
    for (double x : w) tot += x;
    if (!(tot > 0)) throw std::invalid_argument("Rng::categorical: weights must have a positive sum");
    double r = uniform() * tot;
    for (size_t i = 0; i < w.size(); ++i) {
      if (r < w[i]) return i;
      r -= w[i];
    }
    for (size_t i = w.size(); i-- > 0;)
      if (w[i] > 0) return i;
    return w.size() - 1;
  }

 private:
  static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  uint64_t s_[4];
};

} // namespace kg

#endif
