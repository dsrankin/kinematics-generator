#include "kg/Distribution.h"

#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace kg {

namespace {

double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

std::string ctx(const std::string &c) { return c.empty() ? "distribution" : c; }

} // namespace

double sampleRelativisticBW(Rng &rng, double mean, double width, double lo, double hi) {
  if (!(hi >= lo)) throw std::invalid_argument("sampleRelativisticBW: empty range");
  if (width <= 0) {
    if (mean < lo || mean > hi) throw std::invalid_argument("sampleRelativisticBW: zero width and mean outside range");
    return mean;
  }
  const double m2 = mean * mean, mg = mean * width;
  const double slo = lo > 0 ? lo * lo : 0.0;
  const double ulo = std::atan((slo - m2) / mg);
  const double uhi = std::isinf(hi) ? M_PI / 2 : std::atan((hi * hi - m2) / mg);
  const double s = m2 + mg * std::tan(rng.uniform(ulo, uhi));
  double m = std::sqrt(std::max(s, 0.0));
  // Guard against round-off at the edges.
  if (m < lo) m = lo;
  if (m > hi) m = hi;
  return m;
}

Distribution Distribution::fixed(double v) {
  Distribution d;
  d.type = Type::Fixed;
  d.value = v;
  return d;
}
Distribution Distribution::uniform(double lo, double hi) {
  Distribution d;
  d.type = Type::Uniform;
  d.min = lo;
  d.max = hi;
  return d;
}
Distribution Distribution::logUniform(double lo, double hi) {
  Distribution d = uniform(lo, hi);
  d.type = Type::LogUniform;
  return d;
}
Distribution Distribution::gaussian(double mu, double sd) {
  Distribution d;
  d.type = Type::Gaussian;
  d.mean = mu;
  d.sigma = sd;
  return d;
}
Distribution Distribution::truncatedGaussian(double mu, double sd, double lo, double hi) {
  Distribution d = gaussian(mu, sd);
  d.type = Type::TruncatedGaussian;
  d.min = lo;
  d.max = hi;
  return d;
}
Distribution Distribution::breitWigner(double mu, double gamma, double lo, double hi, bool rel) {
  Distribution d;
  d.type = Type::BreitWigner;
  d.mean = mu;
  d.width = gamma;
  d.min = lo;
  d.max = hi;
  d.relativistic = rel;
  return d;
}
Distribution Distribution::categorical(std::vector<double> vals, std::vector<double> w) {
  Distribution d;
  d.type = Type::Categorical;
  d.values = std::move(vals);
  d.weights = std::move(w);
  return d;
}
Distribution Distribution::exponential(double s, double off) {
  Distribution d;
  d.type = Type::Exponential;
  d.scale = s;
  d.offset = off;
  return d;
}

Distribution Distribution::poisson(double mu) {
  Distribution d;
  d.type = Type::Poisson;
  d.mean = mu;
  return d;
}

double Distribution::sample(Rng &rng) const {
  switch (type) {
    case Type::Fixed: return value;
    case Type::Uniform: return rng.uniform(min, max);
    case Type::LogUniform: return std::exp(rng.uniform(std::log(min), std::log(max)));
    case Type::Gaussian: return rng.gauss(mean, sigma);
    case Type::TruncatedGaussian: {
      const double a = (min - mean) / sigma, b = (max - mean) / sigma;
      const double pa = normCdf(a), pb = normCdf(b);
      if (pb - pa > 0.05) {
        while (true) {
          const double x = rng.gauss();
          if (x >= a && x <= b) return mean + sigma * x;
        }
      }
      // Narrow/tail window: invert the CDF by bisection.
      const double target = pa + rng.uniform() * (pb - pa);
      double lo = std::isinf(a) ? -40 : a, hi = std::isinf(b) ? 40 : b;
      for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        (normCdf(mid) < target ? lo : hi) = mid;
      }
      return mean + sigma * 0.5 * (lo + hi);
    }
    case Type::BreitWigner: {
      if (relativistic) return sampleRelativisticBW(rng, mean, width, std::max(min, 0.0), max);
      if (width <= 0) return mean;
      const double hw = 0.5 * width;
      const double ulo = std::isinf(min) ? -M_PI / 2 : std::atan((min - mean) / hw);
      const double uhi = std::isinf(max) ? M_PI / 2 : std::atan((max - mean) / hw);
      double x = mean + hw * std::tan(rng.uniform(ulo, uhi));
      return std::min(std::max(x, min), max);
    }
    case Type::Categorical: {
      if (weights.empty()) return values[static_cast<size_t>(rng.uniform() * values.size()) % values.size()];
      return values[rng.categorical(weights)];
    }
    case Type::Exponential: return offset + rng.exponential(scale);
    case Type::Poisson: return rng.poisson(mean);
  }
  return value;
}

void Distribution::validate(const std::string &c) const {
  auto fail = [&](const std::string &m) { throw std::invalid_argument(ctx(c) + ": " + m); };
  switch (type) {
    case Type::Fixed:
      if (std::isnan(value)) fail("fixed value is NaN");
      break;
    case Type::Uniform:
      if (!std::isfinite(min) || !std::isfinite(max) || !(min <= max)) fail("uniform requires finite min <= max");
      break;
    case Type::LogUniform:
      if (!(min > 0) || !std::isfinite(max) || !(min <= max)) fail("loguniform requires 0 < min <= max");
      break;
    case Type::Gaussian:
      if (!(sigma >= 0)) fail("gaussian requires sigma >= 0");
      break;
    case Type::TruncatedGaussian:
      if (!(sigma > 0)) fail("truncated_gaussian requires sigma > 0");
      if (!(min < max)) fail("truncated_gaussian requires min < max");
      break;
    case Type::BreitWigner:
      if (!(width >= 0)) fail("breit_wigner requires width >= 0");
      if (!(min <= max)) fail("breit_wigner requires min <= max");
      if (width == 0 && (mean < min || mean > max)) fail("breit_wigner with zero width has mean outside [min, max]");
      if (relativistic && !(mean > 0)) fail("relativistic breit_wigner requires mean > 0");
      break;
    case Type::Categorical: {
      if (values.empty()) fail("categorical requires at least one value");
      if (!weights.empty()) {
        if (weights.size() != values.size()) fail("categorical weights and values differ in length");
        double tot = 0;
        for (double w : weights) {
          if (!(w >= 0)) fail("categorical weights must be >= 0");
          tot += w;
        }
        if (!(tot > 0)) fail("categorical weights must have a positive sum");
      }
      break;
    }
    case Type::Exponential:
      if (!(scale > 0)) fail("exponential requires scale > 0");
      break;
    case Type::Poisson:
      if (!(mean >= 0)) fail("poisson requires mean >= 0");
      break;
  }
}

std::string Distribution::typeName() const {
  switch (type) {
    case Type::Fixed: return "fixed";
    case Type::Uniform: return "uniform";
    case Type::LogUniform: return "loguniform";
    case Type::Gaussian: return "gaussian";
    case Type::TruncatedGaussian: return "truncated_gaussian";
    case Type::BreitWigner: return "breit_wigner";
    case Type::Categorical: return "categorical";
    case Type::Exponential: return "exponential";
    case Type::Poisson: return "poisson";
  }
  return "fixed";
}

Distribution Distribution::fromYAML(const YAML::Node &n, const std::string &c, bool parameter_keys) {
  if (!n.IsMap() || !n["dist"]) throw std::invalid_argument(ctx(c) + ": expected a map with a 'dist' key");
  const std::string t = n["dist"].as<std::string>();
  // Keys a distribution map may contain (plus `scope` and `ratio_of` for parameters).
  static const std::map<std::string, std::set<std::string>> keys = {
      {"fixed", {"value"}}, {"uniform", {"min", "max"}}, {"loguniform", {"min", "max"}}, {"log_uniform", {"min", "max"}},
      {"gaussian", {"mean", "sigma"}}, {"normal", {"mean", "sigma"}}, {"truncated_gaussian", {"mean", "sigma", "min", "max"}},
      {"breit_wigner", {"mean", "width", "min", "max", "relativistic"}}, {"categorical", {"values", "weights"}},
      {"exponential", {"scale", "offset"}}, {"poisson", {"mean"}}};
  auto known = keys.find(t);
  if (known != keys.end())
    for (const auto &kv : n) {
      const std::string k = kv.first.as<std::string>();
      if (k != "dist" && !(parameter_keys && (k == "scope" || k == "ratio_of")) && !known->second.count(k))
        throw std::invalid_argument(ctx(c) + ": unknown key '" + k + "' for distribution '" + t + "'");
    }
  auto req = [&](const char *k) {
    if (!n[k]) throw std::invalid_argument(ctx(c) + ": '" + t + "' requires '" + k + "'");
    return n[k].as<double>();
  };
  auto opt = [&](const char *k, double def) { return n[k] ? n[k].as<double>() : def; };
  const double inf = std::numeric_limits<double>::infinity();
  Distribution d;
  if (t == "fixed") d = fixed(req("value"));
  else if (t == "uniform") d = uniform(req("min"), req("max"));
  else if (t == "loguniform" || t == "log_uniform") d = logUniform(req("min"), req("max"));
  else if (t == "gaussian" || t == "normal") d = gaussian(req("mean"), req("sigma"));
  else if (t == "truncated_gaussian") d = truncatedGaussian(req("mean"), req("sigma"), opt("min", -inf), opt("max", inf));
  else if (t == "breit_wigner") d = breitWigner(req("mean"), req("width"), opt("min", -inf), opt("max", inf), n["relativistic"] && n["relativistic"].as<bool>());
  else if (t == "categorical") {
    if (!n["values"] || !n["values"].IsSequence()) throw std::invalid_argument(ctx(c) + ": categorical requires a 'values' list");
    d = categorical(n["values"].as<std::vector<double>>(), n["weights"] ? n["weights"].as<std::vector<double>>() : std::vector<double>{});
  } else if (t == "exponential") d = exponential(req("scale"), opt("offset", 0));
  else if (t == "poisson") d = poisson(req("mean"));
  else throw std::invalid_argument(ctx(c) + ": unknown distribution '" + t + "'");
  d.validate(c);
  return d;
}

void Distribution::emitFields(YAML::Emitter &out) const {
  out << YAML::Key << "dist" << YAML::Value << typeName();
  auto kv = [&](const char *k, double v) { out << YAML::Key << k << YAML::Value << v; };
  switch (type) {
    case Type::Fixed: kv("value", value); break;
    case Type::Uniform:
    case Type::LogUniform: kv("min", min); kv("max", max); break;
    case Type::Gaussian: kv("mean", mean); kv("sigma", sigma); break;
    case Type::TruncatedGaussian: kv("mean", mean); kv("sigma", sigma); kv("min", min); kv("max", max); break;
    case Type::BreitWigner:
      kv("mean", mean); kv("width", width);
      if (!std::isinf(min)) kv("min", min);
      if (!std::isinf(max)) kv("max", max);
      if (relativistic) out << YAML::Key << "relativistic" << YAML::Value << true;
      break;
    case Type::Categorical:
      out << YAML::Key << "values" << YAML::Value << YAML::Flow << values;
      if (!weights.empty()) out << YAML::Key << "weights" << YAML::Value << YAML::Flow << weights;
      break;
    case Type::Exponential: kv("scale", scale); kv("offset", offset); break;
    case Type::Poisson: kv("mean", mean); break;
  }
}

} // namespace kg
