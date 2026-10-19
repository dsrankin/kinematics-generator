#include "kg/Legacy.h"

#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <yaml-cpp/yaml.h>

namespace kg {

namespace {

struct CardParticle {
  std::string name;
  double mmin = 0, mmax = 0;
  int pid = 0;
  bool stable = true, visible = true;
  std::string d1, d2;
  int line = 0;
};

} // namespace

LegacyCard convertLegacyCardString(const std::string &text, const LegacyOptions &opt) {
  LegacyCard card;
  std::vector<CardParticle> parts;
  std::map<std::string, size_t> index;
  std::string prodA = "A", prodB = "B";
  std::istringstream in(text);
  std::string line;
  int lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::string content = line.substr(0, line.find('#'));
    if (content.find_first_not_of(" \t\r\n") == std::string::npos) continue;
    std::istringstream iss(content);
    std::string key;
    iss >> key;
    const std::string where = " at line " + std::to_string(lineNo);
    if (key == "nevents") {
      iss >> card.nevents;
    } else if (key == "outfile") {
      iss >> card.outfile;
    } else if (key == "production") {
      iss >> prodA >> prodB;
      if (iss.fail()) throw std::invalid_argument("invalid production definition" + where);
    } else if (key == "particle") {
      CardParticle p;
      int stable = 0, visible = 0;
      iss >> p.name >> p.mmin >> p.mmax >> p.pid >> stable >> visible >> p.d1 >> p.d2;
      if (iss.fail()) throw std::invalid_argument("invalid particle definition" + where);
      if (p.mmax < p.mmin) throw std::invalid_argument("particle '" + p.name + "' has mass_max < mass_min" + where);
      if (index.count(p.name)) throw std::invalid_argument("particle '" + p.name + "' defined twice" + where);
      p.stable = stable != 0;
      p.visible = visible != 0;
      p.line = lineNo;
      index[p.name] = parts.size();
      parts.push_back(p);
    } else {
      throw std::invalid_argument("unknown key '" + key + "'" + where);
    }
  }
  auto get = [&](const std::string &n, const CardParticle &from) -> const CardParticle & {
    auto it = index.find(n);
    if (it == index.end())
      throw std::invalid_argument("particle '" + from.name + "' (line " + std::to_string(from.line) + ") decays to undefined '" + n + "'");
    return parts[it->second];
  };

  YAML::Emitter out;
  out.SetDoublePrecision(17);
  out << YAML::Comment("Converted from a legacy process card by kg-convert-card");
  out << YAML::BeginMap;
  out << YAML::Key << "name" << YAML::Value << "legacy_card";
  out << YAML::Key << "production" << YAML::Value << YAML::Flow << std::vector<std::string>{prodA, prodB};
  out << YAML::Key << "initial_state" << YAML::Value << YAML::BeginMap;
  out << YAML::Key << "sqrt_s" << YAML::Value << YAML::Flow << YAML::BeginMap << YAML::Key << "model" << YAML::Value
      << "exponential_pdf" << YAML::Key << "ebeam" << YAML::Value << 6500.0 << YAML::Key << "tau" << YAML::Value << 0.05
      << YAML::Key << "tau_slope" << YAML::Value << 0.02 << YAML::Key << "threshold_offset" << YAML::Value << 1.0 << YAML::EndMap;
  out << YAML::Key << "rapidity" << YAML::Value << "from_pdf";
  out << YAML::Key << "recoil" << YAML::Value << YAML::Flow << YAML::BeginMap << YAML::Key << "model" << YAML::Value
      << "exponential" << YAML::Key << "min_pt" << YAML::Value << 10.0 << YAML::Key << "mean_pt" << YAML::Value << 15.0 << YAML::EndMap;
  out << YAML::EndMap;
  out << YAML::Key << "particles" << YAML::Value << YAML::BeginMap;
  for (const auto &p : parts) {
    out << YAML::Key << p.name << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "pdg" << YAML::Value << p.pid;
    out << YAML::Key << "mass" << YAML::Value;
    if (p.mmin == p.mmax) {
      out << p.mmin;
    } else {
      out << YAML::Flow << YAML::BeginMap << YAML::Key << "dist" << YAML::Value << "uniform" << YAML::Key << "min" << YAML::Value
          << p.mmin << YAML::Key << "max" << YAML::Value << p.mmax;
      if (!opt.task_level_masses) out << YAML::Key << "scope" << YAML::Value << "event";
      out << YAML::EndMap;
    }
    if (p.stable) {
      if (!p.visible) out << YAML::Key << "visible" << YAML::Value << false;
      else if (defaultObjectType(p.pid, true).empty()) out << YAML::Key << "object" << YAML::Value << "other";
    } else {
      const CardParticle &a = get(p.d1, p), &b = get(p.d2, p);
      // Legacy rule: two-body unless m(d1)+m(d2) > m(parent); decided per event in test.cc,
      // so it must not change within the mass ranges.
      const bool lo3 = a.mmin + b.mmin > p.mmin, hi3 = a.mmax + b.mmax > p.mmax;
      const bool loM = a.mmin + b.mmin > p.mmax, hiM = a.mmax + b.mmax > p.mmin;
      if (lo3 != hi3 || lo3 != loM || lo3 != hiM)
        throw std::invalid_argument("particle '" + p.name + "': two-/three-body choice depends on where masses fall within their ranges");
      out << YAML::Key << "decays" << YAML::Value << YAML::BeginSeq;
      if (!lo3) {
        out << YAML::Flow << YAML::BeginMap << YAML::Key << "daughters" << YAML::Value << YAML::Flow
            << std::vector<std::string>{a.name, b.name} << YAML::EndMap;
      } else {
        if (a.mmin != a.mmax || b.mmin != b.mmax) {
          if ((a.mmin > b.mmin) != (a.mmax > b.mmax))
            throw std::invalid_argument("particle '" + p.name + "': legacy three-body ordering depends on sampled masses");
        }
        const bool order12 = a.mmin > b.mmin;
        const CardParticle &heavy = order12 ? a : b;
        if (heavy.stable)
          throw std::invalid_argument("particle '" + p.name + "': daughters too heavy for a two-body decay and '" + heavy.name + "' is stable");
        const CardParticle &h1 = get(heavy.d1, heavy), &h2 = get(heavy.d2, heavy);
        std::vector<std::string> d = order12 ? std::vector<std::string>{h2.name, b.name, h1.name}
                                             : std::vector<std::string>{a.name, h2.name, h1.name};
        double sum = 0;
        for (const auto &n : d) sum += parts[index[n]].mmax;
        if (sum > p.mmin) throw std::invalid_argument("particle '" + p.name + "': three-body decay is kinematically forbidden");
        out << YAML::Flow << YAML::BeginMap << YAML::Key << "daughters" << YAML::Value << YAML::Flow << d << YAML::Key << "model"
            << YAML::Value << "legacy_sequential" << YAML::EndMap;
      }
      out << YAML::EndSeq;
    }
    out << YAML::EndMap;
  }
  out << YAML::EndMap;
  out << YAML::EndMap;
  card.yaml = std::string(out.c_str()) + "\n";
  return card;
}

LegacyCard convertLegacyCard(const std::string &path, const LegacyOptions &opt) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open process card '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return convertLegacyCardString(ss.str(), opt);
}

PriorConfig priorFromLegacyCard(const std::string &path, const LegacyOptions &opt) {
  return parsePriorString(convertLegacyCard(path, opt).yaml);
}

} // namespace kg
