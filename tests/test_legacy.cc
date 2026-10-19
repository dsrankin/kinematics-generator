// Compatibility adapter for legacy process cards.
#include <set>

#include "helpers.h"
#include "test_framework.h"

using namespace kg;

TEST(legacy_card_conversion) {
  const LegacyCard card = convertLegacyCard(std::string(KG_SOURCE_DIR) + "/process.card");
  CHECK(card.nevents == 1000000);
  CHECK(card.outfile == "test_TbqqTblv.root");
  const PriorConfig p = parsePriorString(card.yaml);
  // A (5.5 GeV) -> W b is below threshold, so the legacy rule turns it into
  // A -> q b q with W replaced by its daughters, using the sequential model.
  const ParticleDef &A = p.particle("A");
  CHECK(A.modes.size() == 1 && A.modes[0].daughters == (std::vector<std::string>{"q", "b", "q"}));
  CHECK(A.modes[0].model == ThreeBodyModel::LegacySequential);
  const ParticleDef &B = p.particle("B");
  CHECK(B.modes[0].daughters == (std::vector<std::string>{"w", "b"}));
  CHECK(p.particle("v").object == "invisible" && p.particle("N").object == "invisible");
  CHECK(p.particle("b").object == "bjet" && p.particle("l").object == "electron");
  CHECK(p.params.isEvent(p.params.at("B.mass")));
  CHECK(!p.params.isEvent(p.params.at("A.mass")));
}

TEST(legacy_card_generation) {
  const PriorConfig p = priorFromLegacyCard(std::string(KG_SOURCE_DIR) + "/process.card");
  const TaskConfig t = sampleTask(p, 0, 0);
  std::set<double> massesB;
  for (int i = 0; i < 500; ++i) {
    const TruthEvent ev = generateTruthEvent(t, i);
    const auto prod = ev.producedParticles();
    CHECK(prod.size() == 2);
    const TruthParticle &a = ev.particles[prod[0]], &b = ev.particles[prod[1]];
    CHECK(a.name == "A" && b.name == "B");
    CHECK(b.mass >= 165 && b.mass < 180);
    massesB.insert(b.mass);
    CHECK(a.children.size() == 3);
    // B -> w b, w -> l v
    std::multiset<std::string> fs;
    for (int f : ev.finalState())
      if (ev.particles[f].branch == 1) fs.insert(ev.particles[f].name);
    CHECK(fs == (std::multiset<std::string>{"b", "l", "v"}));
    // Legacy threshold: sqrt(s) >= mA + mB + 1.
    CHECK(ev.sqrt_s_hat >= a.mass + b.mass + 1);
    const P4 in = ev.incomingSum(), out = ev.finalSum();
    CHECK_CLOSE(in.e, out.e, 1e-8 * in.e);
  }
  CHECK(massesB.size() == 500); // re-sampled per event, as in test.cc
  // Task-level masses: one value per task.
  LegacyOptions opt;
  opt.task_level_masses = true;
  const TaskConfig u = sampleTask(priorFromLegacyCard(std::string(KG_SOURCE_DIR) + "/process.card", opt), 0, 0);
  std::set<double> mb;
  for (int i = 0; i < 20; ++i) mb.insert(generateTruthEvent(u, i).particles[generateTruthEvent(u, i).producedParticles()[1]].mass);
  CHECK(mb.size() == 1);
}

TEST(legacy_card_errors) {
  CHECK_THROWS_AS(convertLegacyCardString("particle A 10 10 1 0 1 B q\nparticle q 0 0 1 1 1 x x\n"), std::invalid_argument,
                  "undefined 'B'");
  CHECK_THROWS_AS(convertLegacyCardString("frobnicate 3\n"), std::invalid_argument, "unknown key");
  // Two-/three-body choice must not depend on where a mass falls in its range.
  CHECK_THROWS_AS(convertLegacyCardString("production A A\nparticle A 50 150 1 0 1 W q\nparticle W 80.4 80.4 24 0 0 q q\n"
                                          "particle q 0 0 1 1 1 x x\n"),
                  std::invalid_argument, "depends on where masses fall");
}
