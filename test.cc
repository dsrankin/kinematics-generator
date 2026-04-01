#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "myGenerator.h"

using namespace std;

struct ProcessConfig {
  std::map<char, bool> pstable;
  std::map<char, bool> pvisible;
  std::map<char, std::pair<char,char>> pdecay;
  std::map<char, double> pmassMin;
  std::map<char, double> pmassMax;
  std::map<char, int> pid;
  int nevents = 1000000;
  std::string outfile = "test.root";
  char productionA = 'A';
  char productionB = 'B';
};

std::map<char, double> sampleMasses(const ProcessConfig &cfg){
  std::map<char, double> sampledMass;
  for (auto const& x : cfg.pmassMin) {
    const char p = x.first;
    const double minMass = x.second;
    const double maxMass = cfg.pmassMax.at(p);
    if (maxMass < minMass) {
      throw std::runtime_error("Invalid mass range for particle '" + std::string(1, p) + "'");
    }
    if (maxMass == minMass) {
      sampledMass[p] = minMass;
    } else {
      sampledMass[p] = gRandom->Uniform(minMass, maxMass);
    }
  }
  return sampledMass;
}

ProcessConfig loadProcessCard(const std::string &cardFile){
  ProcessConfig cfg;
  std::ifstream fin(cardFile);
  if (!fin.is_open()) {
    throw std::runtime_error("Cannot open process card: " + cardFile);
  }

  std::string line;
  int lineNo = 0;
  while (std::getline(fin, line)) {
    ++lineNo;
    std::string content = line.substr(0, line.find('#'));
    if (content.find_first_not_of(" \t\r\n") == std::string::npos) continue;

    std::istringstream iss(content);
    std::string key;
    iss >> key;
    if (key == "nevents") {
      iss >> cfg.nevents;
    } else if (key == "outfile") {
      iss >> cfg.outfile;
    } else if (key == "production") {
      std::string pa, pb;
      iss >> pa >> pb;
      if (pa.size()!=1 || pb.size()!=1) {
        throw std::runtime_error("Invalid production definition at line " + std::to_string(lineNo));
      }
      cfg.productionA = pa[0];
      cfg.productionB = pb[0];
    } else if (key == "particle") {
      std::string name;
      double massMin;
      double massMax;
      int pdgId;
      int stable;
      int visible;
      std::string d1;
      std::string d2;
      iss >> name >> massMin >> massMax >> pdgId >> stable >> visible >> d1 >> d2;
      if (name.size()!=1 || d1.size()!=1 || d2.size()!=1) {
        throw std::runtime_error("Invalid particle symbol at line " + std::to_string(lineNo));
      }
      char p = name[0];
      if (massMax < massMin) {
        throw std::runtime_error("particle '" + name + "' has mass_max < mass_min at line " + std::to_string(lineNo));
      }
      cfg.pmassMin[p] = massMin;
      cfg.pmassMax[p] = massMax;
      cfg.pid[p] = pdgId;
      cfg.pstable[p] = stable;
      cfg.pvisible[p] = visible;
      cfg.pdecay[p] = std::make_pair(d1[0], d2[0]);
    } else {
      throw std::runtime_error("Unknown key '" + key + "' at line " + std::to_string(lineNo));
    }
  }

  return cfg;
}

///////////////////////////////////////////////////////
//
// Generate pp->AB->...
//
///////////////////////////////////////////////////////
void test(const char* processCard = "process.card"){

  ProcessConfig cfg = loadProcessCard(processCard);

  // Output file
  TFile *fout = new TFile(cfg.outfile.c_str(),"recreate");
  TTree* tout = new TTree("test","test");
  TTree* tinfo = new TTree("info","info");

  float T1Pt, T1Eta, T1Phi, T1M;
  float T2Pt, T2Eta, T2Phi, T2M;
  float METPt, METEta, METPhi;

  std::vector<float> truthM;
  std::vector<float> truthMMin;
  std::vector<float> truthMMax;
  std::vector<int> truthId;

  std::vector<float> P1Pt, P1Eta, P1Phi;
  std::vector<float> P2Pt, P2Eta, P2Phi;
  std::vector<int> P1Id;
  std::vector<int> P2Id;

  tout->Branch("P1Pt", &P1Pt);
  tout->Branch("P1Eta", &P1Eta);
  tout->Branch("P1Phi", &P1Phi);
  tout->Branch("P1Id", &P1Id);
  tout->Branch("P2Pt", &P2Pt);
  tout->Branch("P2Eta", &P2Eta);
  tout->Branch("P2Phi", &P2Phi);
  tout->Branch("P2Id", &P2Id);

  tout->Branch("T1Pt", &T1Pt, "T1Pt/F");
  tout->Branch("T1Eta", &T1Eta, "T1Eta/F");
  tout->Branch("T1Phi", &T1Phi, "T1Phi/F");
  tout->Branch("T1M", &T1M, "T1M/F");
  tout->Branch("T2Pt", &T2Pt, "T2Pt/F");
  tout->Branch("T2Eta", &T2Eta, "T2Eta/F");
  tout->Branch("T2Phi", &T2Phi, "T2Phi/F");
  tout->Branch("T2M", &T2M, "T2M/F");
  tout->Branch("METPt", &METPt, "METPt/F");
  tout->Branch("METEta", &METEta, "METEta/F");
  tout->Branch("METPhi", &METPhi, "METPhi/F");

  tinfo->Branch("truthM", &truthM);
  tinfo->Branch("truthMMin", &truthMMin);
  tinfo->Branch("truthMMax", &truthMMax);
  tinfo->Branch("truthId", &truthId);

  const int seed=0;
  gRandom = new TRandom3(seed);

  for (auto const& x : cfg.pmassMin) {
    truthM.push_back(x.second);
    truthMMin.push_back(x.second);
    truthMMax.push_back(cfg.pmassMax[x.first]);
    truthId.push_back(cfg.pid[x.first]);
  }
  tinfo->Fill();

  for(int kk=0; kk<cfg.nevents; ++kk){

    if (kk%10000==0) {
      std::cout<<"Completed "<<kk<<"/"<<cfg.nevents<<std::endl;
    }

    // Reset
    P1Pt.clear();
    P1Eta.clear();
    P1Phi.clear();
    P1Id.clear();
    P2Pt.clear();
    P2Eta.clear();
    P2Phi.clear();
    P2Id.clear();

    T1Pt = 0.;
    T1Eta = 0.;
    T1Phi = 0.;
    T1M = 0.;
    T2Pt = 0.;
    T2Eta = 0.;
    T2Phi = 0.;
    T2M = 0.;
    METPt = 0.;
    METEta = 0.;
    METPhi = 0.;

    const std::map<char, double> eventMass = sampleMasses(cfg);

    // Generate event
    tlv pCM = myGenerator::getRandCM_softRad(eventMass.at(cfg.productionA)+eventMass.at(cfg.productionB)+1);
    tlv pA(0,0,0,0), pB(0,0,0,0);
    myGenerator::getRand2BodyDecay(pA,pB, pCM, eventMass.at(cfg.productionA), eventMass.at(cfg.productionB));
    T1Pt = pA.Pt();
    T1Eta = pA.Eta();
    T1Phi = pA.Phi();
    T1M = eventMass.at(cfg.productionA);
    T2Pt = pB.Pt();
    T2Eta = pB.Eta();
    T2Phi = pB.Phi();
    T2M = eventMass.at(cfg.productionB);

    tlv met(0,0,0,0);

    std::vector<std::pair<tlv,char>> intA;
    std::vector<std::pair<tlv,char>> stableA;
    std::vector<std::pair<tlv,char>> intB;
    std::vector<std::pair<tlv,char>> stableB;
    if (cfg.pstable[cfg.productionA]) {
      stableA.emplace_back(pA,cfg.productionA);
      if (!cfg.pvisible[cfg.productionA]) {
        met = met + pA;
      }
    } else {
      intA.emplace_back(pA, cfg.productionA);
    }
    if (cfg.pstable[cfg.productionB]) {
      stableB.emplace_back(pB, cfg.productionB);
      if (!cfg.pvisible[cfg.productionB]) {
        met = met + pB;
      }
    } else {
      intB.emplace_back(pB, cfg.productionB);
    }

    while (intA.size()>0) {
      for(int ip = intA.size()-1; ip > -1; ip--) {
        tlv p1(0,0,0,0), p2(0,0,0,0), p3(0,0,0,0);
        char c1 = cfg.pdecay[intA[ip].second].first;
        char c2 = cfg.pdecay[intA[ip].second].second;
        char c3;
        bool is3 = eventMass.at(c1)+eventMass.at(c2)>eventMass.at(intA[ip].second);
        double mass1 = eventMass.at(c1);
        double mass2 = eventMass.at(c2);
        double mass3 = 0;
        bool order12 = mass1>mass2;
	if (is3) {
          if (order12) {
            c3 = cfg.pdecay[c1].first;
            c1 = cfg.pdecay[c1].second;
          } else {
            c3 = cfg.pdecay[c2].first;
            c2 = cfg.pdecay[c2].second;
          }
          mass1 = eventMass.at(c1);
          mass2 = eventMass.at(c2);
          mass3 = eventMass.at(c3);
          myGenerator::getRand3BodyDecay(p1, p2, p3, intA[ip].first, mass1, mass2, mass3);
        } else {
          myGenerator::getRand2BodyDecay(p1, p2, intA[ip].first, mass1, mass2);
        }
        if (cfg.pstable[c1]) {
          stableA.emplace_back(p1, c1);
          if (!cfg.pvisible[c1]) {
            met = met + p1;
          }
        } else {
          intA.emplace_back(p1, c1);
        }
        if (cfg.pstable[c2]) {
          stableA.emplace_back(p2, c2);
          if (!cfg.pvisible[c2]) {
            met = met + p2;
          }
        } else {
          intA.emplace_back(p2, c2);
        }
        if (is3) {
          if (cfg.pstable[c3]) {
            stableA.emplace_back(p3, c3);
            if (!cfg.pvisible[c3]) {
              met = met + p3;
            }
          } else {
            intA.emplace_back(p3, c3);
          }
        }
        intA.erase(intA.begin()+ip);
      }
    }

    while (intB.size()>0) {
      for(int ip = intB.size()-1; ip > -1; ip--) {
        tlv p1(0,0,0,0), p2(0,0,0,0), p3(0,0,0,0);
        char c1 = cfg.pdecay[intB[ip].second].first;
        char c2 = cfg.pdecay[intB[ip].second].second;
        char c3;
        bool is3 = eventMass.at(c1)+eventMass.at(c2)>eventMass.at(intB[ip].second);
        double mass1 = eventMass.at(c1);
        double mass2 = eventMass.at(c2);
        double mass3 = 0;
        bool order12 = mass1>mass2;
	if (is3) {
          if (order12) {
            c3 = cfg.pdecay[c1].first;
            c1 = cfg.pdecay[c1].second;
          } else {
            c3 = cfg.pdecay[c2].first;
            c2 = cfg.pdecay[c2].second;
          }
          mass1 = eventMass.at(c1);
          mass2 = eventMass.at(c2);
          mass3 = eventMass.at(c3);
          myGenerator::getRand3BodyDecay(p1, p2, p3, intB[ip].first, mass1, mass2, mass3);
        } else {
          myGenerator::getRand2BodyDecay(p1, p2, intB[ip].first, mass1, mass2);
        }
        if (cfg.pstable[c1]) {
          stableB.emplace_back(p1, c1);
          if (!cfg.pvisible[c1]) {
            met = met + p1;
          }
        } else {
          intB.emplace_back(p1, c1);
        }
        if (cfg.pstable[c2]) {
          stableB.emplace_back(p2, c2);
          if (!cfg.pvisible[c2]) {
            met = met + p2;
          }
        } else {
          intB.emplace_back(p2, c2);
        }
        if (is3) {
          if (cfg.pstable[c3]) {
            stableB.emplace_back(p3, c3);
            if (!cfg.pvisible[c3]) {
              met = met + p3;
            }
          } else {
            intB.emplace_back(p3, c3);
          }
        }
        intB.erase(intB.begin()+ip);
      }
    }
    
    METPt = met.Pt();
    METEta = met.Eta();
    METPhi = met.Phi();

    for(unsigned int ip = 0; ip < stableA.size(); ip++) {
      P1Pt.push_back(stableA[ip].first.Pt());
      P1Eta.push_back(stableA[ip].first.Eta());
      P1Phi.push_back(stableA[ip].first.Phi());
      P1Id.push_back(cfg.pid[stableA[ip].second]);
    }
    for(unsigned int ip = 0; ip < stableB.size(); ip++) {
      P2Pt.push_back(stableB[ip].first.Pt());
      P2Eta.push_back(stableB[ip].first.Eta());
      P2Phi.push_back(stableB[ip].first.Phi());
      P2Id.push_back(cfg.pid[stableB[ip].second]);
    }

    tout->Fill();
  }

  // Print out the first event in the TTree
  std::cout << "------ Print out the first event in TTree... " << std::endl;
  tout->Show(0);

  tout->SetDirectory(fout);
  fout->Write();
  fout->Close();

  return ;
}
