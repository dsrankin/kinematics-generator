#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include <map>
#include <vector>
#include <iostream>

#include "myGenerator.h"

using namespace std;

std::map<char, bool> pstable = {
    { 'A', false },
    { 'B', false },
    { 'C', false },
    { 'D', false },
    { 'E', false },
    { 'W', false },
    { 'w', false },
    { 'Z', false },
    { 'z', false },
    { 'h', false },
    { 'T', false },
    { 't', false },
    { 'b', true },
    { 'q', true },
    { 'l', true },
    { 'v', true },
    { 'N', true }
};

std::map<char, bool> pvisible = {
    { 'A', false },
    { 'B', false },
    { 'C', false },
    { 'D', false },
    { 'E', false },
    { 'W', false },
    { 'w', false },
    { 'Z', false },
    { 'z', false },
    { 'h', false },
    { 'T', false },
    { 't', false },
    { 'b', true },
    { 'q', true },
    { 'l', true },
    { 'v', false },
    { 'N', false }
};

std::map<char, std::pair<char,char>> pdecay = {
    { 'A', {'W', 'b'} },
    { 'B', {'w', 'b'} },
    { 'C', {'D', 'q'} },
    { 'D', {'E', 'q'} },
    { 'E', {'N', 'q'} },
    { 'W', {'q', 'q' } },
    { 'w', {'l', 'v' } },
    { 'Z', {'q', 'q'} },
    { 'z', {'l', 'l'} },
    { 'h', {'b', 'b'} },
    { 'T', {'W', 'b'} },
    { 't', {'w', 'b'} }
};

std::map<char, double> pmass = {
    //{ 'A', 1000. },
    //{ 'B', 1000. },
    { 'A', 5.5 },
    { 'B', 172.5 },
    { 'C', 400. },
    { 'D', 350. },
    { 'E', 300. },
    { 'W', 80.4 },
    { 'w', 80.4 },
    { 'Z', 91.2 },
    { 'z', 91.2 },
    { 'h', 125. },
    { 'T', 172.5 },
    { 't', 172.5 },
    { 'b', 4.2 },
    { 'q', 0. },
    { 'l', 0. },
    { 'v', 0. },
    { 'N', 200. }
};

std::map<char, int> pid = {
    { 'A', 9999 },
    { 'B', 9998 },
    { 'C', 9997 },
    { 'D', 9996 },
    { 'E', 9995 },
    { 'W', 23 },
    { 'w', 23 },
    { 'Z', 24 },
    { 'z', 24 },
    { 'h', 25 },
    { 'T', 6 },
    { 't', 6 },
    { 'b', 5 },
    { 'q', 1 },
    { 'l', 11 },
    { 'v', 12 },
    { 'N', 5000 }
};
int nevents = 1000000;

///////////////////////////////////////////////////////
//
// Generate pp->AB->...
//
///////////////////////////////////////////////////////
void test(){

  // Output file
  //TFile *fout = new TFile("test_WqqWqqNN.root","recreate");
  TFile *fout = new TFile("test_TbqqTblv.root","recreate");
  //TFile *fout = new TFile("test_AqqqqNBqqqqN.root","recreate");
  TTree* tout = new TTree("test","test");
  TTree* tinfo = new TTree("info","info");

  float T1Pt, T1Eta, T1Phi, T1M;
  float T2Pt, T2Eta, T2Phi, T2M;
  float METPt, METEta, METPhi;

  std::vector<float> truthM;
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
  tinfo->Branch("truthId", &truthId);

  const int seed=0;
  TRandom3 *rand3 = new TRandom3(seed);

  for (auto const& x : pmass) {
    truthM.push_back(x.second);
    truthId.push_back(pid[x.first]);
  }
  tinfo->Fill();

  for(int kk=0; kk<nevents; ++kk){

    if (kk%10000==0) {
      std::cout<<"Completed "<<kk<<"/"<<nevents<<std::endl;
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

    // Generate event
    tlv pCM = myGenerator::getRandCM_softRad(pmass['A']+pmass['B']+1);
    tlv pA(0,0,0,0), pB(0,0,0,0);
    myGenerator::getRand2BodyDecay(pA,pB, pCM, pmass['A'], pmass['B']);
    T1Pt = pA.Pt();
    T1Eta = pA.Eta();
    T1Phi = pA.Phi();
    T1M = pmass['A'];
    T2Pt = pB.Pt();
    T2Eta = pB.Eta();
    T2Phi = pB.Phi();
    T2M = pmass['B'];

    tlv met(0,0,0,0);

    std::vector<std::pair<tlv,char>> intA;
    std::vector<std::pair<tlv,char>> stableA;
    std::vector<std::pair<tlv,char>> intB;
    std::vector<std::pair<tlv,char>> stableB;
    if (pstable['A']) {
      stableA.emplace_back(pA,'A');
      if (!pvisible['A']) {
        met = met + pA;
      }
    } else {
      intA.emplace_back(pA, 'A');
    }
    if (pstable['B']) {
      stableB.emplace_back(pB, 'B');
      if (!pvisible['B']) {
        met = met + pB;
      }
    } else {
      intB.emplace_back(pB, 'B');
    }

    while (intA.size()>0) {
      for(int ip = intA.size()-1; ip > -1; ip--) {
        tlv p1(0,0,0,0), p2(0,0,0,0), p3(0,0,0,0);
        char c1 = pdecay[intA[ip].second].first;
        char c2 = pdecay[intA[ip].second].second;
        char c3;
        bool is3 = pmass[c1]+pmass[c2]>pmass[intA[ip].second];
        double mass1 = pmass[c1];
        double mass2 = pmass[c2];
        double mass3 = 0;
        bool order12 = mass1>mass2;
	if (is3) {
          if (order12) {
            c3 = pdecay[c1].first;
            c1 = pdecay[c1].second;
          } else {
            c3 = pdecay[c2].first;
            c2 = pdecay[c2].second;
          }
          mass1 = pmass[c1];
          mass2 = pmass[c2];
          mass3 = pmass[c3];
          myGenerator::getRand3BodyDecay(p1, p2, p3, intA[ip].first, mass1, mass2, mass3);
        } else {
          myGenerator::getRand2BodyDecay(p1, p2, intA[ip].first, mass1, mass2);
        }
        if (pstable[c1]) {
          stableA.emplace_back(p1, c1);
          if (!pvisible[c1]) {
            met = met + p1;
          }
        } else {
          intA.emplace_back(p1, c1);
        }
        if (pstable[c2]) {
          stableA.emplace_back(p2, c2);
          if (!pvisible[c2]) {
            met = met + p2;
          }
        } else {
          intA.emplace_back(p2, c2);
        }
        if (is3) {
          if (pstable[c3]) {
            stableA.emplace_back(p3, c3);
            if (!pvisible[c3]) {
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
        char c1 = pdecay[intB[ip].second].first;
        char c2 = pdecay[intB[ip].second].second;
        char c3;
        bool is3 = pmass[c1]+pmass[c2]>pmass[intB[ip].second];
        double mass1 = pmass[c1];
        double mass2 = pmass[c2];
        double mass3 = 0;
        bool order12 = mass1>mass2;
	if (is3) {
          if (order12) {
            c3 = pdecay[c1].first;
            c1 = pdecay[c1].second;
          } else {
            c3 = pdecay[c2].first;
            c2 = pdecay[c2].second;
          }
          mass1 = pmass[c1];
          mass2 = pmass[c2];
          mass3 = pmass[c3];
          myGenerator::getRand3BodyDecay(p1, p2, p3, intB[ip].first, mass1, mass2, mass3);
        } else {
          myGenerator::getRand2BodyDecay(p1, p2, intB[ip].first, mass1, mass2);
        }
        if (pstable[c1]) {
          stableB.emplace_back(p1, c1);
          if (!pvisible[c1]) {
            met = met + p1;
          }
        } else {
          intB.emplace_back(p1, c1);
        }
        if (pstable[c2]) {
          stableB.emplace_back(p2, c2);
          if (!pvisible[c2]) {
            met = met + p2;
          }
        } else {
          intB.emplace_back(p2, c2);
        }
        if (is3) {
          if (pstable[c3]) {
            stableB.emplace_back(p3, c3);
            if (!pvisible[c3]) {
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
      P1Id.push_back(pid[stableA[ip].second]);
    }
    for(unsigned int ip = 0; ip < stableB.size(); ip++) {
      P2Pt.push_back(stableB[ip].first.Pt());
      P2Eta.push_back(stableB[ip].first.Eta());
      P2Phi.push_back(stableB[ip].first.Phi());
      P2Id.push_back(pid[stableB[ip].second]);
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

