# ifndef MYGENERATOR_H
# define MYGENERATOR_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"

using namespace std;
typedef TLorentzVector tlv;

namespace myGenerator{
   
  //
  // Utility functions 
  //

  inline double BreitWigner(const double &mean ,const double &width ,const double &x){  
    return 1. / (pow( pow(x,2)-pow(mean,2), 2) + pow(mean,2)*pow(width,2));  
  }

  inline TVector3 PtVect(const tlv &a){ 
    TVector3 b(a[0], a[1], 0.);    return b;  
  }
  
  inline TVector3 getRandUnitVectFlat(){ 
    
    double costh = gRandom -> Uniform(-1,1);
    double sinth = sqrt(1-pow(costh,2));
    double phi = gRandom -> Uniform(0,2*M_PI);
    
    TVector3 out(sinth*cos(phi), sinth*sin(phi), costh);  
    return out;
  }

  inline TVector3 getRandUnitVectFlat_2D(){ 
    double phi = gRandom -> Uniform(0,2*M_PI);
    TVector3 out(cos(phi), sin(phi), 0.);  
    return out;
  }
  
  inline TVector3 getRandUnitVectWdecay(){ 
    double costh = 0;
    while(1){
      costh = gRandom -> Uniform(-1,1);
      double prob = 1-costh;
      if(gRandom->Uniform(0,2)<prob) break;
    }
    double sinth = sqrt(1-pow(costh,2));
    double phi = gRandom -> Uniform(0,2*M_PI);
    
    TVector3 out(sinth*cos(phi), sinth*sin(phi), costh);  
    return out;
  }
  
  //
  inline TVector3 getRandUnitVect(TString mode="flat"){ // mode-- "W": leptonic W-decay (1-cos(theta))
    if(mode=="W") return getRandUnitVectWdecay();
    return getRandUnitVectFlat();  
  }
  
  //////////////////////////////
  //
  // Center of the mass generator
  //
  //////////////////////////////

  // ----------------------
  // Generate a center of mass assuming an exponential proton PDF: f(x) = exp(-tau*x) / tau
  // The dump constant tau is obtained from a fit of sqrt(s^) from a random unskimmed MadGraph samples.
  //
  inline double getRandCM(tlv &pq1, tlv &pq2,  // generated initial parton 4-vector
			  double min_mCM=0.,   // minimum center of mass (=sum of the generated particles)
			  double Ebeam=6500.){ // in GeV

    const double tau=0.05 + 0.02*(min_mCM/1000.);   // dump constant for the exponential proton PDF 
    double mCM=0.; // := sqrt(s^)
    double x1=0, x2=0;
    while(mCM<min_mCM){
      x1 = gRandom->Exp(tau);  if(x1>1.) continue;
      x2 = gRandom->Exp(tau);  if(x2>1.) continue;
      mCM = 2*Ebeam*sqrt(x1*x2);
    }

    pq1.SetXYZM(0, 0, Ebeam*x1, 0);
    pq2.SetXYZM(0, 0, -Ebeam*x2, 0);

    // Weight to apply when one wants to cancel this exponential PDF
    return 1. / ( exp(-x1/tau)*exp(-x2/tau)/pow(tau,2) );  

  }

  // --------------------
  // 
  inline tlv getRandCM(double min_mCM=0, double Ebeam=6500){
    tlv pq1,pq2; getRandCM(pq1,pq2,min_mCM,Ebeam);
    return pq1+pq2;
  }
  
  // -------------------
  // Smear the center of mass system to mimic the effect of (soft) radiation
  inline tlv getRandCM_softRad(double min_mCM=0, double aveRadPt=15., double Ebeam=6500){

    double pRad = 10. + gRandom->Exp(aveRadPt-10.); // dumping factor tau=15GeV (arbitrary!)
    double phi = gRandom->Uniform(0,2*M_PI);
    TVector3 vCM(pRad*cos(phi),pRad*sin(phi),0);  

    tlv pCM =  getRandCM(min_mCM, Ebeam);
    TVector3 betaCM_z = pCM.BoostVector();
    tlv pCM_rest;  pCM_rest.SetXYZM(0,0,0,pCM.M());
    TVector3 betaCM_t = -vCM * (1./pCM.M());

    tlv pCM_recoiled = pCM_rest;   
    pCM_recoiled.Boost(betaCM_t);
    pCM_recoiled.Boost(betaCM_z);

    return pCM_recoiled;
  }


  // -------------------
  // Generate a mass assuming to a Breit-Wigner width
  inline double getRandMass_BW(double mean, double width){

    double x=0;
    double prob_max = BreitWigner(mean,width,mean);

    int ig=0;
    while(1){
      x = gRandom->Gaus(mean,width*100);
      if(x<0) continue;

      double prob_offset = pow((x-mean)/width,2);
      double prob = BreitWigner(mean,width,x) / prob_offset;

      if(gRandom->Uniform(0,prob_max) < prob) return x;

      if(ig>1000) {
	cout << "myGenerator::getRandMass_BW  ERROR  Do not converge. Please set valid mean/width which is now " << mean << "/" << width << ". Return 0.  " << endl;
	return 0.;
      }

    }
  }

  //////////////////////////////
  //
  // Decay generator
  //
  //////////////////////////////
 
  // -------------------
  inline tlv getRand4Mom(const double &mom, const double &mass){
    
    tlv tl;
    tl.SetVectM( mom*getRandUnitVect(), mass );
    
    return tl;
  }
  
  // -------------------
  inline double get2BodyDecayMom(const double &m0, const double &m1, const double &m2){ 
    double M0 = pow(m0,2);    double M1 = pow(m1,2);    double M2 = pow(m2,2);
    double delta = pow(M0 - M1 + M2, 2);
    return sqrt( delta/(4*M0) - M2);
  }

  // -------------------  
  // Generate a 0 -> 1 2 type decay. (For case where particle 0 is at rest.)
  inline void getRand2BodyDecay(tlv &p1, tlv &p2, const double &m0, double m1=0, double m2=0, TString mode="flat"){
    if(m0<m1+m2) {
      cout << "kineFunc::getRand2BodyDecay() ---- error! m0(" << m0 << ")<m1(" << m1 << ")+m2(" << m2 << "). exit." << endl;
      return;
    }  
    TVector3 mom = get2BodyDecayMom(m0, m1, m2)*getRandUnitVect(mode);
    p1.SetVectM(mom , m1);
    p2.SetVectM(-mom , m2);
    
    return;
  }
  
  // -------------------
  // Generate a 0 -> 1 2 type decay. (For case where particle 0 is not at rest.)
  inline void getRand2BodyDecay(tlv &p1, tlv &p2, const tlv  &p0, double m1=0, double m2=0, TString mode="flat"){
    
    getRand2BodyDecay(p1, p2, p0.M(), m1, m2, mode);        
    p1.Boost(p0.BoostVector());
    p2.Boost(p0.BoostVector());

    return;
  }
  // -------------------
  inline void getRand3BodyDecay(tlv &p1, tlv &p2, tlv &p3, const double &m0, double m1=0, double m2=0, double m3=0, double mY=0, TString mode0="flat", TString modeY="flat"){
    
    // Topology:  0 -> 1 Y, Y -> 2 3
    //  mY is treated to uniformly distributed when mY is set to 0.
    
    if(m0<m1+m2+m3) {
      cout << "kineFunc::getRand3BodyDecay  INFO  error! m0<m1+m2+m3." << endl;
      cout << "kineFunc::getRand3BodyDecay  INFO  m0: " << m0 << " m1: " << m1 << " m2: " << m2 << " m3: " << m3 << ". Exit." << endl;
      return;
    }
    if(!mY) mY = gRandom -> Uniform(m2+m3, m0-m1);
    else if(mY<m2+m3 || mY>m0-m1){
      cout << "kineFunc::getRand3BodyDecay  INFO  error! mY<m2+m3 or mY>m0-m1." << endl;
      cout << "kineFunc::getRand3BodyDecay  INFO  m0: " << m0 << " m1: " << m1 << " m2: " << m2 << " m3: " << m3 << " mY: " << mY << ". Exit." << endl;
      return;
    }
    
    TVector3 mom1 = get2BodyDecayMom(m0, m1, mY)*getRandUnitVect(mode0); // 0->1,Y
    p1.SetVectM(mom1 , m1);
    tlv pY;  pY.SetVectM(-mom1 , mY);
    
    TVector3 mom2 = get2BodyDecayMom(mY, m2, m3)*getRandUnitVect(modeY); //Y->2,3
    p2.SetVectM(mom2 , m2);   p2.Boost(pY.BoostVector());
    p3.SetVectM(-mom2 , m3);  p3.Boost(pY.BoostVector());    
    
    //cout << mY << " " << (p2+p3).M() << " " << (p1+p2+p3).M() << endl;
    
    return ;
  }
  
  // -------------------
  // Generate a 0 -> 1 2 type decay. (For case where particle 0 is not at rest.)
  inline void getRand3BodyDecay(tlv &p1, tlv &p2, tlv &p3, const tlv &p0, double m1=0, double m2=0, double m3=0, double mY=0, TString mode0="flat", TString modeY="flat"){
    
    getRand3BodyDecay(p1, p2, p3, p0.M(), m1, m2, m3, mY, mode0, modeY);
    p1.Boost(p0.BoostVector());
    p2.Boost(p0.BoostVector());
    p3.Boost(p0.BoostVector());
    
    return ;
  }

  // --------------------
  inline void smearMET(tlv &pMis){
    double px = pMis.Px()*(1+gRandom->Gaus(0,0.2));
    double py = pMis.Py()*(1+gRandom->Gaus(0,0.2));  
    pMis.SetXYZM(px,py,0,0);
  }

  //////////////////////////////
  //
  // Process generator
  //
  //////////////////////////////
 

  ////////////////////////////////////////////////////////////
  // SM W+ISR production
  inline void getKine_Wjet(tlv &pl, tlv &pv, 
			   double betaT_CM=0. // Transverse recoil of the W
			   ){

    // pp->W
    const double mW = getRandMass_BW(80.4, 2.0);     
    tlv pW;  pW.SetXYZM(0,0,0,mW);     

    // W->Wl Wv
    pl.SetXYZM(0,0,0,0); pv.SetXYZM(0,0,0,0);   
    getRand2BodyDecay(pl, pv, pW, 0, 0, "W");  
        
    TVector3 vbetaT_CM = betaT_CM * getRandUnitVectFlat_2D();
    pl.Boost(vbetaT_CM);
    pv.Boost(vbetaT_CM);
    return;
  }

  ////////////////////////////////////////////////////////////
  // SM WW production
  //   pp->WZ->lvll
  //   Symbols:  W->Wl+Wv,  Z->Zlp+Zlm
  inline void getKine_WZ(tlv &pW, tlv &pWl, tlv &pWv, 
			 tlv &pZ, tlv &pZlp, tlv &pZlm,
			 double mW=80.4, double mZ=91.2){

    if(mW<0){ mW = getRandMass_BW(80.4, 2.0); }
    if(mZ<0){ mZ = getRandMass_BW(91.2, 2.2); }

    // pp->WZ
    tlv pCM = getRandCM_softRad(mW+mZ+1);
    getRand2BodyDecay(pW,pZ, pCM, mW, mZ);
    
    // W -> Wl Wv
    pWl.SetXYZM(0,0,0,0); pWv.SetXYZM(0,0,0,0);   
    getRand2BodyDecay(pWl, pWv, pW, 0, 0, "W");  
    
    // Z -> Zlp Zlm
    pZlp.SetXYZM(0,0,0,0); pZlm.SetXYZM(0,0,0,0);   
    getRand2BodyDecay(pZlp, pZlm, pZ, 0, 0);  
    
    return;
  }

  ////////////////////////////////////////////////////////////
  // WW production
  //   pp->WW->lvlv
  //   Symbols:  W1->l1+v1,  W2->l2+v2
  inline void getKine_WW(tlv &pW1, tlv &pl1, tlv &pv1, 
			 tlv &pW2, tlv &pl2, tlv &pv2, 
			 double mW1=80.4, double mW2=80.4){

    if(mW1<0){ mW1 = getRandMass_BW(80.4, 2.0); }
    if(mW2<0){ mW2 = getRandMass_BW(80.4, 2.0); }
    
    // pp->WW
    tlv pCM = getRandCM_softRad(mW1+mW2+1.);
    getRand2BodyDecay(pW1,pW2, pCM, mW1, mW2);
    
    // W1 -> l1 v1
    pl1.SetXYZM(0,0,0,0); pv1.SetXYZM(0,0,0,0);   
    getRand2BodyDecay(pl1, pv1, pW1, 0, 0, "W");  

    // W2 -> l2 v2
    pl2.SetXYZM(0,0,0,0); pv2.SetXYZM(0,0,0,0);   
    getRand2BodyDecay(pl2, pv2, pW2, 0, 0, "W");  
    
    return;
  }

  ////////////////////////////////////////////////////////////  
  // WW+ISR production
  //   pp->WW+j->lvlv+j
  //   Symbols:  W1->l1+v1,  W2->l2+v2
  inline void getKine_WWj(tlv &pW1, tlv &pW2, 
			  tlv &pl1, tlv &pv1,
			  tlv &pl2, tlv &pv2,
			  double betaT_CM=0. // Transverse recoil of the WW system
			  ){

    const double mW = 80.4;

    // pp->WW    
    tlv pCM = myGenerator::getRandCM(2*mW+1);
    getRand2BodyDecay(pW1,pW2, pCM, mW, mW);

    // W1->l1+v1    
    pl1.SetXYZM(0,0,0,0);  pv1.SetXYZM(0,0,0,0); 
    getRand2BodyDecay(pl1, pv1, pW1, 0, 0);
    
    // W2->l2+v2    
    pl2.SetXYZM(0,0,0,0);  pv2.SetXYZM(0,0,0,0); 
    getRand2BodyDecay(pl2, pv2, pW2, 0, 0);
    
    // Boost the WW system
    TVector3 vbetaT_CM = betaT_CM * getRandUnitVectFlat_2D();
    pW1.Boost(vbetaT_CM);
    pl1.Boost(vbetaT_CM);
    pv1.Boost(vbetaT_CM);
    pW2.Boost(vbetaT_CM);
    pl2.Boost(vbetaT_CM);
    pv2.Boost(vbetaT_CM);    
    return;
  }  

  ////////////////////////////////////////////////////////////  
  // ttbar production
  //  pp->tt->blvblv
  //  Symbols:  t1->b1+W1, W1->l1+v1, t2->b2+W2, W2->l2+v2
  inline void getKine_TT(tlv &pt1, tlv &pb1, tlv &pW1, 
		  tlv &pl1, tlv &pv1,
		  tlv &pt2, tlv &pb2, tlv &pW2, 
		  tlv &pl2, tlv &pv2){

    // pp->tt
    const double mtop = 172.5;
    tlv pCM = getRandCM_softRad(2*mtop+1);    
    getRand2BodyDecay(pt1,pt2, pCM, mtop, mtop);

    // t1->b1+W1
    pb1.SetXYZM(0,0,0,0); pW1.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pb1, pW1, pt1, 5, 80.4);

    // W1->l1+v1    
    pl1.SetXYZM(0,0,0,0);  pv1.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pl1, pv1, pW1, 0, 0);
    
    // t2->b2+W2
    pb2.SetXYZM(0,0,0,0);  pW2.SetXYZM(0,0,0,0);   
    getRand2BodyDecay(pb2, pW2, pt2, 5, 80.4);
    
    // W2->l2+v2    
    pl2.SetXYZM(0,0,0,0);  pv2.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pl2, pv2, pW2, 0, 0);
    
    return;
  }  
  ////////////////////////////////////////////////////////////
  // EWKino production, direct 2-body decay into a SM boson
  //  pp->XX->BB+N1N1 where B:= SM boson (W, Z, h)
  //  Symbols:  X1->B1+N1,  X2->B2+N2
  inline void getKine_XX_BB(double mC1, double mN1, double mB1, double mB2,
			    tlv &pC1, tlv &pB1, tlv &pN1,
			    tlv &pC2, tlv &pB2, tlv &pN2){
    
    // C1,2: gaugino 1,2
    // B1,2: boson from gaugino 1,2
    // N1,2: LSP 1,2
    
    tlv pCM = getRandCM_softRad(2*mC1+1);
    getRand2BodyDecay(pC1,pC2, pCM, mC1, mC1);
    
    pB1.SetXYZM(0,0,0,0); pN1.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pB1, pN1, pC1, mB1, mN1);
    
    pB2.SetXYZM(0,0,0,0); pN2.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pB2, pN2, pC2, mB2, mN1);
    
    return;
  }
  ////////////////////////////////////////////////////////////
  // EWKino production, direct 3-body decay into LSP
  //  pp->XX->qqqq+N1N1
  //  Symbols:  X1->q11+q12+N1,  X2->q21+q22+N2
  inline void getKine_XX(double mX1, double mN1, 		
			 tlv &pX1, tlv &pq11, tlv &pq12, tlv &pN1,
			 tlv &pX2, tlv &pq21, tlv &pq22, tlv &pN2,
			 TString modeX1="flat", TString modeX2="flat"){
    
    // X1,2: EWKino 1,2
    // N1,2: LSP 1,2
    
    TString modeX1_boson = modeX1=="C" ? "W" : "flat";
    TString modeX2_boson = modeX2=="C" ? "W" : "flat";

    // pp->X1+X2
    tlv pCM = getRandCM_softRad(2*mX1+1);
    getRand2BodyDecay(pX1,pX2, pCM, mX1, mX1);
    
    // X1->q11+q12+N1
    pq12.SetXYZM(0,0,0,0); pN1.SetXYZM(0,0,0,0);
    getRand3BodyDecay(pq11, pq12, pN1, pX1, 0, 0, mN1, 0, modeX1, modeX1_boson);
    
    // X2->q21+q22+N2
    pq22.SetXYZM(0,0,0,0); pN2.SetXYZM(0,0,0,0);
    getRand3BodyDecay(pq21, pq22, pN2, pX2, 0, 0, mN1, 0, modeX2, modeX2_boson);
    
    return;
  }

  ////////////////////////////////////////////////////////////
  // Gluino production, one-step decay via EWKino
  //  pp->GG->XX+qqqq->BB+N1N1+qqqq->lvlv+N1N1+qqqq where B:= SM boson (W, Z, h)
  //  Symbols:  G1->X1+q11+q12, X1->B1+N1, B1->l1+v1+N1, G2->X2+q21+q22, X2->B2+N2, B2->l2+v2+N2
  inline void getKine_GG_XX(double mg, double mX1, double dm, 
		     tlv &pg1, tlv &pq11, tlv &pq12, tlv &pX1, 
		     tlv &pg2, tlv &pq21, tlv &pq22, tlv &pX2, 
		     tlv &pl1, tlv &pv1, tlv &pN1, 
		     tlv &pl2, tlv &pv2, tlv &pN2){
    
    const double mN1 = mX1-dm;

    // pp->g1+g2
    tlv pCM = getRandCM_softRad(2*mg+1);
    getRand2BodyDecay(pg1,pg2, pCM, mg, mg);
    
    // g1->X1 q11 q12
    pq11.SetXYZM(0,0,0,0); pq12.SetXYZM(0,0,0,0);   pX1.SetXYZM(0,0,0,0);
    getRand3BodyDecay(pq11, pq12, pX1, pg1, 0, 0, mX1);  
    // X1->l1 v1 N1
    pv1.SetXYZM(0,0,0,0); pN1.SetXYZM(0,0,0,0);
    getRand3BodyDecay(pl1, pv1, pN1, pX1, 0, 0, mN1);
    // g2->X2 q21 q22
    pq21.SetXYZM(0,0,0,0); pq22.SetXYZM(0,0,0,0); pX2.SetXYZM(0,0,0,0);
    getRand3BodyDecay(pq21, pq22, pX2, pg2, 0, 0, mX1);
    // X2->l2 v2 N2
    pv2.SetXYZM(0,0,0,0); pN2.SetXYZM(0,0,0,0);
    getRand3BodyDecay(pl2, pv2, pN2, pX2, 0, 0, mN1);
    
    return;
  }
  ////////////////////////////////////////////////////////////
  // Stop production, direct decay into top+LSP
  //  pp->sTsT->tt+N1N1->bWbW+N1N1->blvblv+N1N1
  //  Symbols:  st1->t1+N1, t1->b1+W1, W1->l1+v1, st2->t2+N2, t2->b2+W2, W2->l2+v2
  inline void getKine_sTsT(double mst, double mN, 
		    tlv &pst1, tlv &pt1, tlv &pN1, 
		    tlv &pst2, tlv &pt2, tlv &pN2, 
		    tlv &pb1, tlv &pW1, 
		    tlv &pb2, tlv &pW2, 
		    tlv &pl1, tlv &pv1, 
		    tlv &pl2, tlv &pv2){
    
    const double mtop = 172.5;

    // pp->st1+st2
    tlv pCM = getRandCM_softRad(2*mst+1);
    getRand2BodyDecay(pst1,pst2, pCM, mst, mst);

    // st1->t1+N1
    pt1.SetXYZM(0,0,0,0);  pN1.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pt1, pN1, pst1, mtop, mN);
    
    // t1->b1+W1
    pb1.SetXYZM(0,0,0,0);  pW1.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pb1, pW1, pt1, 0, 80.4);
    
    // W1->l1+v1
    pl1.SetXYZM(0,0,0,0);  pv1.SetXYZM(0,0,0,0); 
    getRand2BodyDecay(pl1, pv1, pW1, 0, 0);
            
    // st2->t2+N2
    pt2.SetXYZM(0,0,0,0);  pN2.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pt2, pN2, pst2, mtop, mN);
    
    // t2->b2+W2
    pb2.SetXYZM(0,0,0,0);  pW2.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pb2, pW2, pt2, 0, 80.4);
    
    // W2->l2+v2
    pl2.SetXYZM(0,0,0,0);  pv2.SetXYZM(0,0,0,0); 
    getRand2BodyDecay(pl2, pv2, pW2, 0, 0);


    return;
  }
  ////////////////////////////////////////////////////////////
  // Stop production, decay into b+chargino
  //  pp->sTsT->tt+N1N1->bWbW+N1N1->blvblv+N1N1
  //  Symbols:  st1->b1+C1, C1->l1+v1+N1, st2->b2+C2, C2->l2+v2+N2 
  inline void getKine_sTsT_XX(double mst, double mX1, double dm, 
		       tlv &pst1, tlv &pb1, tlv &pX1, 
		       tlv &pst2, tlv &pb2, tlv &pX2, 
		       tlv &pl1, tlv &pv1, tlv &pN1, 
		       tlv &pl2, tlv &pv2, tlv &pN2){
    
    const double mN1 = mX1-dm;
    tlv pCM = getRandCM_softRad(2*mst+1);
    getRand2BodyDecay(pst1,pst2, pCM, mst, mst);
    
    pb1.SetXYZM(0,0,0,0); pX1.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pb1, pX1, pst1, 0, mX1);
    pv1.SetXYZM(0,0,0,0); pN1.SetXYZM(0,0,0,0);
    getRand3BodyDecay(pl1, pv1, pN1, pX1, 0, 0, mN1);
    
    pb2.SetXYZM(0,0,0,0); pX2.SetXYZM(0,0,0,0);
    getRand2BodyDecay(pb2, pX2, pst2, 0, mX1);
    pv2.SetXYZM(0,0,0,0); pN2.SetXYZM(0,0,0,0);
    getRand3BodyDecay(pl2, pv2, pN2, pX2, 0, 0, mN1);
    
    return;
  } 
  ///////////////////////////////////////////////////////////////////////


}
      




#endif

