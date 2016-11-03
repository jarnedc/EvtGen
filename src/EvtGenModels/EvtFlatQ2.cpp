//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtFlatQ2.cc
//
// Description: B->Xu l nu with flat q2 distribution
//
// Modification history:
//
//    David Cote, U. de Montreal, 11/02/2003    Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenModels/EvtFlatQ2.hh"

#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <string>
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
using std::fstream;
#include <iomanip>

double lambda(double q, double m_mu)
{

  double L(1.0);
  if (fabs(q) > 0.0) {
      L -= 4.0*m_mu*m_mu/(q*q);
  }

  return L;

}

EvtFlatQ2::~EvtFlatQ2() {}

std::string EvtFlatQ2::getName(){

    return "FLATQ2";

}

EvtDecayBase* EvtFlatQ2::clone(){

  return new EvtFlatQ2;

}


void EvtFlatQ2::initProbMax(){

  setProbMax(100);

}


void EvtFlatQ2::init(){

  // check that there are 3 daughters
  checkNDaug(3);

  // We expect B -> X lepton lepton events
  checkSpinParent(EvtSpinType::SCALAR);

  EvtSpinType::spintype d1type = EvtPDL::getSpinType(getDaug(1));
  EvtSpinType::spintype d2type = EvtPDL::getSpinType(getDaug(2));

  if (!(d1type == EvtSpinType::DIRAC || d1type == EvtSpinType::NEUTRINO)) {
      EvtGenReport(EVTGEN_ERROR,"EvtGen") << "EvtFlatQ2 expects 2nd daughter to "
					  << "be a lepton" <<std::endl;
      EvtGenReport(EVTGEN_ERROR,"EvtGen") << "Will terminate execution!"<<std::endl;
      ::abort();
  }

  if (!(d2type == EvtSpinType::DIRAC || d2type == EvtSpinType::NEUTRINO)) {
      EvtGenReport(EVTGEN_ERROR,"EvtGen") << "EvtFlatQ2 expects 3rd daughter to "
					  << "be a lepton" <<std::endl;
      EvtGenReport(EVTGEN_ERROR,"EvtGen") << "Will terminate execution!"<<std::endl;
      ::abort();
  }

  // Specify if we want to use the phase space factor
  _usePhsp = false;
  if (getNArg() > 0) {
      if (getArg(0) != 0) {_usePhsp = true;}
  }

  EvtGenReport(EVTGEN_INFO,"EvtGen") <<"EvtFlatQ2 usePhsp = "<<int(_usePhsp)<<std::endl;

}


void EvtFlatQ2::decay( EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtVector4R p4Xu = p->getDaug(0)->getP4();

  EvtVector4R p4ell1 = p->getDaug(1)->getP4();
  EvtVector4R p4ell2 = p->getDaug(2)->getP4();

  double pXu_x2 = p4Xu.get(1)*p4Xu.get(1);
  double pXu_y2 = p4Xu.get(2)*p4Xu.get(2);
  double pXu_z2 = p4Xu.get(3)*p4Xu.get(3);
  double pXu = sqrt(pXu_x2+pXu_y2+pXu_z2);
  double prob(0.0);
  if (fabs(pXu) > 0.0) {prob = 1/pXu;}

  // Include the phase space factor if requested
  if (_usePhsp) {

    double Lambda = lambda((p4ell1+p4ell2).mass(), p4ell1.mass());
    if (Lambda > 0.0) {prob=prob/sqrt(Lambda);}
	  
  }

  if (pXu > 0.01) {setProb(prob);}

  return;

}


