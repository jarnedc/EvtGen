//
//  Sample test program for running EvtGen
//  

#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtKine.hh"

#ifdef EVTGEN_EXTERNAL
#include "EvtGenExternal/EvtExternalGenList.hh"
#endif

#include <iostream>
#include <string>
#include <list>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TApplication.h"
#include "TROOT.h"

int main(int argc, char** argv) {

  //EvtParticle* parent(0);

  // Define the random number generator
  EvtRandomEngine* eng = 0;

#ifdef EVTGEN_CPP11
  // Use the Mersenne-Twister generator (C++11 only)
  eng = new EvtMTRandomEngine();
#else
  eng = new EvtSimpleRandomEngine();
#endif

  EvtRandom::setRandomEngine(eng);

  EvtAbsRadCorr* radCorrEngine = 0;
  std::list<EvtDecayBase*> extraModels;

#ifdef EVTGEN_EXTERNAL
  bool convertPythiaCodes(false);
  bool useEvtGenRandom(true);
  EvtExternalGenList genList(convertPythiaCodes, "", "gamma", useEvtGenRandom);
  radCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();
#endif

  //Initialize the generator - read in the decay table and particle properties
  EvtGen myGenerator("../DECAY.DEC","../evt.pdl", eng,
  		     radCorrEngine, &extraModels);

/*
  //If I wanted a user decay file, I would read it in now.
  //myGenerator.readUDecay("../user.dec");

  static EvtId UPS4 = EvtPDL::getId(std::string("Upsilon(4S)"));
*/
  int nEvents(100);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
  static EvtId Xi1820=EvtPDL::getId(std::string("Xi(1820)"));

   int count;
  // myGenerator.readUDecay("exampleFiles/GENERIC.DEC");
   //myGenerator.readUDecay("exampleFiles/JPSITOLL.DEC");
   TFile *file=new TFile("Xi_1820.root", "RECREATE");

   TH1F* coshel  = new TH1F("h1","cos hel",50,-1.0,1.0);

   count=1;

   do{
     std::cout << "Running for " << count  << std::endl;
     EvtVector4R p_init(EvtPDL::getMass(Xi1820),0.0,0.0,0.0);
     EvtParticle* root_part=EvtParticleFactory::particleFactory(Xi1820,
                                                                p_init);

  //   root_part->setVectorSpinDensity();
     myGenerator.generateDecay(root_part);
     EvtParticle *p = root_part;

     do{
       if (p->getId()==Xi1820) {
         EvtVector4R p4Xi1820=p->getP4Lab();
         EvtVector4R p4Daug0=p->getDaug(0)->getP4Lab();
         EvtVector4R p4Daug1=p->getDaug(1)->getP4Lab();
         double dcostheta=EvtDecayAngle(p_init,p4Xi1820,p4Daug0);
         coshel->Fill(dcostheta);
       }
       p=p->nextIter(root_part);
     }while(p!=0);

     root_part->deleteTree();
   }while (count++<nEvents);
   file->Write(); file->Close();

   EvtGenReport(EVTGEN_INFO,"EvtGen") << "SUCCESS\n";
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // Loop to create nEvents, starting from an Upsilon(4S)
/*  int i;
  for (i = 0; i < nEvents; i++) {

    std::cout<<"Event number "<<i<<std::endl;

    // Set up the parent particle
    EvtVector4R pInit(EvtPDL::getMass(UPS4), 0.0, 0.0, 0.0);
    parent = EvtParticleFactory::particleFactory(UPS4, pInit);
    parent->setVectorSpinDensity();      

    // Generate the event
    myGenerator.generateDecay(parent);    
    
    // Write out the results
    EvtHepMCEvent theEvent;
    theEvent.constructEvent(parent);
    HepMC::GenEvent* genEvent = theEvent.getEvent();
    genEvent->print(std::cout);

    parent->deleteTree();

  }

  delete eng;*/
  return 0;

}
