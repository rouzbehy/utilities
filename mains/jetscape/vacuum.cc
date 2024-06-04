/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
// ------------------------------------------------------------
// JetScape Framework jet in hydro from file Test Program
// (use either shared library (need to add paths; see setup.csh)
// (or create static library and link in)
// -------------------------------------------------------------

#include <iostream>
#include <time.h>
#include <chrono>
#include <thread>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

// User modules derived from jetscape framework clasess
// to be used to run Jetscape ...
#include "Matter.h"
#include "LBT.h"
#include "Brick.h"
#include "HydroFromFile.h"
#include "PythiaGun.h"
#include "PartonPrinter.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColorlessHadronization.h"

#ifdef USE_HDF5
#include "InitialFromFile.h"
#endif
// using namespace std;
// Add initial state module for test
#include "TrentoInitial.h"

#include <chrono>
#include <thread>

using namespace Jetscape;

// Forward declaration
void Show();

// -------------------------------------

int main(int argc, char** argv)
{
  clock_t t; t = clock();
  time_t start, end; time(&start);

  int Nevents = atoi(argv[1]); 
  string directory = argv[2];
  string XMLname   = directory + "jetscape_init.xml";
  string outputname= directory + "evolution_result.dat";
  string sigmaname = directory + "sigmaGen.dat";
  

  cout<<endl;
  JetScapeLogger::Instance()->SetInfo(true);
  //JetScapeLogger::Instance()->SetDebug(false);
  //JetScapeLogger::Instance()->SetRemark(false);
  //JetScapeLogger::Instance()->SetVerboseLevel(0);
   
  Show();

  auto jetscape = make_shared<JetScape>(XMLname, Nevents);
  jetscape->SetId("primary");
  jetscape->SetReuseHydro (true);
  jetscape->SetNReuseHydro (Nevents);

  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();
  //auto hydro = make_shared<HydroFromFile> ();
  auto matter = make_shared<Matter> ();
  //auto martini = make_shared<Martini> ();
  auto pythiaGun= make_shared<PythiaGun> ();
  //auto printer = make_shared<PartonPrinter> ();
  auto hadroMgr = make_shared<HadronizationManager> ();
  auto hadro = make_shared<Hadronization> ();
  //auto hadroModule = make_shared<ColoredHadronization> ();
  auto colorless = make_shared<ColorlessHadronization> ();
  auto writer= make_shared<JetScapeWriterAscii> (outputname);
#ifdef USE_HDF5
  auto initial = make_shared<InitialFromFile>();
  jetscape->Add(initial);
#endif
  jetscape->Add(pythiaGun);
  //jetscape->Add(hydro);
  jloss->Add(matter);
  //jloss->Add(martini);
  jlossmanager->Add(jloss);
  jetscape->Add(jlossmanager);
  //jetscape->Add(printer);
  hadro->Add(colorless);
  hadroMgr->Add(hadro);
  jetscape->Add(hadroMgr);
  jetscape->Add(writer);

  // Intialize all modules tasks
  jetscape->Init();
  jetscape->Exec();
  jetscape->Finish();
  
  INFO_NICE<<"Finished!";
  cout<<endl;

// Some information is only known after the full run,
  // Therefore store information at the end of the file, in a footer
  writer->WriteComment ( "EVENT GENERATION INFORMATION" );
  Pythia8::Info& info = pythiaGun->info;
  std::ostringstream oss;
  oss.str(""); oss << "nTried    = " << info.nTried();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "nSelected = " << info.nSelected();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "nAccepted = " << info.nAccepted();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "sigmaGen  = " << info.sigmaGen();  
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "sigmaErr  = " << info.sigmaErr();
  writer->WriteComment ( oss.str() );

  oss.str(""); oss << "eCM  = " << info.eCM();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "pTHatMin  = " << pythiaGun->GetpTHatMin();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "pTHatMax  = " << pythiaGun->GetpTHatMax();
  writer->WriteComment ( oss.str() );
  writer->WriteComment ( "/EVENT GENERATION INFORMATION" );

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));

  ofstream sigmaGen(sigmaname, std::ios::out);
  sigmaGen << info.sigmaGen() << "\t" << info.sigmaErr() << endl;
  sigmaGen.close();

  // Print pythia statistics
  // pythiaGun->stat();

  // Demonstrate how to work with pythia statistics
  //Pythia8::Info& info = pythiaGun->info;
  cout << " nTried    = " << info.nTried() << endl;
  cout << " nSelected = " << info.nSelected() << endl;
  cout << " nAccepted = " << info.nAccepted() << endl;
  cout << " sigmaGen  = " << info.sigmaGen()  << endl;  
  cout << " sigmaErr  = " << info.sigmaErr()  << endl;
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"------------------------------------------------------";
  INFO_NICE<<"| Jet in hydro from file Test JetScape Framework ... |";
  INFO_NICE<<"------------------------------------------------------";
  INFO_NICE;
}
