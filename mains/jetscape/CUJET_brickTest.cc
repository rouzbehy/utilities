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
// JetScape Framework Brick Test Program
// (use either shared library (need to add paths; see setup.csh)
// (or create static library and link in)
// -------------------------------------------------------------

#include <iostream>
#include <time.h>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

// User modules derived from jetscape framework clasess
#include "TrentoInitial.h"
#include "AdSCFT.h"
#include "Matter.h"
#include "LBT.h"
#include "Martini.h"
#include "CUJET.h"
#include "Brick.h"
#include "GubserHydro.h"
#include "PGun.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"
#include "CausalLiquefier.h"

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
  cout<<endl;
  JetScapeLogger::Instance()->SetInfo(true);
   
  Show();

  auto jetscape = make_shared<JetScape>("./jetscape_init.xml", 150000);
  jetscape->SetId("primary");
  
  // Initial conditions and hydro
  auto trento = make_shared<TrentoInitial>();
  jetscape->Add(trento);
  auto hydro = make_shared<Brick> ();
  jetscape->Add(hydro);
  auto pGun= make_shared<PGun> ();
  jetscape->Add(pGun);

  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();

  //auto myliquefier = make_shared<CausalLiquefier> ();
  //jloss->add_a_liquefier(myliquefier);
  // Note: if you use Matter, it MUST come first (to set virtuality)
  //auto matter = make_shared<Matter> ();
  //jloss->Add(matter);

  // auto adscft = make_shared<AdSCFT> ();
  //jloss->Add(adscft);

  //auto lbt = make_shared<LBT> ();
  //jloss->Add(lbt);

  //auto martini = make_shared<Martini> ();
  //jloss->Add(martini);

  auto cujet = make_shared<CUJET> ();
  jloss->Add(cujet);

  jlossmanager->Add(jloss);  
  jetscape->Add(jlossmanager);

  auto printer = make_shared<PartonPrinter> () ;
  jetscape->Add(printer);
    
  // Hadronization
 // auto hadroMgr = make_shared<HadronizationManager> ();
 // auto hadro = make_shared<Hadronization> ();
 // auto hadroModule = make_shared<ColoredHadronization> ();
// hadro->Add(hadroModule);

  // auto colorless = make_shared<ColorlessHadronization> ();
  // hadro->Add(colorless);
 // hadroMgr->Add(hadro);
 // jetscape->Add(hadroMgr);

  // Output
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  // same as JetScapeWriterAscii but gzipped
  // auto writer= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");
  // HEPMC3
#ifdef USE_HEPMC
  // auto writer= make_shared<JetScapeWriterHepMC> ("test_out.hepmc");
#endif
  jetscape->Add(writer);

  // Intialize all modules tasks
  jetscape->Init();

  // Run JetScape with all task/modules as specified
  jetscape->Exec();

  // For the future, cleanup is mostly already done in write and clear
  jetscape->Finish();
  
  INFO_NICE<<"Finished!";
  cout<<endl;

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"------------------------------------";
  INFO_NICE<<"| Brick Test JetScape Framework ... |";
  INFO_NICE<<"------------------------------------";
  INFO_NICE;
}
