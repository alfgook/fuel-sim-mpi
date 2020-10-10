//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file rdecay01.cc
/// \brief Main program of the radioactivedecay/rdecay01 example
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "G4MPImanager.hh"
#include "G4MPIsession.hh"
#include "G4MPIextraWorker.hh"

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "ParallelWorldForSplitting.hh"
#include "G4ParallelWorldPhysics.hh"


#include "FTFP_BERT.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "G4GenericBiasingPhysics.hh"

#include "G4PhysListFactory.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
//#include "G4ScoringManager.hh"


#include <unistd.h>
#define GetCurrentDir getcwd
#include <stdio.h>  /* defines FILENAME_MAX */

namespace {
    void PrintUsage() {
        G4cerr << " Usage: " << G4endl;
        G4cerr << " FuelPinSF [-m macro ] [-af filename] [-mf filename] [-t nThreads]" << G4endl;
        G4cerr << "   option [-af filename]: file for the activity of the fuel material " << G4endl;
        G4cerr << "                          if not specified a simple G4ParticleGun will used " << G4endl;
        G4cerr << "   option [-mf filename]: file for the fuel material in MCNP format" << G4endl;
        G4cerr << "   option [-biasing keyword]: if keyword=off the geometrical biasing is ignored" << G4endl;
        G4cerr << "   option [-np keyword]: if keyword=off the neutronphysics is not registered" << G4endl;
        G4cerr << "   note: -t option is available only for multi-threaded mode."
        << G4endl;
    }
}

int main(int argc,char** argv) {

  int argcMacro = 0;
  G4String macro;
  G4String materialFile;
  G4String activityFile;
  G4String onOffBiasing = "on";
  G4bool RegisterNP = true;
#ifdef G4MULTITHREADED
  G4int nThreads = 1;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
      if      ( G4String(argv[i]) == "-m" ) {
        macro = argv[i+1];
        argcMacro = i+1;
      }
      else if ( G4String(argv[i]) == "-af" ) activityFile = argv[i+1];
      else if ( G4String(argv[i]) == "-mf" ) materialFile = argv[i+1];
      else if ( G4String(argv[i]) == "-biasing" ) onOffBiasing = argv[i+1];
      else if ( G4String(argv[i]) == "-np" ) {
        if(G4String(argv[i+1]) == "off") RegisterNP = false;
      }
#ifdef G4MULTITHREADED
      else if ( G4String(argv[i]) == "-t" ) {
          nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
      }
#endif
      else {
          PrintUsage();
          return 1;
      }
  }

  G4cout << "input dir = " << std::getenv("FUEL_SIM_INPUT_DIR") << G4endl;

  G4cout << "macro = " << macro << G4endl;
  G4cout << "biasing = " << onOffBiasing << G4endl;
  G4cout << "materialFile = " << materialFile << G4endl;
  G4cout << "activityFile = " << activityFile << G4endl;

  int argcMPI = 1;
  if(argcMacro) ++argcMPI;
  char **argvMPI = new char*[argcMPI];
  argvMPI[0] = argv[0];
  if (argcMacro) argvMPI[1] = argv[argcMacro];


  // --------------------------------------------------------------------
  // MPI session
  // --------------------------------------------------------------------
  // At first, G4MPImanager/G4MPIsession should be created.
  G4int nofExtraWorkers = 0;
//#ifndef G4MULTITHREADED
//  if ( mergeNtuple ) nofExtraWorkers = 1;
//#endif
  G4MPImanager* g4MPI = new G4MPImanager(argcMPI, argvMPI, nofExtraWorkers);
  g4MPI->SetVerbose(1);
  
  // MPI session (G4MPIsession) instead of G4UIterminal
  // Terminal availability depends on your MPI implementation.
  G4MPIsession* session = g4MPI-> GetMPIsession();

  // LAM/MPI users can use G4tcsh.
  G4String prompt = "[40;01;33m";
  prompt += "G4MPI";
  prompt += "[40;31m(%s)[40;36m[%/][00;30m:";
  session-> SetPrompt(prompt);

  // --------------------------------------------------------------------
  // user application setting
  // --------------------------------------------------------------------

  char cCurrentPath[FILENAME_MAX];
  G4cout << "GetCurrentDir : " << GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)) << G4endl;

  G4cout << "G4PARTICLEHPDATA = " << std::getenv("G4PARTICLEHPDATA") << G4endl;
  //G4cout << "G4ALPHAHPDATA = " << std::getenv("G4ALPHAHPDATA") << G4endl;

  //detect interactive mode (if no arguments) and define UI session
  /*G4UIExecutive* ui = 0;
  if (!macro.size()) {
    G4cout << "UI session = " << macro << G4endl;
    ui = new G4UIExecutive(argc,argv);
  }*/

  //choose the Random engine
  //CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  //G4Random::setTheEngine( &defaultEngine ); 
  G4int seed = time( NULL ); 
  //G4Random::setTheSeed( seed );

  //construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  nThreads = std::min(nThreads,G4Threading::G4GetNumberOfCores());
  runManager->SetNumberOfThreads(nThreads);
  //G4ScoringManager::GetScoringManager();
#else
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
  G4RunManager* runManager = new G4RunManager;
#endif
  G4cout << "number of threads = " << runManager->GetNumberOfThreads() << G4endl;

  //set mandatory initialization classes
  //
  DetectorConstruction *theDetector = new DetectorConstruction(materialFile);

  G4String parallelWorldName = "parallelWorldForSplitting";
  //ParallelWorldForSplitting* parallelWorldForSplitting =
  // -- and "augment" detector geometry with the parallel worlds:
  theDetector->RegisterParallelWorld( new ParallelWorldForSplitting(parallelWorldName) );

  runManager->SetUserInitialization(theDetector);

  //----  My own physics list ----
  //G4VModularPhysicsList* thePhysicsList = new FTFP_BERT;
  PhysicsList *thePhysicsList = new PhysicsList(RegisterNP);
  
  //----  One of the reference physics lists ----

  thePhysicsList->RegisterPhysics(new G4ParallelWorldPhysics(parallelWorldName));
  
  // -- and augment it with biasing facilities:
  G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
  biasingPhysics->BeVerbose();
  if ( onOffBiasing == "off" ) {
      G4cout << "      ************************************************* " << G4endl;
      G4cout << "      ********** processes are not wrapped ************ " << G4endl;
      G4cout << "      ************************************************* " << G4endl;
  } else {

      biasingPhysics->Bias("alpha");
      //biasingPhysics->Bias("gamma");
      //biasingPhysics->Bias("neutron");

      thePhysicsList->RegisterPhysics(biasingPhysics);
      G4cout << "      ********************************************************* "
             << G4endl;
      G4cout << "      ********** processes are wrapped for biasing ************ "
             << G4endl;
      G4cout << "      ********************************************************* "
             << G4endl;
  }
  runManager->SetUserInitialization(thePhysicsList);

  runManager->SetUserInitialization(new ActionInitialization(nThreads,activityFile));

  //initialize G4 kernel
  runManager->Initialize();
/*
  G4VisExecutive* visManager = new G4VisExecutive;
  visManager-> Initialize();
  G4cout << G4endl;
*/
  //initialize visualization
  //G4VisManager* visManager = nullptr;

  //get the pointer to the User Interface manager
  //G4UImanager* UImanager = G4UImanager::GetUIpointer();

  /*if (ui)  {
   //interactive mode
   visManager = new G4VisExecutive;
   visManager->Initialize();
   UImanager->ApplyCommand("/control/execute macros/vis.mac");
   ui->SessionStart();
   delete ui;
  }
  else  {
   //batch mode
   G4String command = "/control/execute ";
   G4String fileName = macro;
   UImanager->ApplyCommand(command+fileName);
  }*/

    // --------------------------------------------------------------------
  // ready for go
  // MPIsession treats both interactive and batch modes.
  // Just start your session as below.
  // --------------------------------------------------------------------
  G4cout << "SessionStart" << G4endl;
  session-> SessionStart();

  // --------------------------------------------------------------------
  // termination
  // --------------------------------------------------------------------
  //delete visManager;
  delete g4MPI;
  delete runManager;

  G4cout << "this is the end" << G4endl;
  return EXIT_SUCCESS;

}

