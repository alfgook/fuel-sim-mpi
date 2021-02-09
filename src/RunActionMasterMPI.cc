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
/// \file RunActionMasterMPI.cc
/// \brief Implementation of the RunActionMasterMPI class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMasterMPI.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisMPI.hh"

#ifndef NOT_USING_MPI
#include "G4MPImanager.hh" 
#include "G4MPIntupleMerger.hh"
#endif

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMasterMPI::RunActionMasterMPI()
:G4UserRunAction()
{
  #ifndef NOT_USING_MPI
	if ( G4MPImanager::GetManager()->GetTotalSize() >= 2 ) {
    // Activate MPI ntuple merging
    // The merger must be created before creating G4AnalysisManager:
    // (= the first call to G4AnalysisManager::Instance())
    // and deleted only at the end of program
    G4int nofReducedNtupleFiles = 0;  
       // Multiple reduced ntuple files are not yet supported 
    G4bool rowWise = false;
    G4bool rowMode = true;
    fMPIntupleMerger = new G4MPIntupleMerger(nofReducedNtupleFiles, rowWise, rowMode);
  }
  #endif

  AnalysisMPI::GetAnalysis()->Book();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMasterMPI::~RunActionMasterMPI()
{ 
  #ifndef NOT_USING_MPI
  delete fMPIntupleMerger;
  #endif
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMasterMPI::BeginOfRunAction(const G4Run*)
{ 
//  G4cout << G4MPImanager::GetManager()->GetRank() << " : " << "RunActionMasterMPI::BeginOfRunAction()" << G4endl;
  AnalysisMPI::GetAnalysis()->OpenFile();
  //
  runTimer.Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMasterMPI::EndOfRunAction(const G4Run*)
{    
 //save histograms
  auto analysisManager = AnalysisMPI::GetAnalysis();
  //if(rank==0) analysisManager->OpenFile();
  analysisManager->Write();
  analysisManager->CloseFile();

  runTimer.Stop();
  #ifndef NOT_USING_MPI
  G4cout << G4MPImanager::GetManager()->GetRank() << " : " << "EndOfRunAction: total run time " << runTimer.GetRealElapsed() << " seconds" << G4endl;
  #else
  G4cout << "EndOfRunAction: total run time " << runTimer.GetRealElapsed() << " seconds" << G4endl;
  #endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
