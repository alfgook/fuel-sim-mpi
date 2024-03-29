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
/// \file RunActionMPI.cc
/// \brief Implementation of the RunActionMPI class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMPI.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisMPI.hh"

#ifndef NOT_USING_MPI
#include "G4MPImanager.hh"
#endif

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMPI::RunActionMPI(G4bool aMergeNtuple)
:G4UserRunAction()
{
//  G4cout << G4MPImanager::GetManager()->GetRank() << " : " <<"RunActionMPI::RunActionMPI()" << G4endl;
  AnalysisMPI::GetAnalysis()->Book(aMergeNtuple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMPI::~RunActionMPI()
{ 

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMPI::BeginOfRunAction(const G4Run*)
{ 
//  G4cout << G4MPImanager::GetManager()->GetRank() << " : " << "RunActionMPI::BeginOfRunAction()" << G4endl;
  AnalysisMPI::GetAnalysis()->OpenFile();
  //
  runTimer.Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMPI::EndOfRunAction(const G4Run*)
{
            
 //save histograms
  auto analysisManager = AnalysisMPI::GetAnalysis();
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
