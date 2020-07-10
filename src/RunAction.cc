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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Analysis.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
:G4UserRunAction()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName("FuelPin");

//----------------histograms---------------------------------------------------
  analysisManager->CreateH1("sGamma","Gamma-ray source",2000,0,20); // 0
  analysisManager->CreateH1("sNeutron","Neutron source",2000,0,20); // 1
  analysisManager->CreateH1("sAlpha","alpha source",1000,0,10);     // 2
  
  analysisManager->CreateH1("sGammaOut","Gamma-rays exiting cask",2000,0,20); // 3
  analysisManager->CreateH1("sNeutronOut","Neutrons exiting cask",2000,0,20); // 4
//----------------root tree / ntuple---------------------------------------
  analysisManager->CreateNtuple("SplitPoints", "SplitPoints");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleIColumn("PDGcode");
  analysisManager->CreateNtupleDColumn("Time");
  analysisManager->CreateNtupleDColumn("Energy");
  analysisManager->CreateNtupleDColumn("posX");
  analysisManager->CreateNtupleDColumn("posY");
  analysisManager->CreateNtupleDColumn("posZ");
  analysisManager->CreateNtupleDColumn("dirX");
  analysisManager->CreateNtupleDColumn("dirY");
  analysisManager->CreateNtupleDColumn("dirZ");
  analysisManager->CreateNtupleDColumn("weight");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("FUEL-PIN", "ResponseNtuple");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleIColumn("DetectorNbr");
  analysisManager->CreateNtupleIColumn("PDGcode");
  analysisManager->CreateNtupleDColumn("TimeOfEvent");
  analysisManager->CreateNtupleDColumn("TimeInEvent");
  analysisManager->CreateNtupleDColumn("Light");
  analysisManager->CreateNtupleDColumn("Weight");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ 
  delete G4AnalysisManager::Instance();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  G4cout << "RunAction::BeginOfRunAction()" << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->OpenFile();
  //
  runTimer.Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
            
 //save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  runTimer.Stop();
  G4cout << "EndOfRunAction: total run time " << runTimer.GetRealElapsed() << " seconds" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
