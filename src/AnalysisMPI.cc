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
//
/// @file AnalysisMPI.cc
/// @brief Define histograms and ntuples

#include "G4AutoDelete.hh"
#include "G4SystemOfUnits.hh"
#include "AnalysisMPI.hh"
#include "G4MPImanager.hh"

//Select format of output here
//Note: ntuple merging is supported only with Root format
#include "g4root.hh"
#include "G4RootAnalysisManager.hh"


G4ThreadLocal G4int AnalysisMPI::fincidentFlag = false;
G4ThreadLocal AnalysisMPI* the_analysis = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AnalysisMPI*
AnalysisMPI::GetAnalysis()
{
  if (!the_analysis)
    {
      the_analysis = new AnalysisMPI();
      G4AutoDelete::Register(the_analysis);
    }
  return the_analysis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AnalysisMPI::AnalysisMPI()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
AnalysisMPI::Book()
{
  G4cout << "AnalysisMPI::Book start" << G4endl;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  #ifdef G4MULTITHREADED
  // MT ntuple merging
  analysisManager->SetNtupleMerging(true);
  #endif

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

  G4cout << "AnalysisMPI::Book finished " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AnalysisMPI::~AnalysisMPI()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
AnalysisMPI::OpenFile()
{
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  G4String filename = mgr->GetFileName();
  G4int rank = G4MPImanager::GetManager()->GetRank();
  filename += "_r" + std::to_string(rank);
  mgr->OpenFile(filename.c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
AnalysisMPI::Write()
{
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  mgr->Write();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
AnalysisMPI::CloseFile(G4bool reset)
{
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  mgr->CloseFile(reset);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
AnalysisMPI::FillScintillatorHit(G4int eventID, G4int copyNbr, G4int PDGcode, G4double time, G4double light, G4double weight)
{
  /*G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  const G4int ntupleID = 1;
  analysisManager->FillNtupleIColumn(ntupleID, 0, eventID);
  analysisManager->FillNtupleIColumn(ntupleID, 1, copyNbr);
  analysisManager->FillNtupleIColumn(ntupleID, 2, PDGcode);
  analysisManager->FillNtupleDColumn(ntupleID, 3, 0.);
  analysisManager->FillNtupleDColumn(ntupleID, 4, time);
  analysisManager->FillNtupleDColumn(ntupleID, 5, light);
  analysisManager->FillNtupleDColumn(ntupleID, 6, weight);
  analysisManager->AddNtupleRow(ntupleID);*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
AnalysisMPI::FillSplitEvent(G4int eventID, G4int PDGcode, G4double time, G4double KinE, G4double posX, G4double posY, G4double posZ, G4double dirX, G4double dirY, G4double dirZ, G4double weight)
{
  /*const G4int ntupleID = 0;
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleIColumn(ntupleID, 0, eventID);
  analysisManager->FillNtupleIColumn(ntupleID, 1, PDGcode);
  analysisManager->FillNtupleDColumn(ntupleID, 2, time);
  analysisManager->FillNtupleDColumn(ntupleID, 3, KinE);
  analysisManager->FillNtupleDColumn(ntupleID, 4, posX);
  analysisManager->FillNtupleDColumn(ntupleID, 5, posY);
  analysisManager->FillNtupleDColumn(ntupleID, 6, posZ);
  analysisManager->FillNtupleDColumn(ntupleID, 7, dirX);
  analysisManager->FillNtupleDColumn(ntupleID, 8, dirY);
  analysisManager->FillNtupleDColumn(ntupleID, 9, dirZ);
  analysisManager->FillNtupleDColumn(ntupleID, 10, weight);
  analysisManager->AddNtupleRow(ntupleID);*/
}