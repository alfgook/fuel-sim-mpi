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
#ifndef NOT_USING_MPI
#include "G4MPImanager.hh"
#endif

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
AnalysisMPI::Book(G4bool aMergeNtuple)
{
  G4cout << "AnalysisMPI::Book start" << G4endl;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName("fuel-sim");

  #ifdef G4MULTITHREADED
  // MT ntuple merging
  if(aMergeNtuple) analysisManager->SetNtupleMerging(true);
  #endif

//----------------histograms---------------------------------------------------
  analysisManager->CreateH1("sGamma","Gamma-ray source",2000,0,20); // 0
  analysisManager->CreateH1("sNeutron","Neutron source",2000,0,20); // 1
  analysisManager->CreateH1("sAlpha","alpha source",1000,0,10);     // 2
  
  analysisManager->CreateH1("sGammaOut","Gamma-rays exiting cask",2000,0,20); // 3
  analysisManager->CreateH1("sNeutronOut","Neutrons exiting cask",2000,0,20); // 4

  analysisManager->CreateH1("hNeutronScattered","Spectrum of scattered neutrons",2000,0,20); // 5
  analysisManager->CreateH1("hNeutronNoScattered","Spectrum of non-scattered neutrons",2000,0,20); // 6
  analysisManager->CreateH1("hGammaScattered","Spectrum of scattered gammas",2000,0,20); // 7
  analysisManager->CreateH1("hGammaNonScattered","Spectrum of non-scattered gammas",2000,0,20); // 8
  analysisManager->CreateH1("hToFgn_NonScatter","tof Spectrum of non-scattered events",300,-10,290); // 9
  analysisManager->CreateH1("hToFgn_Scatter","tof Spectrum of scattered events",300,-10,290); // 10
  analysisManager->CreateH1("hNeutronMultiplication","hNeutronMultiplication",2,-0.5,1.5); // 11

  analysisManager->CreateH2("InitPosNonScattered","initial source positions for non-scattered neutrons",380,-380,380,380,-380,380); //0
//----------------root tree / ntuple---------------------------------------
  analysisManager->CreateNtuple("SplitPoints", "SplitPoints");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleIColumn("PDGcode");
  analysisManager->CreateNtupleFColumn("Time");
  analysisManager->CreateNtupleFColumn("Energy");
  analysisManager->CreateNtupleFColumn("posX");
  analysisManager->CreateNtupleFColumn("posY");
  analysisManager->CreateNtupleFColumn("posZ");
  analysisManager->CreateNtupleFColumn("dirX");
  analysisManager->CreateNtupleFColumn("dirY");
  analysisManager->CreateNtupleFColumn("dirZ");
  analysisManager->CreateNtupleFColumn("InitX");
  analysisManager->CreateNtupleFColumn("InitY");
  analysisManager->CreateNtupleFColumn("InitZ");
  analysisManager->CreateNtupleFColumn("weight");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("FUEL-PIN", "ResponseNtuple");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleIColumn("DetectorNbr");
  analysisManager->CreateNtupleIColumn("PDGcode");
  analysisManager->CreateNtupleFColumn("TimeOfEvent");
  analysisManager->CreateNtupleFColumn("TimeInEvent");
  analysisManager->CreateNtupleFColumn("Light");
  analysisManager->CreateNtupleFColumn("Weight");
  analysisManager->CreateNtupleFColumn("InitX");
  analysisManager->CreateNtupleFColumn("InitY");
  analysisManager->CreateNtupleFColumn("InitZ");
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
  G4cout << "AnalysisMPI::Openfile start" << G4endl;
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  G4String filename = mgr->GetFileName();
  //G4String filename = "test";
  G4cout << "AnalysisMPI::filename " << filename << G4endl;
  #ifndef NOT_USING_MPI
  G4int rank = G4MPImanager::GetManager()->GetRank();
  filename += "_r" + std::to_string(rank);
  #endif
  //if(rank==0) mgr->OpenFile(filename.c_str());
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
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  const G4int ntupleID = 1;
  analysisManager->FillNtupleIColumn(ntupleID, 0, eventID);
  analysisManager->FillNtupleIColumn(ntupleID, 1, copyNbr);
  analysisManager->FillNtupleIColumn(ntupleID, 2, PDGcode);
  analysisManager->FillNtupleFColumn(ntupleID, 3, 0.);
  analysisManager->FillNtupleFColumn(ntupleID, 4, time);
  analysisManager->FillNtupleFColumn(ntupleID, 5, light);
  analysisManager->FillNtupleFColumn(ntupleID, 6, weight);
  analysisManager->FillNtupleFColumn(ntupleID, 7, initX);
  analysisManager->FillNtupleFColumn(ntupleID, 8, initY);
  analysisManager->FillNtupleFColumn(ntupleID, 9, initZ);
  analysisManager->AddNtupleRow(ntupleID);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
AnalysisMPI::FillSplitEvent(G4int eventID, G4int PDGcode, G4double time, G4double KinE, G4double posX, G4double posY, G4double posZ, G4double dirX, G4double dirY, G4double dirZ, G4double weight)
{
  const G4int ntupleID = 0;
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleIColumn(ntupleID, 0, eventID);
  analysisManager->FillNtupleIColumn(ntupleID, 1, PDGcode);
  analysisManager->FillNtupleFColumn(ntupleID, 2, time);
  analysisManager->FillNtupleFColumn(ntupleID, 3, KinE);
  analysisManager->FillNtupleFColumn(ntupleID, 4, posX);
  analysisManager->FillNtupleFColumn(ntupleID, 5, posY);
  analysisManager->FillNtupleFColumn(ntupleID, 6, posZ);
  analysisManager->FillNtupleFColumn(ntupleID, 7, dirX);
  analysisManager->FillNtupleFColumn(ntupleID, 8, dirY);
  analysisManager->FillNtupleFColumn(ntupleID, 9, dirZ);
  analysisManager->FillNtupleFColumn(ntupleID, 10, initX);
  analysisManager->FillNtupleFColumn(ntupleID, 11, initY);
  analysisManager->FillNtupleFColumn(ntupleID, 12, initZ);
  analysisManager->FillNtupleFColumn(ntupleID, 13, weight);
  analysisManager->AddNtupleRow(ntupleID);
}