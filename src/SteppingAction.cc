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
// $Id: SteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "AnalysisMPI.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Tubs.hh"
#include "G4TwoVector.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4TrackVector.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{
	fCaskPV = G4PhysicalVolumeStore::GetInstance()->GetVolume("ShieldingWallPV");
  fWorldPV = G4PhysicalVolumeStore::GetInstance()->GetVolume("DetectionVolumePV");

  fSplittingActive = true;
	//timer.Start();

  fMessenger = new G4GenericMessenger(this,"/EventSplitting/", "...doc...");

  fMessenger->DeclareProperty("on-off",
                            fSplittingActive,
                            "1 = splitting is on : 0 = splitting is off");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
	//timer.Stop();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{   
  auto postStepPoint = aStep->GetPostStepPoint();
  auto CurrentProcess = postStepPoint->GetProcessDefinedStep();

  if(!CurrentProcess) {
    return;
  }

  if(CurrentProcess->GetProcessName()=="parallelWorldForSplitting" && fSplittingActive) {

    //G4ThreeVector PostStepPos = postStepPoint->GetPosition();
    // save the track information to file

    G4double kinE = postStepPoint->GetKineticEnergy();
    if(kinE) { //ions that stop inside the parallel geometry also invokes this process

      AnalysisMPI *analysis = AnalysisMPI::GetAnalysis();
      analysis->FillSplitEvent(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(),
                                aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding(),
                                postStepPoint->GetGlobalTime(),
                                postStepPoint->GetKineticEnergy(),
                                postStepPoint->GetPosition().x(),
                                postStepPoint->GetPosition().y(),
                                postStepPoint->GetPosition().z(),
                                postStepPoint->GetMomentumDirection().x(),
                                postStepPoint->GetMomentumDirection().y(),
                                postStepPoint->GetMomentumDirection().z(),
                                postStepPoint->GetWeight());
      /*const G4int ntupleID = 0;
      auto analysisManager = G4AnalysisManager::Instance();
      analysisManager->FillNtupleIColumn(ntupleID, 0, G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
      analysisManager->FillNtupleIColumn(ntupleID, 1, aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
      analysisManager->FillNtupleDColumn(ntupleID, 2, postStepPoint->GetGlobalTime());
      analysisManager->FillNtupleDColumn(ntupleID, 3, postStepPoint->GetKineticEnergy());
      analysisManager->FillNtupleDColumn(ntupleID, 4, postStepPoint->GetPosition().x() );
      analysisManager->FillNtupleDColumn(ntupleID, 5, postStepPoint->GetPosition().y() );
      analysisManager->FillNtupleDColumn(ntupleID, 6, postStepPoint->GetPosition().z() );
      analysisManager->FillNtupleDColumn(ntupleID, 7, postStepPoint->GetMomentumDirection().x() );
      analysisManager->FillNtupleDColumn(ntupleID, 8, postStepPoint->GetMomentumDirection().y() );
      analysisManager->FillNtupleDColumn(ntupleID, 9, postStepPoint->GetMomentumDirection().z() );
      analysisManager->FillNtupleDColumn(ntupleID, 10, postStepPoint->GetWeight() );
      analysisManager->AddNtupleRow(ntupleID);*/

    }
    G4Track *track = aStep->GetTrack();
    track->SetTrackStatus(fStopAndKill);

  }

//------- scoring gammas/neutrons that make it through the cask wall -------------
	/*G4int PDGcode = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
	if(PDGcode!=22 && PDGcode!=2112) return;
	G4VPhysicalVolume *preStepVolume = aStep->GetPreStepPoint()->GetPhysicalVolume();
	G4VPhysicalVolume *postStepVolume = aStep->GetPostStepPoint()->GetPhysicalVolume();

	if(postStepVolume==fWorldPV && preStepVolume==fCaskPV) {
		//G4cout << "exit" << G4endl;
		G4double Energy = aStep->GetPostStepPoint()->GetKineticEnergy();
		G4double Weight = aStep->GetTrack()->GetWeight();
		//G4cout << "  energy = " << Energy << G4endl;
		//G4cout << "  weight = " << Weight << G4endl;

		//G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		if(PDGcode==22) analysisManager->FillH1(3,Energy,Weight);
		if(PDGcode==2112) analysisManager->FillH1(4,Energy,Weight);
	}*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
