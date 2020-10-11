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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "EventAction.hh"
#include "TrackingMessenger.hh"
#include "Analysis.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessType.hh"

#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(EventAction* EA)
:G4UserTrackingAction(),
 fEventAction(EA),fTrackMessenger(0),
 fFullChain(true)
 
{
  fPrimaryTime = 0.;
  fTrackMessenger = new TrackingMessenger(this);   
  
  fTimeWindow1 = fTimeWindow2 = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
   G4cout << "~TrackingAction1" << G4endl;
  delete fTrackMessenger;
   G4cout << "~TrackingAction2" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::SetTimeWindow(G4double t1, G4double dt)
{
  fTimeWindow1 = t1;
  fTimeWindow2 = fTimeWindow1 + dt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  initWeight = track->GetWeight();
  //G4cout << track->GetParticleDefinition()->GetParticleName() << G4endl;
  //fTimer.Start();
  //G4cout << "initWeight = " << initWeight << G4endl;

  if(!track->GetParentID() && track->GetTrackID()==1) {
    // this is to account for the initial weight of the event,
    // should be executed only once. The result of accumulated
    // weight for this track will be w_0*accumulated weight of
    // this track, where w_0 is the initial event weight
    if(fEventAction->GetAccumulatedWeight()==1.) {
      initWeight = 1.;
      fEventAction->SetPrimaryWeight(track->GetWeight());
    }
  }/* else {
    G4double Charge = track->GetDefinition()->GetPDGCharge();
    if(Charge>2) {
      G4cout << "parent track: " << track->GetParentID() << G4endl;
      G4cout << " " << track->GetParticleDefinition()->GetParticleName() << G4endl;
      //G4cout << "creator process: " << track->GetCreatorProcess() << G4endl;
      //G4cout << "creator process type: " << track->GetCreatorProcess()->GetProcessType() << G4endl;
      //G4cout << "creator process type of decay: " << G4ProcessType::fDecay << G4endl;
      G4cout << "creator process name: " << track->GetCreatorProcess()->GetProcessName() << G4endl;
      if(track->GetCreatorProcess()->GetProcessType()!=G4ProcessType::fDecay) {
        // Do not decay nuclei that are not created in a decay chain
        // Radio active nuclei may be knocked out by collision from a neutron
        // And I don't want to decay these
        G4cout << "Killing the secondary: " << track->GetParticleDefinition()->GetParticleName() << G4endl;
        G4Track* tr = (G4Track*) track;
        tr->SetTrackStatus(fStopAndKill); 
      }
      
    }
  }*/

  

  G4int PDGcode = track->GetParticleDefinition()->GetPDGEncoding();
  if(PDGcode==2112 && track->GetKineticEnergy()>19.99*MeV) {
    G4Track *tr = (G4Track*) track;
    tr->SetTrackStatus(fStopAndKill); //kill neutron above 19.99 MeV otherwise the code crashes because I have no other model than neutronHP
  }

  G4double initialEnergy = track->GetVertexKineticEnergy();
  G4double Weight = track->GetWeight();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if(PDGcode==22) analysisManager->FillH1(0,initialEnergy,Weight);          //gammas
  if(PDGcode==2112) analysisManager->FillH1(1,initialEnergy,Weight);        //neutrons
  if(PDGcode==1000020040) analysisManager->FillH1(2,initialEnergy,Weight);  //alphas
	
  #ifdef DEBUG_ALPHA_N
  //this is for setting up the (a,n) biasing, only alpha particles will be tracked
    if(track->GetParticleDefinition()->GetPDGEncoding()!=1000020040) {
      G4Track *tr = (G4Track*) track;
      tr->SetTrackStatus(fStopAndKill);
      if(fpTrackingManager->GetVerboseLevel()) {
        G4cout << "UserTrackingAction: " << tr->GetParticleDefinition()->GetParticleName() << " is not an alpha-particle : killing it!" << G4endl;
      }
    }
  #endif

  //fTimer.Stop();
   //G4cout << "PreUserTrackingAction time (real): " << fTimer.GetRealElapsed() << G4endl;
  //fTimer.Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  G4double finalWeight = track->GetWeight();
  //G4cout << "finalWeight = " << finalWeight << G4endl;
  //G4cout << "----------------------------" << G4endl;
  if(initWeight) accumulatedWeight = finalWeight/initWeight;
  fEventAction->AccumulateWeight(accumulatedWeight);

  //fTimer.Stop();
  //G4cout << "Total tracking time (real): " << fTimer.GetRealElapsed() << G4endl;
  //G4cout << "Total tracking time (user): " << fTimer.GetUserElapsed() << G4endl;
  //G4cout << "Total tracking time (syst): " << fTimer.GetSystemElapsed() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

