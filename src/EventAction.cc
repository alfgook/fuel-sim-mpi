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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "AnalysisMPI.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "PrimaryGeneratorAction.hh"

#include "Randomize.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction()
{
	fScintillatorHCID = -1;
  // Set default print level 
  //G4RunManager::GetRunManager()->SetPrintProgress(10000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ScintilatorHitsCollection* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  ScintilatorHitsCollection* hitsCollection =
		static_cast<ScintilatorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));


  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("B4cEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}   
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  //G4cout << "EventAction::BeginOfEventAction()" << G4endl;
	fAccumulatedEventWeight = 1.;

	/*G4int eventID = evt->GetEventID();
	if((eventID%100000)==0) {
		G4Random::saveEngineStatus();
		G4Random::showEngineStatus();
	}*/

	//fTimer.Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{

	//G4cout << "final weight of event = " << fAccumulatedEventWeight << G4endl;

	G4int eventID = evt->GetEventID();
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	// Get hits collections IDs (only once)
	if ( fScintillatorHCID == -1 ) {
		G4HCtable *HCtable = G4SDManager::GetSDMpointer()->GetHCtable();
		G4String HCname = HCtable->GetHCname(0);
		
	    fScintillatorHCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHC");
	}

	//------Sort the data-----------------------------
	ScintilatorHitsCollection* fHC = GetHitsCollection(fScintillatorHCID, evt);

	for(size_t i=0;i<fHC->entries();++i) { //loop over the hits in the detector
		G4int copyNbr = (*fHC)[i]->GetVolCopyNo();
		if(copyNbr<1 || copyNbr>96) {
			G4cout << "Error!! Found a hit in a detector that does not exist: copyNbr = " << copyNbr << G4endl;
			//continue;
		}
		G4int PDGcode = (*fHC)[i]->GetPDGcode(); 
		G4double light = (*fHC)[i]->GetLight(); 
		G4double weight = (*fHC)[i]->GetWeight();
		G4double TimeInEvent = (*fHC)[i]->GetTime();

		if(light>0.) {
    		AnalysisMPI *analysis = AnalysisMPI::GetAnalysis();
    		analysis->FillScintillatorHit(eventID, copyNbr, PDGcode, TimeInEvent, light, weight);
		}

	}

//==========================================================================================================

	//fTimer.Stop();
  //G4cout << "Total event time : " << fTimer.GetRealElapsed() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


