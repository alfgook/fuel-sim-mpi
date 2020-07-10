#include "MyStackingAction.hh"

#include "EventAction.hh"
#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyStackingAction::MyStackingAction(EventAction *aEventAction)
{
	fEventAction = aEventAction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyStackingAction::~MyStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
MyStackingAction::ClassifyNewTrack(const G4Track* track)
{
	/*G4cout << "MyStackingAction::ClassifyNewTrack(...) : " <<fEventAction->GetAccumulatedWeight()<< G4endl;
	if(fEventAction->GetAccumulatedWeight()<1.E-15) { //play russian roulette on events below a threshold weight
			G4cout << "russian roulette " << fEventAction->GetAccumulatedWeight() << G4endl;
		if(G4UniformRand()>=0.9) {
			G4cout << "aborting event with weight " << fEventAction->GetAccumulatedWeight() << G4endl;
			stackManager->clear(); // remove all tracks from the stacks
			return fKill;
		} else {
			fEventAction->AccumulateWeight(10.);
		}
	}*/

	//kill secondary neutrino
	if (track->GetDefinition() == G4NeutrinoE::NeutrinoE()) return fKill;

	if (track->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE()) return fKill;
	//keep primary particle
	if (track->GetParentID() == 0) return fUrgent;

	if(track->GetDefinition() == G4Gamma::Gamma()) return fUrgent;

	if(track->GetDefinition() == G4Neutron::Neutron()) return fUrgent;

	else return fWaiting;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
