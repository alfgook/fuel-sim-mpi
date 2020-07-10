
#ifndef MyStackingAction_h
#define MyStackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"

/// Stacking action class : manage the newly generated particles
///
/// One wishes do not track secondary neutrino.Therefore one kills it 
/// immediately, before created particles will  put in a stack.

class EventAction;

class MyStackingAction : public G4UserStackingAction
{
  public:
    MyStackingAction(EventAction *aEventAction);
    virtual ~MyStackingAction();
     
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*); 

  private:
    EventAction *fEventAction;     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
