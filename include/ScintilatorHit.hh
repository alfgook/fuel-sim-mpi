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
// $Id: ScintilatorHit.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file ScintilatorHit.hh
/// \brief Definition of the ScintilatorHit class

#ifndef ScintilatorHit_h
#define ScintilatorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class ScintilatorHit : public G4VHit
{
  public:
    ScintilatorHit();
    ScintilatorHit(const ScintilatorHit&);
    ScintilatorHit(G4int, G4int, G4int, G4double);
    virtual ~ScintilatorHit();

    // operators
    const ScintilatorHit& operator=(const ScintilatorHit&);
    G4int operator==(const ScintilatorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

	//
	//ScintilatorHit *clone(G4int, G4int)

    // Set methods
    void AddLight(G4double);
    void SetPDGcode(G4int aPDGcode) {fPDGcode=aPDGcode;};
    void SetTime(G4double aTime) { fTime = aTime; };

    // Get methods
    G4int GetParentID() const     { return fParentID; };
    G4int GetPDGcode() const     { return fPDGcode; };
    G4double GetLight() const     { return fLight; };
    //G4double GetTime() const     { if(fLight){ return fTime/fLight; }else{ return 0.; } };
    //G4double GetTime() const     { if(fLight){ return fTime; }else{ return 0.; } };
    G4double GetTime() const     { return fTime; };
	G4double GetWeight() const  { return fWeight; };
	G4int GetVolCopyNo() const { return fVolumeCopyNo; };

  private:

	G4int				fVolumeCopyNo;
	G4int				fParentID;
	G4int         		fPDGcode;
	G4double      		fLight;
	G4double      		fTime;
	G4double      		fWeight;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<ScintilatorHit> ScintilatorHitsCollection;

extern G4ThreadLocal G4Allocator<ScintilatorHit>* ScintilatorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* ScintilatorHit::operator new(size_t)
{
  if(!ScintilatorHitAllocator)
      ScintilatorHitAllocator = new G4Allocator<ScintilatorHit>;
  return (void *) ScintilatorHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void ScintilatorHit::operator delete(void *hit)
{
  ScintilatorHitAllocator->FreeSingle((ScintilatorHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
