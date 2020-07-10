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
/// \file ParticleSourceFromFile.hh
/// \brief Definition of the ParticleSourceFromFile class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ParticleSourceFromFile_h
#define ParticleSourceFromFile_h 1

#include "G4VPrimaryGenerator.hh"
#include <vector>
#include "G4ThreeVector.hh"


class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct TreeEntry
{
        G4int eventID;
        G4int PDGcode;
        G4double Time;
        G4double Energy;
        G4double posX;
        G4double posY;
        G4double posZ;
        G4double dirX;
        G4double dirY;
        G4double dirZ;
        G4double weight;
};

class ParticleSourceFromFile : public G4VPrimaryGenerator
{
  public:
    ParticleSourceFromFile(G4String);    
   ~ParticleSourceFromFile();

  public:
    virtual void GeneratePrimaryVertex(G4Event*);
    G4bool GetNextEvent();
    G4bool EventHasGamma();
    void OpenFile(G4String);
    G4String GetFileName() const;
    void Print(G4int) const;
    void SetSplitFactor(G4int);

  private:
    G4int fNbrGenerated;
    G4int fSplitFactor;
    G4double fWeight;

    TreeEntry currentEntry;
    std::vector<TreeEntry> TreeEntries;

    G4bool eofReached;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
