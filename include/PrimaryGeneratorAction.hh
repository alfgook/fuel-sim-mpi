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
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGenerator.hh"
#include "globals.hh"
#include "MyRadioactiveDecayBase.hh"
#include "ActivityTable.hh"
#include "G4Timer.hh"

class G4Event;
class G4DecayProducts;
class G4GeneralParticleSource;
class G4GenericMessenger;
class ParticleSourceFromFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(G4String);    
   ~PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);


    void SetNumberOfWorkers(G4int nn) { nWorkers=nn; };
    //void UseSimpleGun() { fUseSimpleGun = true; }
    void RestrictDecayToSF() {
        if(!fActivityTable) fActivityTable = new ActivityTable(fFile, fRadDecay);
        fActivityTable->RestrictTo("SF decay");
        restrictedDecay = true;
        fUseSimpleGun = false;
    }
    void RestrictDecayToAlpha() {
        if(!fActivityTable) fActivityTable = new ActivityTable(fFile, fRadDecay);
        fActivityTable->RestrictTo("alpha decay");
        restrictedDecay = true;
        fUseSimpleGun = false;
    }
    void ExcludeAphaAndSF() {
        if(!fActivityTable) fActivityTable = new ActivityTable(fFile, fRadDecay);
        fActivityTable->ExcludeAphaAndSF();
        restrictedDecay = true;
        fUseSimpleGun = false;
    }
    void UseSpectrumFile() {
        fParticleGun->ReadSpectrum(fSpectrumFileName);
        fUseSimpleGun = false;
        fUseSpectrum = true;
    }
    void ShootFromFile(G4String);
    void SetSplitFactor(G4String);

    G4double GetTheTime() { return theTime; };
    G4int GetNumberOfWorkers() { return nWorkers; };

    //G4double GetWeight() { return fWeight; };
            
  private:
    G4DecayProducts* DoDecay(const G4ParticleDefinition& theParticleDef);
    G4DecayProducts* DoDecay(const G4ParticleDefinition& theParticleDef, G4int channelNbr);
    G4DecayProducts* DoDecayBRbias(const G4ParticleDefinition& theParticleDef, G4double &weigth);

    G4GenericMessenger *fMessenger;
    
    MyRadioactiveDecayBase *fRadDecay;
    PrimaryGenerator *fParticleGun;

    G4DecayTable *theDecayTable;

    ActivityTable *fActivityTable;
    G4String fFile;

    G4double theTime;
    G4int nWorkers; //must be taken into account when calculating the time between events

    G4GeneralParticleSource* fSimpleGun;
    G4bool fUseSimpleGun;

    ParticleSourceFromFile* fParticleSourceFromFile;
    G4bool fShootFromFile = false;

    G4String fSpectrumFileName;
    G4bool fUseSpectrum;

    G4Timer fTimer; //for debuging

    G4bool restrictedDecay;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
