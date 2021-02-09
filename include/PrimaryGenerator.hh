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
/// \file PrimaryGenerator.hh
/// \brief Definition of the PrimaryGenerator class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGenerator_h
#define PrimaryGenerator_h 1

#include "G4VPrimaryGenerator.hh"
#include <vector>


class G4Event;
class G4SPSPosDistribution;
class G4SPSRandomGenerator;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGenerator : public G4VPrimaryGenerator
{
  public:
    PrimaryGenerator();    
   ~PrimaryGenerator();

  public:
    virtual void GeneratePrimaryVertex(G4Event*);
    virtual void GeneratePrimaryVertexFromSpectrum(G4Event*);
    void GenerateDecayPos();
    void SetWeight(G4double aWeight) { fWeight = aWeight; };
    void AddParticle(G4ParticleDefinition* aParticle, G4double kinEnergy, G4ThreeVector direction);
    void SetTime(G4double aTime) { fTime = aTime; };
    void ReadSpectrum(G4String FileName);

  private:
    G4SPSPosDistribution	*fPosDist;
    G4SPSRandomGenerator	*fPosGenerator;
    G4ThreeVector decayPos;
    G4int nAssemblies;

    G4double        fTime;
  	G4double				fWeight;
  	G4int					nParticles;
  	std::vector<G4double>				fKinEnergy;
  	std::vector<G4ThreeVector>			fDirection;
    std::vector<G4ParticleDefinition*>	theParticles;
    G4ThreeVector assemblyPos[12];

    G4ParticleDefinition  *theSpectrumParticle;
    std::vector<G4double> theSpectrumI;
    std::vector<G4double> theSpectrumE;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
