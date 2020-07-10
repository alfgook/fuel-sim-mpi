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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4SaG4nParticleHPFSFissionFS_h
#define G4SaG4nParticleHPFSFissionFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4DynamicParticleVector.hh"
#include "G4SaG4nParticleHPFinalState.hh"
#include "G4SaG4nParticleHPNames.hh"
#include "G4SaG4nParticleHPParticleYield.hh"
#include "G4SaG4nParticleHPVector.hh"
#include "G4SaG4nParticleHPFissionERelease.hh"
#include "G4SaG4nParticleHPEnergyDistribution.hh"
#include "G4SaG4nParticleHPPhotonDist.hh"
#include "G4SaG4nParticleHPAngular.hh"
#include "G4Cache.hh"

class G4SaG4nParticleHPFSFissionFS : public G4SaG4nParticleHPFinalState
{
  struct toBeCached
  {
    const G4ReactionProduct* theNeutronRP;
    const G4ReactionProduct* theTarget;
    toBeCached() : theNeutronRP(0),theTarget(0) {}
  };

public:
  
  G4SaG4nParticleHPFSFissionFS() { hasXsec = true; }
  ~G4SaG4nParticleHPFSFissionFS() {}
  
  void Init (G4double A, G4double Z, G4int M,
             G4String & dirName, G4String & aFSType, G4ParticleDefinition*);
  
  G4DynamicParticleVector * ApplyYourself(G4int Prompt, G4int delayed,
                                          G4double *decayconst);
  
  G4SaG4nParticleHPFinalState* New()
  {
    G4SaG4nParticleHPFSFissionFS * theNew = new G4SaG4nParticleHPFSFissionFS;
    return theNew;
  }
  
  inline G4double GetMass()
  {
    return theFinalStateNeutrons.GetTargetMass();
  }
  
  void SampleNeutronMult(G4int&all, 
	  		 G4int&Prompt, 
			 G4int&delayed, 
			 G4double energy,
			 G4int off);
						 
  inline void SetNeutronRP(const G4ReactionProduct & aNeutron)
  { 
    fCache.Get().theNeutronRP = &aNeutron; 
    theNeutronAngularDis.SetProjectileRP(aNeutron);
  }
  
  inline void SetTarget(const G4ReactionProduct & aTarget)
  { 
    fCache.Get().theTarget = &aTarget; 
    theNeutronAngularDis.SetTarget(aTarget);
  }
    
  G4DynamicParticleVector * GetPhotons();
  
  inline G4SaG4nParticleHPFissionERelease * GetEnergyRelease()
  {
    return &theEnergyRelease;
  }
  
private:

  G4HadFinalState * ApplyYourself(const G4HadProjectile & ) { return 0; }  
  
  G4SaG4nParticleHPParticleYield theFinalStateNeutrons;
  G4SaG4nParticleHPEnergyDistribution thePromptNeutronEnDis;
  G4SaG4nParticleHPEnergyDistribution theDelayedNeutronEnDis;
  G4SaG4nParticleHPAngular theNeutronAngularDis;
  
  G4SaG4nParticleHPPhotonDist theFinalStatePhotons;
  G4SaG4nParticleHPFissionERelease theEnergyRelease;
  
  G4Cache<toBeCached> fCache;
  G4SaG4nParticleHPNames theNames;  
};

#endif
