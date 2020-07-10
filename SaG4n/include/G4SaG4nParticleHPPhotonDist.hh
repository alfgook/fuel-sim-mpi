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
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
//
// 070606 fix for Valgrind error by T. Koi
// 070612 fix memory leaking by T. Koi
// 070615 fix memory leaking by T. Koi
// 080625 fix memory leaking by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4SaG4nParticleHPPhotonDist_h
#define G4SaG4nParticleHPPhotonDist_h 1

#include <fstream>

#include "globals.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4SaG4nParticleHPVector.hh"
#include "G4SaG4nParticleHPLegendreTable.hh"
#include "G4SaG4nParticleHPAngularP.hh"
#include "G4SaG4nParticleHPPartial.hh"
#include "G4SaG4nParticleHPFastLegendre.hh"
#include "G4SaG4nParticleHPInterpolator.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4Gamma.hh"
#include "G4InterpolationManager.hh"
#include "G4Cache.hh"

class G4SaG4nParticleHPPhotonDist
{
public:

  G4SaG4nParticleHPPhotonDist()
      : repFlag( 0 ) 
      , targetMass( 0.0 ) 
      , nDiscrete( 0 ) 
      , isoFlag( 0 )
      , tabulationType( 0 )
      , nDiscrete2( 0 ) 
      , nIso( 0 ) 
      , nPartials( 0 )
      , theInternalConversionFlag( 0 )
      , nGammaEnergies( 0 )
      , theBaseEnergy( 0.0 )
  {

     disType = 0;
     energy = 0;
     theYield = 0;
     thePartialXsec = 0;
     isPrimary = 0;
     theShells = 0;
     theGammas = 0;
     nNeu = 0;
     theLegendre = 0;
     theAngular = 0;
     distribution = 0;
     probs = 0;
     partials = 0;
     actualMult.Put( 0 );

     theLevelEnergies = 0;
     theTransitionProbabilities = 0;
     thePhotonTransitionFraction = 0;

  }

  ~G4SaG4nParticleHPPhotonDist()
  {
     delete [] disType;
     delete [] energy;
     delete [] theYield;
     delete [] thePartialXsec;
     delete [] isPrimary;
     delete [] theShells;
     delete [] theGammas;
     delete [] nNeu;
     delete [] theAngular;
     delete [] distribution;
     delete [] probs;

     if ( theLegendre != 0 )
     {
        for ( G4int i = 0 ; i < (nDiscrete2-nIso) ; i++ )
           if ( theLegendre[i] != 0 ) delete[] theLegendre[i]; 

        delete [] theLegendre;
     }

     if ( partials != 0 ) 
     {
        for ( G4int i = 0 ; i < nPartials ; i++ )
           { delete partials[i]; }

        delete [] partials;
     }

     delete [] theLevelEnergies;
     delete [] theTransitionProbabilities;
     delete [] thePhotonTransitionFraction;
     if (actualMult.Get() != 0) delete actualMult.Get();
  }
  
  G4bool InitMean(std::istream & aDataFile);
    
  void InitAngular(std::istream & aDataFile);
  
  void InitEnergies(std::istream & aDataFile);
  
  void InitPartials(std::istream & aDataFile);
  
  G4ReactionProductVector * GetPhotons(G4double anEnergy);
  
  inline G4double GetTargetMass() {return targetMass;}
  
  inline G4bool NeedsCascade() {return repFlag==2;}
  
  inline G4double GetLevelEnergy() {return theBaseEnergy;}

private:

   G4int repFlag;  //representation as multiplicities or transition probability arrays.
   G4double targetMass;
   
   G4int nDiscrete;  //number of discrete photons 
   G4int * disType;  // discrete, or continuum photons
   G4double * energy;  // photon energies
   G4SaG4nParticleHPVector * theYield; // multiplicity as a function of neutron energy.
   G4SaG4nParticleHPVector theTotalXsec;
   G4SaG4nParticleHPVector * thePartialXsec;
   G4int * isPrimary;
  
   G4int isoFlag; // isotropic or not?
   G4int tabulationType;
   G4int nDiscrete2;
   G4int nIso;
   G4double * theShells;
   G4double * theGammas;
   G4int * nNeu;
   G4InterpolationManager theLegendreManager;
   G4SaG4nParticleHPLegendreTable ** theLegendre;
   G4SaG4nParticleHPAngularP ** theAngular;
   
   G4int * distribution; // not used for the moment.                                 
   G4int nPartials;
   G4SaG4nParticleHPVector *  probs; // probabilities for the partial distributions.
   G4SaG4nParticleHPPartial ** partials; // the partials, parallel to the above

   G4Cache< std::vector<G4int>* > actualMult;
   
    // for transition prob arrays start
   G4int theInternalConversionFlag;
   G4int nGammaEnergies;
   G4double theBaseEnergy;
   G4double * theLevelEnergies;
   G4double * theTransitionProbabilities;
   G4double * thePhotonTransitionFraction;
    // for transition prob arrays end

   G4SaG4nParticleHPFastLegendre theLegend; // fast look-up for leg-integrals
   G4SaG4nParticleHPInterpolator theInt; // interpolation
};

#endif
