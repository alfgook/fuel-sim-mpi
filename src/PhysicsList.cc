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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//


#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4DecayPhysics.hh"
#include "G4NuclideTable.hh"
#include "G4StepLimiterPhysics.hh"
//#include "BiasedRDPhysics.hh"
#include "MyRadioactiveDecayPhysics.hh"
//#include "G4RadioactiveDecayPhysics.hh"
#include "G4GenericBiasingPhysics.hh"

#include "NeutronHPphysics.hh"

#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonPhysicsPHP_new.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsPHP.hh"

//#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4RegionStore.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"


PhysicsList::PhysicsList(G4bool np=true)
:G4VModularPhysicsList()
{
  G4int verb = 0;
  SetVerboseLevel(verb);
  
  //add new units for radioActive decays
  //
  new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
  // 
  const G4double minute = 60*second;
  const G4double hour   = 60*minute;
  const G4double day    = 24*hour;
  const G4double year   = 365*day;
  new G4UnitDefinition("minute", "min", "Time", minute);
  new G4UnitDefinition("hour",   "h",   "Time", hour);
  new G4UnitDefinition("day",    "d",   "Time", day);
  new G4UnitDefinition("year",   "y",   "Time", year);


  G4int verbose = 0;

  // Mandatory for G4NuclideTable
  // Half-life threshold must be set small or many short-lived isomers 
  // will not be assigned life times (default to 0) 
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

  // Radioactive decay: I don't use it decay is done in the Primary Generator Action

#ifndef DEBUG_ALPHA_N

  // Hadron Physics
  //RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verbose)); //try to add on top

  // Neutron Physics
  if(np) RegisterPhysics( new NeutronHPphysics("neutronHP")); //using FREYA
#endif
  // Stopping Physics
  //RegisterPhysics( new G4StoppingPhysics());

  // Ion Physics
  //RegisterPhysics( new G4IonPhysicsPHP_new(verbose)); //this is the physics list from the SaG4n application-package from ciemat

  // EM Physics
  //RegisterPhysics( new G4EmStandardPhysics_option4(verbose) ); // photons, electrons et al
  //RegisterPhysics( new G4EmExtraPhysics(verbose) ); // neutrinos et al
  //RegisterPhysics( new G4EmStandardPhysics() );
  //RegisterPhysics( new G4EmPenelopePhysics() );
  //RegisterPhysics( new G4EmLivermorePhysics() );

  //step limit
  //RegisterPhysics( new G4StepLimiterPhysics(verbose) );
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
 
  SetCutsWithDefault();

  // Production thresholds for detector regions
  G4Region* region;
  G4String regName;
  G4ProductionCuts* cuts;
  regName = "FuelRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(1.*cm); // same cuts for gamma, e- and e+ region->SetProductionCuts(cuts);
  region->SetProductionCuts(cuts);

  G4ProductionCuts* cutsDetectors;
  regName = "DetectorRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cutsDetectors = new G4ProductionCuts;
  cutsDetectors->SetProductionCut(0.01*mm); // same cuts for gamma, e- and e+ region->SetProductionCuts(cuts);
  region->SetProductionCuts(cutsDetectors);
  

  if (verboseLevel > 0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
