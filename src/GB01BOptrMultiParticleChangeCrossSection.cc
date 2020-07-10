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
/// \file GB01/src/GB01BOptrMultiParticleChangeCrossSection.cc
/// \brief Implementation of the GB01BOptrMultiParticleChangeCrossSection class
//
#include "GB01BOptrMultiParticleChangeCrossSection.hh"
#include "G4BiasingProcessInterface.hh"

#include "GB01BOptrChangeCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB01BOptrMultiParticleChangeCrossSection::GB01BOptrMultiParticleChangeCrossSection(G4String aName)
  : G4VBiasingOperator("TestManyExponentialTransform")
{
  fName = aName;
  //cmd.SetParameterName("factor",false);
  //cmd.GetParameter(1)->SetParameterName("biasfactor");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB01BOptrMultiParticleChangeCrossSection::AddParticle(G4String particleName, G4double biasFactor)
{
  const G4ParticleDefinition* particle =
    G4ParticleTable::GetParticleTable()->FindParticle( particleName );
  
  if ( particle == 0 )
    {
      G4ExceptionDescription ed;
      ed << "Particle `" << particleName << "' not found !" << G4endl;
      G4Exception("GB01BOptrMultiParticleChangeCrossSection::AddParticle(...)",
                  "exGB01.02",
                  JustWarning,
                  ed);
      return;
    }
  
  GB01BOptrChangeCrossSection* optr = new GB01BOptrChangeCrossSection(particleName,biasFactor);
  fParticlesToBias.push_back( particle );
  fBOptrForParticle[ particle ] = optr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB01BOptrMultiParticleChangeCrossSection::ChangeBiasForParticle(G4String particleName, G4double biasFactor)
{
  G4ParticleDefinition* particle =
    G4ParticleTable::GetParticleTable()->FindParticle( particleName );
  
  if ( particle == 0 )
    {
      G4ExceptionDescription ed;
      ed << "Particle `" << particleName << "' not found !" << G4endl;
      G4Exception("GB01BOptrMultiParticleChangeCrossSection::ChangeBiasForParticle(...)",
                  "exGB01.03",
                  JustWarning,
                  ed);
      return;
    }

  G4bool found = false;
  for( auto it=fBOptrForParticle.begin();it!=fBOptrForParticle.end(); ++ it) {
    //G4ParticleDefinition* key = it->first;
    if(it->first==particle) {
      found = true;
      break;
    }
  }
  
  //if(fBOptrForParticle.contains(particle)) { //why doesn't std::map::contains() work ???
  if(found) {
    GB01BOptrChangeCrossSection* optr = fBOptrForParticle.at(particle);
    optr->SetBiasFactor(biasFactor);
  } else {
    G4ExceptionDescription ed;
      ed << "XS bias for particle '" << particleName << "' has not been initialised!" << G4endl;
      G4Exception("GB01BOptrMultiParticleChangeCrossSection::ChangeBiasForParticle(...)",
                  "exGB01.04",
                  JustWarning,
                  ed);
      return;
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VBiasingOperation*
GB01BOptrMultiParticleChangeCrossSection::
ProposeOccurenceBiasingOperation(const G4Track* track,
                                 const G4BiasingProcessInterface* callingProcess)
{
  // -- examples of limitations imposed to apply the biasing:
  // -- limit application of biasing to primary particles only:
  //if ( track->GetParentID() != 0 ) return 0;
  // -- limit to at most 5 biased interactions:
  //if ( fnInteractions > 4 )        return 0;
  // -- and limit to a weight of at least 0.05:
  /*const G4double minWeight = 1.E-08;
  if ( track->GetWeight() < 0.05 ) {
    G4cout << "track weight " << track->GetWeight() << " is below "<< minWeight<< "biasing not applied" << G4endl;
    return 0;
  }*/

  //G4cout << "ProposeOccurenceBiasingOperation track weight " << track->GetWeight() << G4endl;
  
  if ( fCurrentOperator ) return fCurrentOperator->
                            GetProposedOccurenceBiasingOperation(track, callingProcess);
  else                    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB01BOptrMultiParticleChangeCrossSection::StartTracking( const G4Track* track )
{
  // -- fetch the underneath biasing operator, if any, for the current particle type:
  const G4ParticleDefinition* definition = track->GetParticleDefinition();
  std::map < const G4ParticleDefinition*, GB01BOptrChangeCrossSection* > :: iterator
    it = fBOptrForParticle.find( definition );
  fCurrentOperator = 0;
  if ( it != fBOptrForParticle.end() ) fCurrentOperator = (*it).second;

  // -- reset count for number of biased interactions:
  fnInteractions = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
GB01BOptrMultiParticleChangeCrossSection::
OperationApplied( const G4BiasingProcessInterface*               callingProcess, 
                  G4BiasingAppliedCase                              biasingCase,
                  G4VBiasingOperation*                occurenceOperationApplied, 
                  G4double                        weightForOccurenceInteraction,
                  G4VBiasingOperation*               finalStateOperationApplied, 
                  const G4VParticleChange*               particleChangeProduced )
{
  // -- count number of biased interactions:
  fnInteractions++;

  //G4cout << "GB01BOptrMultiParticleChangeCrossSection fnInteractions = " << fnInteractions << G4endl;
  // -- inform the underneath biasing operator that a biased interaction occured:
  if ( fCurrentOperator ) fCurrentOperator->ReportOperationApplied( callingProcess,
                                                                    biasingCase,
                                                                    occurenceOperationApplied,
                                                                    weightForOccurenceInteraction,
                                                                    finalStateOperationApplied,
                                                                    particleChangeProduced );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
