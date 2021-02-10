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
/// \file PrimaryGenerator.cc
/// \brief Implementation of the PrimaryGenerator1 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSRandomGenerator.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Geantino.hh"
#include "AnalysisMPI.hh"

#include "G4AutoLock.hh"
namespace {
  G4Mutex PrimaryGeneratorMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::PrimaryGenerator()
: G4VPrimaryGenerator()
{
  fPosDist = new G4SPSPosDistribution();
  fPosGenerator = new G4SPSRandomGenerator();
  fPosDist->SetBiasRndm(fPosGenerator);
  fPosDist->SetPosDisType("Volume");
  fPosDist->SetPosDisShape("Cylinder");

  fPosDist->SetRadius(0.5*8.36*mm); //PWR
  fPosDist->SetHalfZ(0.5*3712.*mm); //PWR

  /*
  fPosDist = new G4SPSPosDistribution();
  fPosGenerator = new G4SPSRandomGenerator();
  fPosDist->SetBiasRndm(fPosGenerator);
  fPosDist->SetPosDisType("Volume");
  fPosDist->SetPosDisShape("Para");
  //fPosDist->SetHalfX(214./2.*mm); //PWR
  //fPosDist->SetHalfY(214./2.*mm); //PWR
  //fPosDist->SetHalfZ(4296./2.*mm); //PWR
  */

  nParticles = 0;
  fKinEnergy.clear();
  fDirection.clear();
  theParticles.clear();
  fKinEnergy.reserve(50);
  fDirection.reserve(50);
  theParticles.reserve(50);
  fWeight =1.;
  fTime = 0.;

  //BWR-assemblies in transpot cask
  nAssemblies = 12;
  G4double InsertDimensionJhalf = 0.5*210.;
  assemblyPos[0] = G4ThreeVector(-InsertDimensionJhalf,3.*InsertDimensionJhalf,0);
  assemblyPos[1] = G4ThreeVector(InsertDimensionJhalf,3.*InsertDimensionJhalf,0);
  assemblyPos[2] = G4ThreeVector(-3.*InsertDimensionJhalf,InsertDimensionJhalf,0);
  assemblyPos[3] = G4ThreeVector(-InsertDimensionJhalf,InsertDimensionJhalf,0);
  assemblyPos[4] = G4ThreeVector(InsertDimensionJhalf,InsertDimensionJhalf,0);
  assemblyPos[5] = G4ThreeVector(3.*InsertDimensionJhalf,InsertDimensionJhalf,0);
  assemblyPos[6] = G4ThreeVector(-3.*InsertDimensionJhalf,-InsertDimensionJhalf,0);
  assemblyPos[7] = G4ThreeVector(-InsertDimensionJhalf,-InsertDimensionJhalf,0);
  assemblyPos[8] = G4ThreeVector(InsertDimensionJhalf,-InsertDimensionJhalf,0);
  assemblyPos[9] = G4ThreeVector(3.*InsertDimensionJhalf,-InsertDimensionJhalf,0);
  assemblyPos[10] = G4ThreeVector(-InsertDimensionJhalf,-3.*InsertDimensionJhalf,0);
  assemblyPos[11] = G4ThreeVector(InsertDimensionJhalf,-3.*InsertDimensionJhalf,0);

  //PWR-assemblies in copper canister
/*  nAssemblies = 4;
  assemblyPos[0] = G4ThreeVector(185,185,0);
  assemblyPos[1] = G4ThreeVector(-185,185,0);
  assemblyPos[2] = G4ThreeVector(-185,-185,0);
  assemblyPos[3] = G4ThreeVector(185,-185,0);
*/
  decayPos = G4ThreeVector(0,0,0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{
  //delete fPosDist;
  //delete fPosGenerator;
}

void PrimaryGenerator::GenerateDecayPos()
{
  const G4double Pitch = 12.4*mm;

  decayPos = fPosDist->GenerateOne();
  G4int rowNbr = G4int(10.*G4UniformRand()) - 4;
  G4int colNbr = G4int(10.*G4UniformRand()) - 4;
  while((rowNbr==1 && colNbr==1) || (rowNbr==0 && colNbr==0)) {
    rowNbr = G4int(10.*G4UniformRand()) - 4;
    colNbr = G4int(10.*G4UniformRand()) - 4;
  }

  //G4int assemblyNbr = 0;
  G4int assemblyNbr = nAssemblies*G4UniformRand();
  while(assemblyNbr>=nAssemblies || assemblyNbr<0) assemblyNbr = nAssemblies*G4UniformRand(); //just for safety

  //G4cout << "rowNbr = " << rowNbr << G4endl;
  //G4cout << "colNbr = " << colNbr << G4endl;
  //G4cout << "assemblyNbr = " << assemblyNbr << G4endl;

  decayPos += G4ThreeVector(rowNbr*Pitch-0.5*Pitch,colNbr*Pitch-0.5*Pitch,0.);
  decayPos += assemblyPos[assemblyNbr];

  AnalysisMPI::GetAnalysis()->SetInitial(decayPos.x(),decayPos.y(),decayPos.z());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{

  //G4double time1 = 0.*s;
  G4PrimaryVertex *theVertex = new G4PrimaryVertex(decayPos,fTime);

  if(!nParticles) {
    //G4cout << "no particles in the event shooting a geantino" << G4endl;
    //G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
    G4ParticleDefinition* particleDefinition = G4Geantino::Geantino();
    this->AddParticle(particleDefinition,10.*MeV,G4ThreeVector(0,0,1));
  }
  
  /* The bug doesn't appear to be in here!!!
  G4cout << "nParticles = " << nParticles << G4endl;
  G4cout << "theParticles.size() = " << theParticles.size() << G4endl;
  G4cout << "fDirection.size() = " << fDirection.size() << G4endl;
  G4cout << "fKinEnergy.size() = " << fKinEnergy.size() << G4endl;
  //---Test to debug-------
  if(nParticles!=theParticles.size()) {
    G4cout << "error nParticles!=theParticles.size()" << G4endl;
  }
  if(nParticles!=fDirection.size()) {
    G4cout << "error nParticles!=fDirection.size()" << G4endl;
  }
  if(nParticles!=fKinEnergy.size()) {
    G4cout << "error nParticles!=fKinEnergy.size()" << G4endl;
  }
  //--end Test to debug-------
  */

  for(G4int i=0;i<nParticles;i++) {
    G4PrimaryParticle* particle = new G4PrimaryParticle(theParticles.at(i));
    particle->SetMomentumDirection(fDirection.at(i));
    particle->SetKineticEnergy(fKinEnergy.at(i));
    particle->SetWeight(fWeight);
    //particle->SetProperTime(fTime);
    theVertex->SetPrimary(particle);
  }
  event->AddPrimaryVertex(theVertex);

  //clear before next call
  nParticles = 0;
  theParticles.clear();
  fKinEnergy.clear();
  fDirection.clear();
  fWeight = 1.;
  fTime = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertexFromSpectrum(G4Event* event)
{
  //G4cout << "PrimaryGenerator::GeneratePrimaryVertexFromSpectrum(...)" << G4endl;
  //G4double time1 = 0.*s;
  G4PrimaryVertex *theVertex = new G4PrimaryVertex(decayPos,fTime);
  G4double xx = G4UniformRand();
  size_t bin = 0;
  for(bin=0;bin<theSpectrumI.size();bin++) {
    if(xx<theSpectrumI.at(bin)) break;
  }

  G4double E1 = theSpectrumE.at(bin);
  G4double E2 = theSpectrumE.at(bin+1);

  G4double Energy = G4UniformRand()*(E2-E1) + E1;

  //direction
  G4double pz = 2.*G4UniformRand() - 1.;
  G4double phi = 360.*G4UniformRand()*degree;
  G4double px = sqrt(1.-pz*pz)*std::cos(phi);
  G4double py = sqrt(1.-pz*pz)*std::sin(phi);

  //G4PrimaryParticle* particle = new G4PrimaryParticle(G4Gamma::Gamma());
  G4PrimaryParticle* particle = new G4PrimaryParticle(theSpectrumParticle);

  particle->SetMomentumDirection(G4ThreeVector(px,py,pz));
  particle->SetKineticEnergy(Energy);
  
  theVertex->SetPrimary(particle);
  
  event->AddPrimaryVertex(theVertex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGenerator::ReadSpectrum(G4String FileName)
{
    G4AutoLock lock(&PrimaryGeneratorMutex);
    if(!FileName.size()) {
      G4cout << "PrimaryGenerator::ReadSpectrum() no file name provided" << G4endl;
      return;
    }

    std::ifstream textfile(FileName);
    if(!textfile.is_open()) {
      G4ExceptionDescription ed;
      ed << "      Could not open the file:" << FileName << G4endl;
      G4Exception("PrimaryGenerator::ReadSpectrum(...)",
                  "PrimaryGenerator.01",
                  JustWarning,
                  ed);
    }

    theSpectrumI.clear();
    theSpectrumE.clear();

    std::string particle;
    getline(textfile,particle);
    if(particle=="gamma") {
      theSpectrumParticle = G4Gamma::Gamma();
    } else if (particle=="neutron") {
      theSpectrumParticle = G4Neutron::Neutron();
    } else {
      theSpectrumParticle = G4Geantino::Geantino();
    }

    G4double sumProb = 0.;
    std::string str;
    while(getline(textfile,str)) {

//    G4cout << str << G4endl;
      if(str[0]!='0' && str[0]!='1' && str[0]!='2' && str[0]!='3' && str[0]!='4' && str[0]!='5' && str[0]!='6' && str[0]!='7' && str[0]!='8' && str[0]!='9') continue;
      

      std::stringstream ss(str);
      G4double Emin, Prob;
      ss >> Emin >> Prob;

      theSpectrumI.push_back(Prob);
      theSpectrumE.push_back(Emin);
      sumProb += Prob;

//      G4cout << Emin << "  " << Prob << G4endl;
    }
    textfile.close();

    if(theSpectrumI.back()!=0) {
      G4cout << "The last entry in the spectrum must have prob 0." << G4endl;
      theSpectrumI.clear();
      theSpectrumE.clear();
      return;
    }

    for(auto it=theSpectrumI.begin();it!=theSpectrumI.end();it++) {
      *it = *it/sumProb;
    }

    G4double sum = 0.;
    for(size_t i=0; i<theSpectrumI.size();i++) {
      sum += theSpectrumI.at(i);
      theSpectrumI.at(i) = sum;
    }
    lock.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGenerator::AddParticle(G4ParticleDefinition* aParticle, G4double kinEnergy, G4ThreeVector direction)
{
  theParticles.emplace_back(aParticle);
  fKinEnergy.emplace_back(kinEnergy);
  fDirection.emplace_back(direction);
  ++nParticles;
}