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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4DynamicParticle.hh"
//#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4SPSPosDistribution.hh"
#include "PrimaryGenerator.hh"
#include "G4GenericMessenger.hh"

#include "Randomize.hh"

#include "MyRadioactiveDecayBase.hh"
#include "ParticleSourceFromFile.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String aFile)
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{

  fSimpleGun = new G4GeneralParticleSource();
  fUseSimpleGun = false;

  if(aFile=="SimpleGun" || aFile=="simpleGun"  || aFile=="simplegun") {
    fUseSimpleGun = true;
    aFile = "input/activities.txt";
  }

  if(!aFile.size()) {
    //fUseSimpleGun = true;
    aFile = "input/activities.txt";
  }

  fParticleGun  = new PrimaryGenerator();

  fRadDecay = new MyRadioactiveDecayBase();

  fMessenger = new G4GenericMessenger(this,"/PrimaryGeneratorAction/", "...doc...");

  //auto& RestrictDecayCmd =
  fMessenger->DeclareMethod("ExcludeAphaAndSF",
                            &PrimaryGeneratorAction::ExcludeAphaAndSF,
                            "Restric decay to exclude all alpha- and SF- decay channels. Cannot be undone!");
  fMessenger->DeclareMethod("RestrictDecayToSF",
                            &PrimaryGeneratorAction::RestrictDecayToSF,
                            "Restric decay to spontaneous fission only. Cannot be undone!");
  fMessenger->DeclareMethod("RestrictDecayToAlphaN",
                            &PrimaryGeneratorAction::RestrictDecayToAlpha,
                            "Restric decay to alpha only. Cannot be undone!");
  fMessenger->DeclareProperty("UseSimpleGun",
                            fUseSimpleGun,
                            "(de-)activate the use of the G4GeneralParticleSource");
  fMessenger->DeclareProperty("SetSpectrumFile",
                            fSpectrumFileName,
                            "set the name of the file to read the spectrum from");
  fMessenger->DeclareMethod("UseSpectrumFile",
                            &PrimaryGeneratorAction::UseSpectrumFile,
                            "Use the spectrum in file to emit particles");
  fMessenger->DeclareMethod("ShootFromFile",
                            &PrimaryGeneratorAction::ShootFromFile,
                            "set the name of the file to read from");
  fMessenger->DeclareMethod("SetSplitFactor",
                            &PrimaryGeneratorAction::SetSplitFactor,
                            "set the splitting factor for the source from file. Must be integer value.");
  fParticleSourceFromFile = nullptr;

  //fActivityTable = new ActivityTable(aFile, fRadDecay);
  fActivityTable = nullptr;
  fFile = aFile;
  restrictedDecay = false;

  //fActivityTable->RestrictTo("SF decay",fRadDecay);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fMessenger;
  delete fParticleGun;
  delete fSimpleGun;
  delete fParticleSourceFromFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::ShootFromFile(G4String name)
{
  if(fParticleSourceFromFile) delete fParticleSourceFromFile;
  fParticleSourceFromFile = new ParticleSourceFromFile(name);
  fShootFromFile = true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::SetSplitFactor(G4String sVal)
{
  G4int factor = atoi(sVal.data());
  if(factor<1 || factor>10) {
    G4cout << "SetSplitFactor(" << sVal << ") ignored:" << G4endl;
    G4cout << "value outside permissible range [1,10]" << G4endl;
  }
  if(fParticleSourceFromFile) {
    fParticleSourceFromFile->SetSplitFactor(factor);
  } else {
    G4cout << "SetSplitFactor(" << sVal << ") ignored:" << G4endl;
    G4cout << "ShootFromFile not yet initilised" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if(fShootFromFile) {
    fParticleSourceFromFile->GeneratePrimaryVertex(anEvent);
    return;
  }

  if(!fActivityTable) fActivityTable = new ActivityTable(fFile, fRadDecay);
  //fTimer.Start();
  if(fUseSimpleGun) {
    fSimpleGun->GeneratePrimaryVertex(anEvent);
    //fTimer.Stop();
    return;
  }

  if(fUseSpectrum) {
    fParticleGun->GenerateDecayPos();
    fParticleGun->GeneratePrimaryVertexFromSpectrum(anEvent);
    return;
  }

  if(!fActivityTable->IsInit()) {
    fParticleGun->GeneratePrimaryVertex(anEvent);
    G4cout << "Error in the activity file shooting geantino!!" << G4endl;
    //fTimer.Stop();
    return;
  }

  //G4cout << "PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)" << G4endl;
  //const G4double rate = fActivityTable->GetTotalActivity()/nWorkers;
  //theTime += G4RandExponential::shoot(1./rate);

  //if(!nPrimaries) nPrimaries = fParticleGun->GetNumberOfParticles();

  G4int Z = 0;
  G4int A = 0;
  G4double excitationEnergy = 0.;
  bool applyBRbias = false;

  fActivityTable->GenerateMotherNuclide(Z,A,excitationEnergy,applyBRbias);
  G4int decayChannelNbr = -1;
  if(restrictedDecay) decayChannelNbr = fActivityTable->GetRestrictedDecayChannelNbr();
  //excitationEnergy = 661.659*keV;

  //G4cout << "G4IonTable::GetIonTable()->GetNucleusEncoding(Z,A,0) = " << G4IonTable::GetIonTable()->GetNucleusEncoding(Z,A,0) << G4endl;

//---------------------------------------------------
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitationEnergy);

  theDecayTable = fRadDecay->GetDecayTable(ion);
  fParticleGun->GenerateDecayPos();

  if (theDecayTable == 0 || theDecayTable->entries() == 0) {
    // No data in the decay table.
    G4cout << "PrimaryGeneratorAction::GeneratePrimaries : "
           << "decay table not defined for "
           << ion->GetParticleName() 
           << G4endl;

    fParticleGun->GeneratePrimaryVertex(anEvent);
    G4cout << "Error in the activity file shooting geantino!!" << G4endl;
    return;

  }

  //const G4double timeCutOff = 1.E30*ns; //this is in ns
  //const G4double timeCutOff = 1.E3*ns; //this is in ns
  //const G4double timeCutOff = DBL_MAX; //this is in ns
  const G4double timeCutOff = 1.*ns; //only take the prompt emission into account
  G4double time = 0.;
  G4bool hasDaughter = true;
  G4int loopCount = 0;
  while(theDecayTable && hasDaughter && time<timeCutOff) {
    if(!theDecayTable->entries()) break;
    if(loopCount>10) break; //do not continue decay chain neyond 10 decays
    // Data found.  Try to decay nucleus
    // theDecayTable->DumpInfo();
    // Decay without variance reduction
    G4DecayProducts* products = 0;
    G4double weight = 1.;

    products = DoDecay(*ion,decayChannelNbr);
    if(!products) break;
    decayChannelNbr = -1; //only applies to the first in the chain

    fParticleGun->SetWeight(weight);
    fParticleGun->SetTime(time);

    hasDaughter = false; //the fission process does not produce any daugthers
    G4int numberOfSecondaries = products->entries();
    for (G4int index = 0; index < numberOfSecondaries; index++) {
      G4DynamicParticle *secondary = (*products)[index];
      G4ParticleDefinition* particleDef = (G4ParticleDefinition*) secondary->GetParticleDefinition();
      G4double kinEnergy = secondary->GetKineticEnergy();

      if(particleDef->GetPDGCharge()>2) { //this is the daugther of the decay => should be decayed further
        ion = particleDef;
        theDecayTable = fRadDecay->GetDecayTable(ion);
        hasDaughter = true;
        G4double LifeTime = ion->GetPDGLifeTime();
        //f(LifeTime>0.) time += G4RandExponential::shoot(LifeTime);
        if(LifeTime>0.) time += -LifeTime*log(G4UniformRand());
        //G4cout << "LifeTime = " << LifeTime << G4endl;
        //G4cout << "time = " << time << G4endl;
        continue; // do not add it to the primary vertex
      }

      /*if(particleDef->GetPDGEncoding()==11) {
        continue; //don't produce any electrons from B- decays
      }*/

      if(kinEnergy>=20.*MeV) {
        G4cout << "Warning: a " << particleDef->GetParticleName()
               << "with kinetic energy >20 MeV has been generated, "
               << "transport will lead to a fatal exception. "
               << "It will therefore be rejected!" << G4endl;
        continue;
      }
      #ifdef DEBUG_ALPHA_N
        //this is for setting up the (a,n) biasing, only alpha particle will be shot
        //if there is no alpha particle to shoot a geantino will be shot in PrimaryGenerator.cc
        if(particleDef->GetPDGEncoding()==1000020040) fParticleGun->AddParticle(particleDef,kinEnergy,secondary->GetMomentumDirection());
      #else
        fParticleGun->AddParticle(particleDef,kinEnergy,secondary->GetMomentumDirection());
      #endif
    }

    fParticleGun->GeneratePrimaryVertex(anEvent); ++loopCount;

    delete products;
    products = nullptr;
  }

  if(!loopCount) { //this should not happen
    G4cout << "Warning no decay products generated shooting a geantino" << G4endl;
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
    //fTimer.Stop();
    //G4cout << "Time spent in PrimaryGeneratorAction " << fTimer.GetRealElapsed() << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4DecayProducts*
PrimaryGeneratorAction::DoDecay(const G4ParticleDefinition& theParticleDef)
{
  G4DecayProducts* products = 0;
  // Choose a decay channel.
  //G4cout << "Select a channel..." << G4endl;

  // G4DecayTable::SelectADecayChannel checks to see if sum of daughter masses
  // exceeds parent mass. Pass it the parent mass + maximum Q value to account
  // for difference in mass defect.
  G4double parentPlusQ = theParticleDef.GetPDGMass() + 30.*MeV;
  G4VDecayChannel* theDecayChannel = theDecayTable->SelectADecayChannel(parentPlusQ);
  //theRadDecayMode = (static_cast<G4NuclearDecay*>(theDecayChannel))->GetDecayMode();

  //G4cout << "DoDecay. theDecayChannel->GetKinematicsName() = " << theDecayChannel->GetKinematicsName() << G4endl;

  if (theDecayChannel == 0) {
    // Decay channel not found.
    G4ExceptionDescription ed;
    ed << " Cannot determine decay channel for " << theParticleDef.GetParticleName() << G4endl;
    G4Exception("G4RadioactiveDecay::DoDecay", "HAD_RDM_013",
                FatalException, ed);
  } else {
    // A decay channel has been identified, so execute the DecayIt.
    //G4cout << "G4RadioactiveDecay::DoIt : selected decay channel addr: "
    //       << theDecayChannel << G4endl;

    products = theDecayChannel->DecayIt(theParticleDef.GetPDGMass() );
  }

  return products;
}

G4DecayProducts*
PrimaryGeneratorAction::DoDecay(const G4ParticleDefinition& theParticleDef, G4int channelNbr)
{

  G4DecayProducts* products = 0;
  if(channelNbr==-1) {
    products = DoDecay(theParticleDef); //decay it normally by choosing radnom channel
    return products;
  }
  // Choose a decay channel.
  //if (GetVerboseLevel() > 1) theDecayTable ->DumpInfo();

  G4VDecayChannel* theDecayChannel = theDecayTable->GetDecayChannel(channelNbr);

  if (theDecayChannel == 0) {
    // Decay channel not found.

      G4cout << " PrimaryGeneratorAction::DoDecayBRbias : cannot determine decay channel ";
      G4cout << " for this nucleus; decay as if no biasing active. ";
      G4cout << G4endl;
      theDecayTable ->DumpInfo();

    products = DoDecay(theParticleDef);  // DHW 6 Dec 2010 - do decay as if no biasing
                                            //           to avoid deref of temppprods = 0
  } else {
    // A decay channel has been identified, so execute the DecayIt.
    G4double tempmass = theParticleDef.GetPDGMass();
    products = theDecayChannel->DecayIt(tempmass);
    //weigth *= theDecayChannel->GetBR();
  }

  return products;
}

G4DecayProducts*
PrimaryGeneratorAction::DoDecayBRbias(const G4ParticleDefinition& theParticleDef, G4double &weigth)
{
  G4DecayProducts* products = 0;
  // Choose a decay channel.
  //if (GetVerboseLevel() > 1) theDecayTable ->DumpInfo();
  
  G4int ndecaych = G4int(theDecayTable->entries()*G4UniformRand());
  G4VDecayChannel* theDecayChannel = theDecayTable->GetDecayChannel(ndecaych);

  if (theDecayChannel == 0) {
    // Decay channel not found.

      G4cout << " PrimaryGeneratorAction::DoDecayBRbias : cannot determine decay channel ";
      G4cout << " for this nucleus; decay as if no biasing active. ";
      G4cout << G4endl;
      theDecayTable ->DumpInfo();

    products = DoDecay(theParticleDef);  // DHW 6 Dec 2010 - do decay as if no biasing
                                            //           to avoid deref of temppprods = 0
  } else {
    // A decay channel has been identified, so execute the DecayIt.
    G4double tempmass = theParticleDef.GetPDGMass();
    products = theDecayChannel->DecayIt(tempmass);
    weigth *= (theDecayChannel->GetBR())*(theDecayTable->entries());
  }

  return products;
}
