#ifdef FISSION_NEW
#include "MyFissionLibrary.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"
#include "G4ParticleHPManager.hh"
#include "G4NeutronHPDataUsed.hh"

MyFissionLibrary::MyFissionLibrary()
  : G4ParticleHPFinalState(), theIsotope(0), theZ(0), theA(0), targetMass(0.0)
{
  hasXsec = false;
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
  theFissionGenerator = FissionGenerator::Instance();
}

MyFissionLibrary::~MyFissionLibrary()
{}

G4ParticleHPFinalState * MyFissionLibrary::New()
{
  MyFissionLibrary * theNew = new MyFissionLibrary;
  return theNew;
}

void MyFissionLibrary::Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String &, G4ParticleDefinition*)
{

  G4String tString = "/FS/";
  G4bool dbool;
  theZ = static_cast<G4int> (Z);
  theA = static_cast<G4int> (A);
  theIsotope = static_cast<G4int>(1000*Z+A);
  G4ParticleHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), M, dirName, tString, dbool);
  G4String filename = aFile.GetName();

  if(!dbool)
  {
    hasAnyData = false;
    hasFSData = false;
    hasXsec = false;
    return;
  }
  //std::ifstream theData(filename, std::ios::in);
  std::istringstream theData(std::ios::in);
  G4ParticleHPManager::GetInstance()->GetDataStream(filename,theData);

  // here it comes
  G4int infoType, dataType;
  hasFSData = false;
  while (theData >> infoType) // Loop checking, 22.03.2015. T. Koi
  {
    hasFSData = true;
    theData >> dataType;
    switch(infoType)
    {
      case 1:
        if(dataType==4) theNeutronAngularDis.Init(theData);
        if(dataType==5) thePromptNeutronEnDis.Init(theData);
        if(dataType==12) theFinalStatePhotons.InitMean(theData);
        if(dataType==14) theFinalStatePhotons.InitAngular(theData);
        if(dataType==15) theFinalStatePhotons.InitEnergies(theData);
        break;
      case 2:
        if(dataType==1) theFinalStateNeutrons.InitMean(theData);
        break;
      case 3:
        if(dataType==1) theFinalStateNeutrons.InitDelayed(theData);
        if(dataType==5) theDelayedNeutronEnDis.Init(theData);
        break;
      case 4:
        if(dataType==1) theFinalStateNeutrons.InitPrompt(theData);
        break;
      case 5:
        if(dataType==1) theEnergyRelease.Init(theData);
        break;
      default:
        G4cout << "MyFissionLibrary::Init: unknown data type"<<dataType<<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "MyFissionLibrary::Init: unknown data type");
        break;
    }
  }
  targetMass = theFinalStateNeutrons.GetTargetMass();
  //theData.close();
}

G4HadFinalState * MyFissionLibrary::ApplyYourself(const G4HadProjectile & theTrack)
{  
  if ( theResult.Get() == NULL ) theResult.Put( new G4HadFinalState );
  // theResult.Clear();     // geant4.10.00
  theResult.Get()->Clear(); // geant4.10.01

  // prepare neutron
  G4double eKinetic = theTrack.GetKineticEnergy();
  const G4HadProjectile* incidentParticle = &theTrack;
  G4ReactionProduct theNeutron( const_cast<G4ParticleDefinition *>(incidentParticle->GetDefinition()) );
  theNeutron.SetMomentum(incidentParticle->Get4Momentum().vect() );
  theNeutron.SetKineticEnergy(eKinetic);

  // prepare target
  G4Nucleus aNucleus;
  G4ReactionProduct theTarget; 
  G4ThreeVector neuVelo = (1./incidentParticle->GetDefinition()->GetPDGMass())*theNeutron.GetMomentum();
  theTarget = aNucleus.GetBiasedThermalNucleus( targetMass, neuVelo, theTrack.GetMaterial()->GetTemperature());

  // set neutron and target in the FS classes 
  // theNeutronAngularDis.SetNeutron(theNeutron);
  theNeutronAngularDis.SetProjectileRP(theNeutron);
  theNeutronAngularDis.SetTarget(theTarget);

  // boost to target rest system
  theNeutron.Lorentz(theNeutron, -1*theTarget);

  eKinetic = theNeutron.GetKineticEnergy();    

  // dice neutron and gamma multiplicities, energies and momenta in Lab. @@
  // no energy conservation on an event-to-event basis. we rely on the data to be ok. @@
  // also for mean, we rely on the consistency of the data. @@

  G4int nPrompt=0, gPrompt=0;
  //SampleMult(theTrack, &nPrompt, &gPrompt, eKinetic);
  //******************************************************
  G4double promptNeutronMulti = 0;
  promptNeutronMulti = theFinalStateNeutrons.GetPrompt(eKinetic); // prompt nubar from Geant
  G4double delayedNeutronMulti = 0;
  delayedNeutronMulti = theFinalStateNeutrons.GetDelayed(eKinetic); // delayed nubar from Geant
  G4double time = theTrack.GetGlobalTime()/second;
  G4double totalNeutronMulti = theFinalStateNeutrons.GetMean(eKinetic);

  if(promptNeutronMulti==0) {
     // no data for prompt and delayed neutrons in Geant
     // but there is perhaps data for the total neutron multiplicity, in which case 
     // we use it for prompt neutron emission
    promptNeutronMulti = totalNeutronMulti - delayedNeutronMulti;
    // if delayed is 0, prompt will be equal to total multiplicity
  }

  //G4cout << "MyFissionLibrary::SampleMult locking mutex" << G4endl;
  //G4cout << "nCalls = " << nCalls << G4endl;

  //G4cout << "MyFissionLibrary::SampleMult locking mutex" << G4endl;
  //G4AutoLock lk(&MyFissionLibrary::MyFissionLibraryMutex);
  
  //G4cout << "(Z,A) = (" << theZ
  //         << "," << theA << ")" << G4endl;
  //G4cout << "ZAID = " << theIsotope << G4endl;
  //G4cout << "time = " << time << G4endl;
  //G4cout << "totalNeutronMulti = " << totalNeutronMulti << G4endl;
  //G4cout << "eKinetic = " << eKinetic << G4endl;

  //if(aFission) G4cout << "error: aFission is allready pointing somewhere: " << aFission << G4endl;
  //aFission = new fissionEvent(theIsotope, time, totalNeutronMulti, eKinetic, 1);
  fissionEvent *aFission = theFissionGenerator->newFissionEvent(theIsotope, time, totalNeutronMulti, eKinetic, 1);
  //fissionEvent aFission(theIsotope, time, totalNeutronMulti, eKinetic, 1);
  //G4cout << "corr opt = " << aFission->getCorrelationOption() << G4endl;

  nPrompt = aFission->getNeutronNu();
  if(nPrompt == -1) nPrompt = 0; // the fission library libFission.a has no data for neutrons
  gPrompt = aFission->getPhotonNu();
  if (gPrompt == -1) gPrompt = 0; // the fission library libFission.a has no data for gammas

  //******************************************************

  // Build neutrons and add them to dynamic particle vector
  G4double momentum;
  for(G4int i=0; i<nPrompt; i++)
  {
    G4DynamicParticle * it = new G4DynamicParticle;
    it->SetDefinition(G4Neutron::Neutron());
    it->SetKineticEnergy(aFission->getNeutronEnergy(i)*MeV);
    momentum = it->GetTotalMomentum();
    G4ThreeVector temp(momentum*aFission->getNeutronDircosu(i), 
                       momentum*aFission->getNeutronDircosv(i), 
                       momentum*aFission->getNeutronDircosw(i));
    it->SetMomentum( temp );
//    it->SetGlobalTime(aFission->getNeutronAge(i)*second);
    // theResult.AddSecondary(it);     // geant4.10.00
    theResult.Get()->AddSecondary(it); // geant4.10.01
//    G4cout <<"MyFissionLibrary::ApplyYourself: energy of prompt neutron " << i << " = " << it->GetKineticEnergy()<<G4endl;
  }

  // Build gammas, lorentz transform them, and add them to dynamic particle vector
  for(G4int i=0; i<gPrompt; i++)
  {
    G4ReactionProduct thePhoton;// = new G4ReactionProduct;
    thePhoton.SetDefinition(G4Gamma::Gamma());
    thePhoton.SetKineticEnergy(aFission->getPhotonEnergy(i)*MeV);
    momentum = thePhoton.GetTotalMomentum();
    G4ThreeVector temp(momentum*aFission->getPhotonDircosu(i), 
                       momentum*aFission->getPhotonDircosv(i), 
                       momentum*aFission->getPhotonDircosw(i));
    thePhoton.SetMomentum( temp );
    thePhoton.Lorentz(thePhoton, -1.*theTarget);
    
    G4DynamicParticle * it = new G4DynamicParticle;
    it->SetDefinition(thePhoton.GetDefinition());
    it->SetMomentum(thePhoton.GetMomentum());
    // theResult.AddSecondary(it);     // geant4.10.00
    theResult.Get()->AddSecondary(it); // geant4.10.01 
  }

  delete aFission;

//  G4cout <<"MyFissionLibrary::ApplyYourself: Number of induced prompt neutron = "<<nPrompt<<G4endl;                         // geant4.10.00
//  G4cout <<"MyFissionLibrary::ApplyYourself: Number of secondaries = "<<theResult.Get()->GetNumberOfSecondaries()<< G4endl; // geant4.10.01
//  G4cout <<"MyFissionLibrary::ApplyYourself: Number of induced prompt photons = "<<gPrompt<<G4endl;

// finally deal with local energy depositions.
  G4double eDepByFragments = theEnergyRelease.GetFragmentKinetic();
  // theResult.SetLocalEnergyDeposit(eDepByFragments);     // geant4.10.00
  theResult.Get()->SetLocalEnergyDeposit(eDepByFragments); // geant4.10.01
//   G4cout << "MyFissionLibrary::local energy deposit" << eDepByFragments<<G4endl;
  // clean up the primary neutron
  // theResult.SetStatusChange(stopAndKill);     // geant4.10.00
  theResult.Get()->SetStatusChange(stopAndKill); // geant4.10.01
  // return &theResult;   // geant4.10.00
  return theResult.Get(); // geant4.10.01
}
#endif
