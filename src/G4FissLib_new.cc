#ifdef FISSION_NEW
#include "G4FissLib_new.hh"
//#include "G4FissionLibrary_new.hh"
#include "MyFissionLibrary.hh"
#include "G4SystemOfUnits.hh"

#include "G4NeutronHPManager.hh"

G4FissLib_new::G4FissLib_new()
:G4HadronicInteraction("FissLib_new")
{
  xSec = 0;
  SetMinEnergy(0.0);
  SetMaxEnergy(20.*MeV);
  if(!getenv("G4NEUTRONHPDATA")) {
     G4cout << "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files." << G4endl;
     throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
  }
  dirName = getenv("G4NEUTRONHPDATA");
  G4String tString = "/Fission/";
  dirName = dirName + tString;

  numEle = G4Element::GetNumberOfElements();
  //G4FissionLibrary_new *theFS = new G4FissionLibrary_new;
  MyFissionLibrary *theFS = new MyFissionLibrary;
  for (G4int i=0; i<numEle; i++)
  {
    theFission.push_back( new G4NeutronHPChannel );
    (*theFission[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
    (*theFission[i]).Register(theFS);
  }
  delete theFS;
  
}
  
G4FissLib_new::~G4FissLib_new()
{
  G4cout << "G4FissLib_new::~G4FissLib_new() 1" << G4endl;
  //delete [] theFission;

  //G4cout << "theFission.size() = " << theFission.size() << G4endl;
  for ( std::vector<G4NeutronHPChannel*>::iterator 
        ite = theFission.begin() ; ite != theFission.end() ; ite++ )
  {
    //G4cout << "ii = " << ii++ << G4endl;
    //G4cout << "ite = " << (*ite) << G4endl;
    //delete *ite;
    //G4cout << "ite = " << (*ite) << G4endl;
  }
  theFission.clear();
 // G4cout << "theFission.size() = " << theFission.size() << G4endl;
  G4cout << "G4FissLib_new::~G4FissLib_new() 2" << G4endl;
}
  
G4HadFinalState*
G4FissLib_new::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus)
{
  if ( numEle < (G4int)G4Element::GetNumberOfElements() ) addChannelForNewElement();
  G4NeutronHPManager::GetInstance()->OpenReactionWhiteBoard();

  const G4Material* theMaterial = aTrack.GetMaterial();
  G4int n = theMaterial->GetNumberOfElements();
  G4int index = theMaterial->GetElement(0)->GetIndex();

  if (n != 1) {
    xSec = new G4double[n];
    G4double sum = 0;
    G4int i;
    //G4int imat;
    const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
    G4double rWeight;    
    G4NeutronHPThermalBoost aThermalE;
    for (i = 0; i < n; i++) {
      //imat = theMaterial->GetElement(i)->GetIndex();
      index = theMaterial->GetElement(i)->GetIndex();
      rWeight = NumAtomsPerVolume[i];
      //xSec[i] = theFission[imat].GetXsec(aThermalE.GetThermalEnergy(aTrack,theMaterial->GetElement(i),theMaterial->GetTemperature()));
      xSec[i] = (*theFission[index]).GetXsec(aThermalE.GetThermalEnergy(aTrack,theMaterial->GetElement(i),theMaterial->GetTemperature()));
      xSec[i] *= rWeight;
      sum+=xSec[i];
    }

    G4double random = G4UniformRand();
    G4double running = 0;
    for (i = 0; i < n; i++) {
      running += xSec[i];
      index = theMaterial->GetElement(i)->GetIndex();
      //if(random<=running/sum) break;
      if( sum == 0 || random <= running/sum ) break;
    }
    if(i==n) i=std::max(0, n-1);
    delete [] xSec;
  }

  //return theFission[index].ApplyYourself(aTrack);
  G4HadFinalState* result = (*theFission[index]).ApplyYourself(aTrack);
  aNucleus.SetParameters(G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
  G4NeutronHPManager::GetInstance()->CloseReactionWhiteBoard();
  return result;
}

const std::pair<G4double, G4double> G4FissLib_new::GetFatalEnergyCheckLevels() const
{
  // max energy non-conservation is mass of heavy nucleus (taken from G4LFission)
  return std::pair<G4double, G4double>(5*perCent,250*GeV);
}

void G4FissLib_new::addChannelForNewElement()
{
  //G4FissionLibrary_new* theFS = new G4FissionLibrary_new;
  MyFissionLibrary* theFS = new MyFissionLibrary;
  for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ ) 
  {
     G4cout << "G4FissionLibrary_new Prepairing Data for the new element of " << (*(G4Element::GetElementTable()))[i]->GetName() << G4endl;
    theFission.push_back( new G4NeutronHPChannel );
    (*theFission[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
    (*theFission[i]).Register(theFS);
  }
  delete theFS;
  numEle = (G4int)G4Element::GetNumberOfElements();
}
#endif
