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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   MyRadioactiveDecayBase.cc                                         //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   9 August 2017                                                     //
//  Description: version the G4RadioactiveDecay process by F. Lei and         //
//               P.R. Truscott with biasing and activation calculations       //
//               removed to a derived class.  It performs alpha, beta,        //
//               electron capture and isomeric transition decays of           //
//               radioactive nuclei.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "MyRadioactiveDecayBase.hh"
#include "MyRadioactiveDecayBaseMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4ParticleChangeForRadDecay.hh"
#include "G4ITDecay.hh"
#include "G4BetaDecayType.hh"
#include "G4BetaMinusDecay.hh"
#include "G4BetaPlusDecay.hh"
#include "G4ECDecay.hh"
#include "G4AlphaDecay.hh"
#include "G4TritonDecay.hh"
#include "G4ProtonDecay.hh"
#include "G4NeutronDecay.hh"
#include "MySFDecay.hh"
#include "G4VDecayChannel.hh"
#include "G4NuclearDecay.hh"
#include "G4RadioactiveDecayMode.hh"
#include "G4Fragment.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4BetaDecayType.hh"
#include "Randomize.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4LevelManager.hh"
#include "G4ThreeVector.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"
#include "G4Triton.hh"
#include "G4Proton.hh"

#include "G4HadronicProcessType.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicException.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4PhotonEvaporation.hh"

#include <vector>
#include <sstream>
#include <algorithm>
#include <fstream>

#include <filesystem>

using namespace CLHEP;

const G4double MyRadioactiveDecayBase::levelTolerance = 10.0*eV;
const G4ThreeVector MyRadioactiveDecayBase::origin(0.,0.,0.);

#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
G4Mutex MyRadioactiveDecayBase::radioactiveDecayMutex = G4MUTEX_INITIALIZER;
DecayTableMap* MyRadioactiveDecayBase::master_dkmap = 0;

G4int& MyRadioactiveDecayBase::NumberOfInstances()
{
  static G4int numberOfInstances = 0;
  return numberOfInstances;
}
#endif

MyRadioactiveDecayBase::MyRadioactiveDecayBase(const G4String& processName)
 : G4VRestDiscreteProcess(processName, fDecay), isInitialised(false),
   forceDecayDirection(0.,0.,0.), forceDecayHalfAngle(0.*deg), dirPath(""),
   verboseLevel(0)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "MyRadioactiveDecayBase constructor: processName = " << processName
           << G4endl;
  }
#endif

  SetProcessSubType(fRadioactiveDecay);

  theRadioactiveDecayBaseMessenger = new MyRadioactiveDecayBaseMessenger(this);
  pParticleChange = &fParticleChangeForRadDecay;

  // Set up photon evaporation for use in G4ITDecay
  photonEvaporation = new G4PhotonEvaporation();
  photonEvaporation->RDMForced(true);
  photonEvaporation->SetICM(true);

  // DHW G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  // DHW deex->SetCorrelatedGamma(true);

  // Check data directory
  char* path_var = std::getenv("G4RADIOACTIVEDATA");
  if (!path_var) {
    G4Exception("G4RadioactiveDecay()","HAD_RDM_200",FatalException,
                "Environment variable G4RADIOACTIVEDATA is not set");
  } else {
    dirPath = path_var;   // convert to string
    std::ostringstream os;
    os << dirPath << "/z1.a3";   // used as a dummy 
    std::ifstream testFile;
    testFile.open(os.str() );
    if (!testFile.is_open() )
      G4Exception("G4RadioactiveDecay()","HAD_RDM_201",FatalException,
                  "Environment variable G4RADIOACTIVEDATA is set, but does not point to correct directory");
  }

  // Reset the list of user defined data files
  theUserRadioactiveDataFiles.clear();

  // Instantiate the map of decay tables
#ifdef G4MULTITHREADED
  G4AutoLock lk(&MyRadioactiveDecayBase::radioactiveDecayMutex);
  //G4cout << "MyRadioactiveDecayBase:: locking mutex 155" << G4endl;
  NumberOfInstances()++;
  if(!master_dkmap) master_dkmap = new DecayTableMap;
#endif
  dkmap = new DecayTableMap;

  // Apply default values
  applyARM = true;
  applyICM = true;  // Always on; keep only for backward compatibility
  applyBRbias = false; //for now this will be false can be changed from UI command
 
  // RDM applies to all logical volumes by default
  isAllVolumesMode = true;
  SelectAllVolumes();
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
  #ifdef G4MULTITHREADED
  //G4cout << "MyRadioactiveDecayBase:: unlocking mutex 172" << G4endl;
  lk.unlock();
#endif
}

void MyRadioactiveDecayBase::ProcessDescription(std::ostream& outFile) const
{
  outFile << "The radioactive decay process (G4RadioactiveDecay) handles the\n"
          << "alpha, beta+, beta-, electron capture and isomeric transition\n"
          << "decays of nuclei (G4GenericIon) with masses A > 4.\n"
          << "The required half-lives and decay schemes are retrieved from\n"
          << "the RadioactiveDecay database which was derived from ENSDF.\n";
}


MyRadioactiveDecayBase::~MyRadioactiveDecayBase()
{
  delete theRadioactiveDecayBaseMessenger;
  delete photonEvaporation;
  for (DecayTableMap::iterator i = dkmap->begin(); i != dkmap->end(); i++) {
    delete i->second;
  }
  dkmap->clear();
  delete dkmap;
#ifdef G4MULTITHREADED
  G4AutoLock lk(&MyRadioactiveDecayBase::radioactiveDecayMutex);
  G4cout << "MyRadioactiveDecayBase:: locking mutex 197" << G4endl;
  
  --NumberOfInstances();
  if(NumberOfInstances()==0)
  {
    for (DecayTableMap::iterator i = master_dkmap->begin(); i != master_dkmap->end(); i++) {
      delete i->second;
    }
    master_dkmap->clear();
    delete master_dkmap;
  }
  lk.unlock();
  G4cout << "MyRadioactiveDecayBase:: unlocking mutex 209" << G4endl;
#endif
}


G4bool MyRadioactiveDecayBase::IsApplicable(const G4ParticleDefinition& aParticle)
{
  // All particles other than G4Ions, are rejected by default
  if (((const G4Ions*)(&aParticle))->GetExcitationEnergy() > 0.) {return true;}
  if (aParticle.GetParticleName() == "GenericIon") {
    return true;
  } else if (!(aParticle.GetParticleType() == "nucleus")
             || aParticle.GetPDGLifeTime() < 0. ) {
    return false;
  }

  // Determine whether the nuclide falls into the correct A and Z range
  G4int A = ((const G4Ions*) (&aParticle))->GetAtomicMass();
  G4int Z = ((const G4Ions*) (&aParticle))->GetAtomicNumber();

  if (A > theNucleusLimits.GetAMax() || A < theNucleusLimits.GetAMin())
    {return false;}
  else if (Z > theNucleusLimits.GetZMax() || Z < theNucleusLimits.GetZMin())
    {return false;}
  return true;
}

G4DecayTable* MyRadioactiveDecayBase::GetDecayTable(const G4ParticleDefinition* aNucleus)
{
  G4String key = aNucleus->GetParticleName();
  DecayTableMap::iterator table_ptr = dkmap->find(key);

  G4DecayTable* theDecayTable = 0;
  if (table_ptr == dkmap->end() ) {                   // If table not there,     
    theDecayTable = LoadDecayTable(*aNucleus);        // load from file and
    if(theDecayTable) (*dkmap)[key] = theDecayTable;  // store in library 
  } else {
    theDecayTable = table_ptr->second;
  }
  return theDecayTable;
}


void MyRadioactiveDecayBase::SelectAVolume(const G4String aVolume)
{
  G4LogicalVolumeStore* theLogicalVolumes;
  G4LogicalVolume* volume;
  theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  for (size_t i = 0; i < theLogicalVolumes->size(); i++) {
    volume = (*theLogicalVolumes)[i];
    if (volume->GetName() == aVolume) {
      ValidVolumes.push_back(aVolume);
      std::sort(ValidVolumes.begin(), ValidVolumes.end());
      // sort need for performing binary_search

      if (GetVerboseLevel() > 0)
	G4cout << " Radioactive decay applied to " << aVolume << G4endl; 

    } else if (i == theLogicalVolumes->size() ) {
      G4cout << " G4RadioactiveDecay::SelectAVolume: " << aVolume
             << " is not a valid logical volume name."
             << " Decay not activated for it." << G4endl; 
    }
  }
}


void MyRadioactiveDecayBase::DeselectAVolume(const G4String aVolume)
{
  G4LogicalVolumeStore* theLogicalVolumes;
  G4LogicalVolume* volume;
  theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  for (size_t i = 0; i < theLogicalVolumes->size(); i++) {
    volume = (*theLogicalVolumes)[i];
    if (volume->GetName() == aVolume) {
      std::vector<G4String>::iterator location;
      location = std::find(ValidVolumes.begin(),ValidVolumes.end(),aVolume);
      if (location != ValidVolumes.end() ) {
        ValidVolumes.erase(location);
        std::sort(ValidVolumes.begin(), ValidVolumes.end());
        isAllVolumesMode = false;
        if (GetVerboseLevel() > 0)
          G4cout << " G4RadioactiveDecay::DeselectAVolume: " << aVolume
                 << " is removed from list " << G4endl;
      } else {
        G4cout << " G4RadioactiveDecay::DeselectAVolume: " << aVolume 
               << " is not in the list.  No action taken. " << G4endl; 
      }
    } else if (i ==  theLogicalVolumes->size()) {
      G4cout << " G4RadioactiveDecay::DeselectVolume:" << aVolume
             << " is not a valid logical volume name.  No action taken." 
             << G4endl; 
    }
  }
}


void MyRadioactiveDecayBase::SelectAllVolumes() 
{
  G4LogicalVolumeStore* theLogicalVolumes;
  G4LogicalVolume* volume;
  theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  ValidVolumes.clear();
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0)
    G4cout << " RDM Applies to all Volumes"  << G4endl;
#endif
  for (size_t i = 0; i < theLogicalVolumes->size(); i++){
    volume = (*theLogicalVolumes)[i];
    ValidVolumes.push_back(volume->GetName());    
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
      G4cout << "       RDM Applies to Volume " << volume->GetName() << G4endl;
#endif
  }
  std::sort(ValidVolumes.begin(), ValidVolumes.end());
  // sort needed in order to allow binary_search
  isAllVolumesMode=true;
}


void MyRadioactiveDecayBase::DeselectAllVolumes() 
{
  ValidVolumes.clear();
  isAllVolumesMode=false;
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) G4cout << "RDM removed from all volumes" << G4endl; 
#endif
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  GetMeanLifeTime (required by the base class)                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4double MyRadioactiveDecayBase::GetMeanLifeTime(const G4Track& theTrack,
                                                 G4ForceCondition*)
{
  G4double meanlife = 0.;
  const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
  const G4ParticleDefinition* theParticleDef = theParticle->GetDefinition();
  G4double theLife = theParticleDef->GetPDGLifeTime();
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << "G4RadioactiveDecay::GetMeanLifeTime() " << G4endl;
    G4cout << "KineticEnergy: " << theParticle->GetKineticEnergy()/GeV
           << " GeV, Mass: " << theParticle->GetMass()/GeV
           << " GeV, Life time: " << theLife/ns << " ns " << G4endl;
  }
#endif
  if (theParticleDef->GetPDGStable()) {meanlife = DBL_MAX;}
  else if (theLife < 0.0) {meanlife = DBL_MAX;}
  else {meanlife = theLife;}
  // Set meanlife to zero for excited istopes which are not in the
  // RDM database
  if (((const G4Ions*)(theParticleDef))->GetExcitationEnergy() > 0. &&
                                        meanlife == DBL_MAX) {meanlife = 0.;}
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout << " mean life time: " << meanlife/s << " s " << G4endl;
#endif

  return meanlife;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  GetMeanFreePath for decay in flight                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4double MyRadioactiveDecayBase::GetMeanFreePath (const G4Track& aTrack, G4double,
                                              G4ForceCondition*)
{
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
  G4double tau = aParticleDef->GetPDGLifeTime();
  G4double aMass = aParticle->GetMass();

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << "G4RadioactiveDecay::GetMeanFreePath() " << G4endl;
    G4cout << "  KineticEnergy: " << aParticle->GetKineticEnergy()/GeV
           << " GeV, Mass: " << aMass/GeV << " GeV, tau: " << tau << " ns "
           << G4endl;
  }
#endif
  G4double pathlength = DBL_MAX;
  if (tau != -1) {
    // Ion can decay

    if (tau < -1000.0) {
      pathlength = DBL_MIN;  // nuclide had very short lifetime or wasn't in table

    } else if (tau < 0.0) {
      G4cout << aParticleDef->GetParticleName() << " has lifetime " << tau << G4endl;
      G4ExceptionDescription ed;
      ed << "Ion has negative lifetime " << tau
         << " but is not stable.  Setting mean free path to DBL_MAX" << G4endl; 
      G4Exception("G4RadioactiveDecay::GetMeanFreePath()", "HAD_RDM_011",
                   JustWarning, ed);
      pathlength = DBL_MAX;

    } else {
      // Calculate mean free path
      G4double betaGamma = aParticle->GetTotalMomentum()/aMass;
      pathlength = c_light*tau*betaGamma;

      if (pathlength < DBL_MIN) {
        pathlength = DBL_MIN;
#ifdef G4VERBOSE
        if (GetVerboseLevel() > 2) {
          G4cout << "G4Decay::GetMeanFreePath: "
                 << aParticleDef->GetParticleName()
                 << " stops, kinetic energy = "
                 << aParticle->GetKineticEnergy()/keV <<" keV " << G4endl;
        }
#endif
      }
    }
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "mean free path: "<< pathlength/m << " m" << G4endl;
  }
#endif
  return  pathlength;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  BuildPhysicsTable - initialization of atomic de-excitation                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void MyRadioactiveDecayBase::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if (!isInitialised) {
    isInitialised = true;
#ifdef G4VERBOSE
    if(G4Threading::IsMasterThread()) { StreamInfo(G4cout, "\n"); }
#endif
  }
  G4HadronicProcessStore::
   Instance()->RegisterParticleForExtraProcess(this,G4GenericIon::GenericIon());
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  StreamInfo - stream out parameters                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void
MyRadioactiveDecayBase::StreamInfo(std::ostream& os, const G4String& endline)
{
  G4DeexPrecoParameters* deex =
    G4NuclearLevelData::GetInstance()->GetParameters();
  G4EmParameters* emparam = G4EmParameters::Instance();

  G4int prec = os.precision(5);
  os << "======================================================================"
     << endline;
  os << "======          Radioactive Decay Physics Parameters           ======="
     << endline;
  os << "======================================================================"
     << endline;
  os << "Max life time                                     "
     << deex->GetMaxLifeTime()/CLHEP::ps << " ps" << endline;
  os << "Internal e- conversion flag                       "
     << deex->GetInternalConversionFlag() << endline;
  os << "Stored internal conversion coefficients           "
     << deex->StoreICLevelData() << endline;
  os << "Enable correlated gamma emission                  "
     << deex->CorrelatedGamma() << endline;
  os << "Max 2J for sampling of angular correlations       "
     << deex->GetTwoJMAX() << endline;
  os << "Atomic de-excitation enabled                      "
     << emparam->Fluo() << endline;
  os << "Auger electron emission enabled                   "
     << emparam->Auger() << endline;
  os << "Auger cascade enabled                             "
     << emparam->AugerCascade() << endline;
  os << "Check EM cuts disabled for atomic de-excitation   "
     << emparam->DeexcitationIgnoreCut() << endline;
  os << "Use Bearden atomic level energies                 "
     << emparam->BeardenFluoDir() << endline;
  os << "======================================================================"
     << endline;
  os.precision(prec);
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  LoadDecayTable loads the decay scheme from the RadioactiveDecay database  // 
//  for the parent nucleus.                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4DecayTable*
MyRadioactiveDecayBase::LoadDecayTable(const G4ParticleDefinition& theParentNucleus)
{
  // Generate input data file name using Z and A of the parent nucleus
  // file containing radioactive decay data.
  G4int A = ((const G4Ions*)(&theParentNucleus))->GetAtomicMass();
  G4int Z = ((const G4Ions*)(&theParentNucleus))->GetAtomicNumber();

  G4double levelEnergy = ((const G4Ions*)(&theParentNucleus))->GetExcitationEnergy();
  G4Ions::G4FloatLevelBase floatingLevel =
    ((const G4Ions*)(&theParentNucleus))->GetFloatLevelBase();

#ifdef G4MULTITHREADED
  G4AutoLock lk(&MyRadioactiveDecayBase::radioactiveDecayMutex);
  G4cout << "LoadDecayTable locking its mutex 526" << G4endl;
  G4cout << "theParentNucleus.GetParticleName() = " << theParentNucleus.GetParticleName() << G4endl;

  G4String key = theParentNucleus.GetParticleName();
  DecayTableMap::iterator master_table_ptr = master_dkmap->find(key);

  if (master_table_ptr != master_dkmap->end() ) {   // If table is there             
    return master_table_ptr->second;
  }
#endif

  //Check if data have been provided by the user
  G4String file = theUserRadioactiveDataFiles[1000*A+Z];

  if (file == "") {
    std::ostringstream os;
    os << dirPath << "/z" << Z << ".a" << A << '\0';
    file = os.str();
  }

  G4DecayTable* theDecayTable = new G4DecayTable();
  G4bool found(false);     // True if energy level matches one in table

  std::ifstream DecaySchemeFile;
  DecaySchemeFile.open(file);

  if (DecaySchemeFile.good()) {
    // Initialize variables used for reading in radioactive decay data
    G4bool floatMatch(false);
    const G4int nMode = 12;
    G4double modeTotalBR[nMode] = {0.0};
    G4double modeSumBR[nMode];
    for (G4int i = 0; i < nMode; i++) {
      modeSumBR[i] = 0.0;
    }

    char inputChars[120]={' '};
    G4String inputLine;
    G4String recordType("");
    G4String floatingFlag("");
    G4String daughterFloatFlag("");
    G4Ions::G4FloatLevelBase daughterFloatLevel;
    G4RadioactiveDecayMode theDecayMode;
    G4double decayModeTotal(0.0);
    G4double parentExcitation(0.0);
    G4double a(0.0);
    G4double b(0.0);
    G4double c(0.0);
    G4double dummy(0.0);
    G4BetaDecayType betaType(allowed);

    // Loop through each data file record until you identify the decay
    // data relating to the nuclide of concern.

    G4bool complete(false);  // bool insures only one set of values read for any
                             // given parent energy level
    G4int loop = 0;
    while (!complete && !DecaySchemeFile.getline(inputChars, 120).eof()) {  /* Loop checking, 01.09.2015, D.Wright */
      loop++;
      if (loop > 100000) {
        G4Exception("G4RadioactiveDecay::LoadDecayTable()", "HAD_RDM_100",
                    JustWarning, "While loop count exceeded");
        break;
      }
 
      inputLine = inputChars;
      inputLine = inputLine.strip(1);
      if (inputChars[0] != '#' && inputLine.length() != 0) {
        std::istringstream tmpStream(inputLine);

        if (inputChars[0] == 'P') {
          // Nucleus is a parent type.  Check excitation level to see if it
          // matches that of theParentNucleus
          tmpStream >> recordType >> parentExcitation >> floatingFlag >> dummy;
          // "dummy" takes the place of half-life
          //  Now read in from ENSDFSTATE in particle category

          if (found) {
            complete = true;
          } else {
            // Take first level which matches excitation energy regardless of floating level
            found = (std::abs(parentExcitation*keV - levelEnergy) < levelTolerance);
            if (floatingLevel != noFloat) {
              // If floating level specificed, require match of both energy and floating level
              floatMatch = (floatingLevel == G4Ions::FloatLevelBase(floatingFlag.back()) );
              if (!floatMatch) found = false;
            }
          }

        } else if (found) {
          // The right part of the radioactive decay data file has been found.  Search
          // through it to determine the mode of decay of the subsequent records.

          // Store for later the total decay probability for each decay mode 
          if (inputLine.length() < 72) {
            tmpStream >> theDecayMode >> dummy >> decayModeTotal;
            switch (theDecayMode) {
              case IT:
                {
                G4ITDecay* anITChannel = new G4ITDecay(&theParentNucleus, decayModeTotal,
                                                       0.0, 0.0, photonEvaporation);
//                anITChannel->SetHLThreshold(halflifethreshold);
                anITChannel->SetARM(applyARM);
                theDecayTable->Insert(anITChannel);
//                anITChannel->DumpNuclearInfo();
                }
                break;
              case BetaMinus:
                modeTotalBR[1] = decayModeTotal; break;
              case BetaPlus:
                modeTotalBR[2] = decayModeTotal; break;
              case KshellEC:
                modeTotalBR[3] = decayModeTotal; break;
              case LshellEC:
                modeTotalBR[4] = decayModeTotal; break;
              case MshellEC:
                modeTotalBR[5] = decayModeTotal; break;
              case NshellEC:
                modeTotalBR[6] = decayModeTotal; break;
              case Alpha:
                modeTotalBR[7] = decayModeTotal; break;
              case Proton:
                modeTotalBR[8] = decayModeTotal; break;
              case Neutron:
                modeTotalBR[9] = decayModeTotal; break;
              case BDProton:
                break;
              case BDNeutron:
                break;
              case Beta2Minus:
                break;
              case Beta2Plus:
                break;
              case Proton2:
                break;
              case Neutron2:
                break;
              case SpFission:
                modeTotalBR[10] = decayModeTotal; break;
              case Triton:
                modeTotalBR[11] = decayModeTotal; break;
              case RDM_ERROR:

              default:
                G4Exception("G4RadioactiveDecay::LoadDecayTable()", "HAD_RDM_000",
                            FatalException, "Selected decay mode does not exist");
            }  // switch

          } else {
            if (inputLine.length() < 84) {
              tmpStream >> theDecayMode >> a >> daughterFloatFlag >> b >> c;
              betaType = allowed;
            } else {
              tmpStream >> theDecayMode >> a >> daughterFloatFlag >> b >> c >> betaType;
            }

            // Allowed transitions are the default. Forbidden transitions are
            // indicated in the last column.
            a /= 1000.;
            c /= 1000.;
            b /= 100.;
            daughterFloatLevel = G4Ions::FloatLevelBase(daughterFloatFlag.back());

            switch (theDecayMode) {
              case BetaMinus:
              {
                G4BetaMinusDecay* aBetaMinusChannel =
                  new G4BetaMinusDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                       daughterFloatLevel, betaType);
//              aBetaMinusChannel->DumpNuclearInfo();
//                aBetaMinusChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(aBetaMinusChannel);
                modeSumBR[1] += b;
              }
              break;

              case BetaPlus:
              {
                G4BetaPlusDecay* aBetaPlusChannel =
                  new G4BetaPlusDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                      daughterFloatLevel, betaType);
//              aBetaPlusChannel->DumpNuclearInfo();
//                aBetaPlusChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(aBetaPlusChannel);
                modeSumBR[2] += b;
              }
              break;

              case KshellEC:  // K-shell electron capture
              {
                G4ECDecay* aKECChannel =
                  new G4ECDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                daughterFloatLevel, KshellEC);
//              aKECChannel->DumpNuclearInfo();
//                aKECChannel->SetHLThreshold(halflifethreshold);
                aKECChannel->SetARM(applyARM);
                theDecayTable->Insert(aKECChannel);
                modeSumBR[3] += b;
              }
              break;

              case LshellEC:  // L-shell electron capture
              {
                G4ECDecay* aLECChannel =
                  new G4ECDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                daughterFloatLevel, LshellEC);
//              aLECChannel->DumpNuclearInfo();
//                aLECChannel->SetHLThreshold(halflifethreshold);
                aLECChannel->SetARM(applyARM);
                theDecayTable->Insert(aLECChannel);
                modeSumBR[4] += b;
              }
              break;

              case MshellEC:  // M-shell electron capture
              {
                G4ECDecay* aMECChannel =
                  new G4ECDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                daughterFloatLevel, MshellEC);
//              aMECChannel->DumpNuclearInfo();
//                aMECChannel->SetHLThreshold(halflifethreshold);
                aMECChannel->SetARM(applyARM);
                theDecayTable->Insert(aMECChannel);
                modeSumBR[5] += b;
              }
              break;

              case NshellEC:  // N-shell electron capture
              {
                G4ECDecay* aNECChannel =
                  new G4ECDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                daughterFloatLevel, NshellEC);
//              aNECChannel->DumpNuclearInfo();
//                aNECChannel->SetHLThreshold(halflifethreshold);
                aNECChannel->SetARM(applyARM);
                theDecayTable->Insert(aNECChannel);
                modeSumBR[6] += b;
              }
              break;

              case Alpha:
              {
                G4AlphaDecay* anAlphaChannel =
                  new G4AlphaDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                   daughterFloatLevel);
//              anAlphaChannel->DumpNuclearInfo();
//                anAlphaChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(anAlphaChannel);
                modeSumBR[7] += b;
              }
              break;

	      case Proton:
              {
                G4ProtonDecay* aProtonChannel =
                  new G4ProtonDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                    daughterFloatLevel);
//              aProtonChannel->DumpNuclearInfo();
//                aProtonChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(aProtonChannel);
                modeSumBR[8] += b;
              }
              break;

              case Neutron:
              {
                G4NeutronDecay* aNeutronChannel =
                  new G4NeutronDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                     daughterFloatLevel);
//              aNeutronChannel->DumpNuclearInfo();
//                aNeutronChannel->SetHLThreshold(halflifethreshold);
                theDecayTable->Insert(aNeutronChannel);
                modeSumBR[9] += b;
              }
              break;

              case BDProton:
                  // Not yet implemented
                  // G4cout << " beta-delayed proton decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case BDNeutron:
                  // Not yet implemented
                  // G4cout << " beta-delayed neutron decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case Beta2Minus:
                  // Not yet implemented
                  // G4cout << " Double beta- decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case Beta2Plus:
                  // Not yet implemented
                  // G4cout << " Double beta+ decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case Proton2:
                  // Not yet implemented
                  // G4cout << " Double proton decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case Neutron2:
                  // Not yet implemented
                  // G4cout << " Double beta- decay, a = " << a << ", b = " << b << ", c = " << c << G4endl;
                  break;
              case SpFission:
              {
                MySFDecay* aSpontFissChannel =
//                  new MySFDecay(&theParentNucleus, decayModeTotal, 0.0, 0.0);
                  new MySFDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                daughterFloatLevel);
                theDecayTable->Insert(aSpontFissChannel);
                modeSumBR[10] += b;
              }
              break;
              case Triton:
                {
                    G4TritonDecay* aTritonChannel =
                    new G4TritonDecay(&theParentNucleus, b, c*MeV, a*MeV,
                                     daughterFloatLevel);
                    //              anAlphaChannel->DumpNuclearInfo();
                    //                anAlphaChannel->SetHLThreshold(halflifethreshold);
                    theDecayTable->Insert(aTritonChannel);
                    modeSumBR[11] += b;
                }
              break;

              case RDM_ERROR:

              default:
                G4Exception("G4RadioactiveDecay::LoadDecayTable()", "HAD_RDM_000",
                            FatalException, "Selected decay mode does not exist");
            }  // switch
          }  // line < 72
        }  // if char == P
      }  // if char != #
    }  // While


    // Go through the decay table and make sure that the branching ratios are
    // correctly normalised.

    G4VDecayChannel* theChannel = 0;
    G4NuclearDecay* theNuclearDecayChannel = 0;
    G4String mode = "";

    G4double theBR = 0.0;
    for (G4int i = 0; i < theDecayTable->entries(); i++) {
      theChannel = theDecayTable->GetDecayChannel(i);
      theNuclearDecayChannel = static_cast<G4NuclearDecay*>(theChannel);
      theDecayMode = theNuclearDecayChannel->GetDecayMode();

      if (theDecayMode != IT) {
	theBR = theChannel->GetBR();
	theChannel->SetBR(theBR*modeTotalBR[theDecayMode]/modeSumBR[theDecayMode]);
      }
    }
  }  // decay file exists

  DecaySchemeFile.close();

  if (!found && levelEnergy > 0) {
    // Case where IT cascade for excited isotopes has no entries in RDM database
    // Decay mode is isomeric transition.
    G4ITDecay* anITChannel = new G4ITDecay(&theParentNucleus, 1.0, 0.0, 0.0,
                                           photonEvaporation);
//    anITChannel->SetHLThreshold(halflifethreshold);
    anITChannel->SetARM(applyARM);
    theDecayTable->Insert(anITChannel);
  }

  if (theDecayTable && GetVerboseLevel() > 1) {
    if(theDecayTable->entries()) theDecayTable->DumpInfo();
  }

#ifdef G4MULTITHREADED
  //(*master_dkmap)[key] = theDecayTable;                  // store in master library
  lk.unlock(); //explicit unlock
  //G4cout << "LoadDecayTable unlocking its mutex" << G4endl;
#endif
  return theDecayTable;
}

void MyRadioactiveDecayBase::LoadAllDecayTables()
{
  std::string path = dirPath;
  for (const auto & entry : std::filesystem::directory_iterator(path)) {
      //G4cout << entry.path().string() << G4endl;
      //G4cout << entry.path().extension().string() << G4endl;

      std::string stem = entry.path().stem().string();
      std::string extension = entry.path().extension().string();

      stem.erase(0,1); //remove the z character
      extension.erase(0,2); //remove the a character

      G4int Z = atoi(stem.data());
      G4int A = atoi(extension.data());

      //G4cout << "Z = " << Z << " : A = " << A << G4endl;

      if(Z && A) {
        G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,0.);
        GetDecayTable(ion);
      }
  }
        

}

void
MyRadioactiveDecayBase::AddUserDecayDataFile(G4int Z, G4int A, G4String filename)
{
  if (Z < 1 || A < 2) G4cout << "Z and A not valid!" << G4endl;

  std::ifstream DecaySchemeFile(filename);
  if (DecaySchemeFile) {
    G4int ID_ion = A*1000 + Z;
    theUserRadioactiveDataFiles[ID_ion] = filename;
  } else {
    G4cout << "The file " << filename << " does not exist!" << G4endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DecayIt                                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

G4VParticleChange*
MyRadioactiveDecayBase::DecayIt(const G4Track& theTrack, const G4Step&)
{
  //G4cout << "MyRadioactiveDecayBase::DecayIt weight : " << theTrack.GetWeight() << G4endl;
  //G4cout << " : " << theTrack.GetParticleDefinition()->GetParticleName() << G4endl;
  // Initialize G4ParticleChange object, get particle details and decay table
  fParticleChangeForRadDecay.Initialize(theTrack);
  fParticleChangeForRadDecay.ProposeWeight(theTrack.GetWeight());
  const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
  const G4ParticleDefinition* theParticleDef = theParticle->GetDefinition();

  // First check whether RDM applies to the current logical volume
  if (!isAllVolumesMode) {
    if (!std::binary_search(ValidVolumes.begin(), ValidVolumes.end(),
                     theTrack.GetVolume()->GetLogicalVolume()->GetName())) {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) {
        G4cout <<"G4RadioactiveDecay::DecayIt : "
               << theTrack.GetVolume()->GetLogicalVolume()->GetName()
               << " is not selected for the RDM"<< G4endl;
        G4cout << " There are " << ValidVolumes.size() << " volumes" << G4endl;
        G4cout << " The Valid volumes are " << G4endl;
        for (size_t i = 0; i< ValidVolumes.size(); i++)
                                  G4cout << ValidVolumes[i] << G4endl;
      }
#endif
      fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

      // Kill the parent particle.
      fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
      fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChangeForRadDecay;
    }
  }

  // Now check if particle is valid for RDM
  if (!(IsApplicable(*theParticleDef) ) ) { 
    // Particle is not an ion or is outside the nucleuslimits for decay

    if (GetVerboseLevel() > 0) {
      G4cout << "G4RadioactiveDecay::DecayIt : "
             << theParticleDef->GetParticleName() 
             << " is not an ion or is outside (Z,A) limits set for the decay. " 
             << " Set particle change accordingly. "
             << G4endl;
    }
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

    // Kill the parent particle
    fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForRadDecay;
  }

  G4DecayTable* theDecayTable = GetDecayTable(theParticleDef);

  if (theDecayTable == 0 || theDecayTable->entries() == 0) {
    // No data in the decay table.  Set particle change parameters
    // to indicate this.
    if (GetVerboseLevel() > 0) {
      G4cout << "G4RadioactiveDecay::DecayIt : "
             << "decay table not defined for "
             << theParticleDef->GetParticleName() 
             << ". Set particle change accordingly. "
             << G4endl;
    }
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);

    // Kill the parent particle.
    fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForRadDecay;

  } else { 
    // Data found.  Try to decay nucleus

    // Decay without variance reduction 
    DecayAnalog(theTrack);
    return &fParticleChangeForRadDecay ;
  }
} 


void MyRadioactiveDecayBase::DecayAnalog(const G4Track& theTrack)
{
  const G4DynamicParticle* theParticle = theTrack.GetDynamicParticle();
  const G4ParticleDefinition* theParticleDef = theParticle->GetDefinition();
  G4DecayProducts* products = 0;
  G4double weight = theTrack.GetWeight();
  //G4cout << "MyRadioactiveDecayBase getweight = " << weight << G4endl;
  if(!applyBRbias) {
    products = DoDecay(*theParticleDef);
  } else {
    products = DoDecayBRbias(*theParticleDef, weight);
  }

  // Check if the product is the same as input and kill the track if
  // necessary to prevent infinite loop (11/05/10, F.Lei)
  if (products->entries() == 1) {
    fParticleChangeForRadDecay.SetNumberOfSecondaries(0);
    fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill);
    fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(0.0);
    ClearNumberOfInteractionLengthLeft();
    return;
  }

  G4double energyDeposit = 0.0;
  G4double finalGlobalTime = theTrack.GetGlobalTime();
  G4double finalLocalTime = theTrack.GetLocalTime();

  // Get parent particle information and boost the decay products to the
  // laboratory frame

  // ParentEnergy used for the boost should be the total energy of the nucleus
  // of the parent ion without the energy of the shell electrons
  // (correction for bug 1359 by L. Desorgher)
  G4double ParentEnergy = theParticle->GetKineticEnergy()
                        + theParticle->GetParticleDefinition()->GetPDGMass();
  G4ThreeVector ParentDirection(theParticle->GetMomentumDirection());

  if (theTrack.GetTrackStatus() == fStopButAlive) {
    // this condition seems to be always True, further investigation is needed (L.Desorgher)

    // The particle is decayed at rest
    // Since the time is for the particle at rest, need to add additional time
    // lapsed between particle coming to rest and the actual decay.  This time
    // is sampled with the mean-life of the particle.  Need to protect the case 
    // PDGTime < 0.  (F.Lei 11/05/10)
    G4double temptime = -std::log(G4UniformRand() ) *
                        theParticleDef->GetPDGLifeTime();
    if (temptime < 0.) temptime = 0.;
    finalGlobalTime += temptime;
    finalLocalTime += temptime;
    energyDeposit += theParticle->GetKineticEnergy();
  }
  products->Boost(ParentEnergy, ParentDirection);

  // Add products in theParticleChangeForRadDecay.
  G4int numberOfSecondaries = products->entries();
  fParticleChangeForRadDecay.SetNumberOfSecondaries(numberOfSecondaries);

  if (GetVerboseLevel() > 1) {
    G4cout << "G4RadioactiveDecay::DecayAnalog: Decay vertex :";
    G4cout << " Time: " << finalGlobalTime/ns << "[ns]";
    G4cout << " X:" << (theTrack.GetPosition()).x() /cm << "[cm]";
    G4cout << " Y:" << (theTrack.GetPosition()).y() /cm << "[cm]";
    G4cout << " Z:" << (theTrack.GetPosition()).z() /cm << "[cm]";
    G4cout << G4endl;
    G4cout << "G4Decay::DecayIt : decay products in Lab. Frame" << G4endl;
    products->DumpInfo();
    products->IsChecked();
  }

  for (G4int index = 0; index < numberOfSecondaries; index++) {
    G4Track* secondary = new G4Track(products->PopProducts(), finalGlobalTime,
                                     theTrack.GetPosition() );
    secondary->SetCreatorModelIndex(theRadDecayMode);
    secondary->SetWeight(secondary->GetWeight()*weight);
    //Change for atomics relaxation
    if (theRadDecayMode == IT  && index>0){
    if (index == numberOfSecondaries-1) secondary->SetCreatorModelIndex(IT);
         else secondary->SetCreatorModelIndex(30);
    }
    else if (theRadDecayMode >= KshellEC && theRadDecayMode <= NshellEC
                     			                             && index <numberOfSecondaries-1){
        secondary->SetCreatorModelIndex(30);
    }
    secondary->SetGoodForTrackingFlag();
    secondary->SetTouchableHandle(theTrack.GetTouchableHandle());
    fParticleChangeForRadDecay.AddSecondary(secondary);
  }

  delete products;

  // Kill the parent particle
  fParticleChangeForRadDecay.ProposeTrackStatus(fStopAndKill) ;
  fParticleChangeForRadDecay.ProposeLocalEnergyDeposit(energyDeposit);
  fParticleChangeForRadDecay.ProposeLocalTime(finalLocalTime);

  // Reset NumberOfInteractionLengthLeft.
  ClearNumberOfInteractionLengthLeft();
}


G4DecayProducts*
MyRadioactiveDecayBase::DoDecay(const G4ParticleDefinition& theParticleDef)
{
  G4DecayProducts* products = 0;
  G4DecayTable* theDecayTable = GetDecayTable(&theParticleDef);
  // Choose a decay channel.
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) G4cout << "Select a channel..." << G4endl;
#endif

  // G4DecayTable::SelectADecayChannel checks to see if sum of daughter masses
  // exceeds parent mass. Pass it the parent mass + maximum Q value to account
  // for difference in mass defect.
  G4double parentPlusQ = theParticleDef.GetPDGMass() + 30.*MeV;
  G4VDecayChannel* theDecayChannel = theDecayTable->SelectADecayChannel(parentPlusQ);
  theRadDecayMode = (static_cast<G4NuclearDecay*>(theDecayChannel))->GetDecayMode();

  if (theDecayChannel == 0) {
    // Decay channel not found.
    G4ExceptionDescription ed;
    ed << " Cannot determine decay channel for " << theParticleDef.GetParticleName() << G4endl;
    G4Exception("G4RadioactiveDecay::DoDecay", "HAD_RDM_013",
                FatalException, ed);
  } else {
    // A decay channel has been identified, so execute the DecayIt.
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cout << "G4RadioactiveDecay::DoIt : selected decay channel addr: "
             << theDecayChannel << G4endl;
    }
#endif
    products = theDecayChannel->DecayIt(theParticleDef.GetPDGMass() );

    // Apply directional bias if requested by user
    CollimateDecay(products);
  }

  return products;
}

G4DecayProducts*
MyRadioactiveDecayBase::DoDecayBRbias(const G4ParticleDefinition& theParticleDef, G4double &weigth)
{
  G4DecayProducts* products = 0;
  G4DecayTable* theDecayTable = GetDecayTable(&theParticleDef);
  // Choose a decay channel.
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) G4cout << "Select a channel..." << G4endl;
#endif

  if (GetVerboseLevel() > 1) theDecayTable ->DumpInfo();
  //Select a decay channel
  //G4DecayTable* decayTable = GetDecayTable1(parentNucleus);
  G4int ndecaych = G4int(theDecayTable->entries()*G4UniformRand());
  G4VDecayChannel* theDecayChannel = theDecayTable->GetDecayChannel(ndecaych);

  if (theDecayChannel == 0) {
    // Decay channel not found.

    if (GetVerboseLevel() > 0) {
      G4cout << " G4RadioactiveDecay::DoIt : cannot determine decay channel ";
      G4cout << " for this nucleus; decay as if no biasing active. ";
      G4cout << G4endl;
      theDecayTable ->DumpInfo();
    }

    products = DoDecay(theParticleDef);  // DHW 6 Dec 2010 - do decay as if no biasing
                                            //           to avoid deref of temppprods = 0
  } else {
    // A decay channel has been identified, so execute the DecayIt.
    G4double tempmass = theParticleDef.GetPDGMass();
    products = theDecayChannel->DecayIt(tempmass);
    weigth *= (theDecayChannel->GetBR())*(theDecayTable->entries());
    
    CollimateDecay(products);

  }

  return products;
}


// Apply directional bias for "visible" daughters (e+-, gamma, n, p, alpha)

void MyRadioactiveDecayBase::CollimateDecay(G4DecayProducts* products) {
  if (origin == forceDecayDirection) return;	// No collimation requested
  if (180.*deg == forceDecayHalfAngle) return;
  if (0 == products || 0 == products->entries()) return;

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 0) G4cout << "Begin of CollimateDecay..." << G4endl;
#endif

  // Particles suitable for directional biasing (for if-blocks below)
  static const G4ParticleDefinition* electron = G4Electron::Definition();
  static const G4ParticleDefinition* positron = G4Positron::Definition();
  static const G4ParticleDefinition* neutron  = G4Neutron::Definition();
  static const G4ParticleDefinition* gamma    = G4Gamma::Definition();
  static const G4ParticleDefinition* alpha    = G4Alpha::Definition();
  static const G4ParticleDefinition* triton  = G4Triton::Definition();
  static const G4ParticleDefinition* proton   = G4Proton::Definition();

  G4ThreeVector newDirection;		// Re-use to avoid memory churn
  for (G4int i=0; i<products->entries(); i++) {
    G4DynamicParticle* daughter = (*products)[i];
    const G4ParticleDefinition* daughterType =
                                  daughter->GetParticleDefinition();
    if (daughterType == electron || daughterType == positron ||
	daughterType == neutron || daughterType == gamma ||
	daughterType == alpha || daughterType == triton || daughterType == proton) CollimateDecayProduct(daughter);
  }
}

void MyRadioactiveDecayBase::CollimateDecayProduct(G4DynamicParticle* daughter) {
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "CollimateDecayProduct for daughter "
	   << daughter->GetParticleDefinition()->GetParticleName() << G4endl;
  }
#endif

  G4ThreeVector collimate = ChooseCollimationDirection();
  if (origin != collimate) daughter->SetMomentumDirection(collimate);
}


// Choose random direction within collimation cone

G4ThreeVector MyRadioactiveDecayBase::ChooseCollimationDirection() const {
  if (origin == forceDecayDirection) return origin;	// Don't do collimation
  if (forceDecayHalfAngle == 180.*deg) return origin;

  G4ThreeVector dir = forceDecayDirection;

  // Return direction offset by random throw
  if (forceDecayHalfAngle > 0.) {
    // Generate uniform direction around central axis
    G4double phi = 2.*pi*G4UniformRand();
    G4double cosMin = std::cos(forceDecayHalfAngle);
    G4double cosTheta = (1.-cosMin)*G4UniformRand() + cosMin;	// [cosMin,1.)
    
    dir.setPhi(dir.phi()+phi);
    dir.setTheta(dir.theta()+std::acos(cosTheta));
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    G4cout << " ChooseCollimationDirection returns " << dir << G4endl;
#endif

  return dir;
}

