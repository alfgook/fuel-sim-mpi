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
/// \file ParticleSourceFromFile.cc
/// \brief Implementation of the PrimaryGenerator1 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ParticleSourceFromFile.hh"

#include "G4Event.hh"
#include "Analysis.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "AnalysisMPI.hh"

#ifndef NOT_USING_MPI
#include "G4MPImanager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ParticleSourceFromFile::ParticleSourceFromFile(G4String FileName)
: G4VPrimaryGenerator()
{

  fNbrGenerated = 0;
  fSplitFactor = 2;
  fWeight = 1./G4double(fSplitFactor);
  
  eofReached = true;
  OpenFile(FileName);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ParticleSourceFromFile::~ParticleSourceFromFile()
{

   G4cout << "~ParticleSourceFromFile1" << G4endl;
  // Delete analysis reader
  delete G4AnalysisReader::Instance(); 
   G4cout << "~ParticleSourceFromFile2" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void  ParticleSourceFromFile::OpenFile(G4String FileName)
{
  eofReached = false;
  #ifndef NOT_USING_MPI
  //G4cout << "THIS IS AN MPI SESSION" << G4endl;
  G4int rank = G4MPImanager::GetManager()->GetRank();
  FileName += "_r" + std::to_string(rank);
  #endif
  // set up the ntuple reader
  G4AnalysisReader* analysisReader = G4AnalysisReader::Instance();
  analysisReader->SetFileName(FileName);
  G4cout << "input filename = " << FileName << G4endl;
  G4int ntupleId = analysisReader->GetNtuple("SplitPoints");
  if(ntupleId==G4Analysis::kInvalidId) {
    G4cout << "*************************************" << G4endl;
    G4cout << "* Invalid ntuple                    *" << G4endl;
    G4cout << "* prob. the file:                   *" << G4endl;
    //G4cout << "* " << FileName << " does not exist *" << G4endl;
    G4cout << "* " << analysisReader->GetFileName() << " does not exist *" << G4endl;
    G4cout << "*************************************" << G4endl;
    eofReached = true;
    return;
  }
  if ( ntupleId >= 0 ) {
    analysisReader->SetNtupleIColumn("eventID", currentEntry.eventID);
    analysisReader->SetNtupleIColumn("PDGcode", currentEntry.PDGcode);
    analysisReader->SetNtupleFColumn("Time", currentEntry.Time);
    analysisReader->SetNtupleFColumn("Energy", currentEntry.Energy);
    analysisReader->SetNtupleFColumn("posX", currentEntry.posX);
    analysisReader->SetNtupleFColumn("posY", currentEntry.posY);
    analysisReader->SetNtupleFColumn("posZ", currentEntry.posZ);
    analysisReader->SetNtupleFColumn("dirX", currentEntry.dirX);
    analysisReader->SetNtupleFColumn("dirY", currentEntry.dirY);
    analysisReader->SetNtupleFColumn("dirZ", currentEntry.dirZ);
    analysisReader->SetNtupleFColumn("InitX", currentEntry.initX);
    analysisReader->SetNtupleFColumn("InitY", currentEntry.initY);
    analysisReader->SetNtupleFColumn("InitZ", currentEntry.initZ);
    analysisReader->SetNtupleFColumn("weight", currentEntry.weight);
  }

  // read the first event
  if(analysisReader->GetNtupleRow()) {
    TreeEntries.emplace_back(currentEntry);
    /*fPDGcodes.emplace_back(PDGcode);
    fTimes.emplace_back(Time);
    fEnergies.emplace_back(Energy);
    fPositions.emplace_back(G4ThreeVector(posX,posY,posZ));
    fDirections.emplace_back(G4ThreeVector(dirX,dirY,dirZ));
    fWeightIn.emplace_back(weightIn);*/
  } else {
    eofReached = true;
  }
  // last entry in the vectors will be from the next event
  G4int currentEventID = currentEntry.eventID;
  do {
    if(analysisReader->GetNtupleRow()) {
      TreeEntries.emplace_back(currentEntry);
      /*fPDGcodes.emplace_back(PDGcode);
      fTimes.emplace_back(Time);
      fEnergies.emplace_back(Energy);
      fPositions.emplace_back(G4ThreeVector(posX,posY,posZ));
      fDirections.emplace_back(G4ThreeVector(dirX,dirY,dirZ));
      fWeightIn.emplace_back(weightIn);*/
    } else {
      eofReached = true;
      break; //end of file reached?
    }
  } while (currentEntry.eventID==currentEventID);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String ParticleSourceFromFile::GetFileName() const
{
  return G4AnalysisReader::Instance()->GetFileName();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool ParticleSourceFromFile::GetNextEvent()
{
  if(eofReached) {
    G4cout << "End of Input file reached: stopping the run" << G4endl;
    G4RunManager::GetRunManager()->AbortRun();
    return false;
    //should work even in multithreaded: G4RunManager::GetRunManager() returns the thread local G4WorkerRunManager
  }

  if(fNbrGenerated>=fSplitFactor) {
    fNbrGenerated = 0;

    TreeEntries[0] = TreeEntries.back();
    TreeEntries.resize(1);

    G4AnalysisReader* analysisReader = G4AnalysisReader::Instance();

    G4int currentEventID = currentEntry.eventID; //this is the last one read from the former call it corresponds to the next event
    do {
      if(analysisReader->GetNtupleRow()) {
        TreeEntries.emplace_back(currentEntry);
      } else {
        eofReached = true;
        break; //end of file reached?
      }
    } while (currentEntry.eventID==currentEventID);

    G4int nonScatteredGamma = -1;
    G4int nonScatteredNeutron = -1;
    G4int ScatteredGamma = -1;
    G4int ScatteredNeutron = -1;
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    for(size_t i=0;i<TreeEntries.size()-1;i++) { // the last elements are from the next event
      if(TreeEntries.at(i).PDGcode==22) {
        G4ThreeVector v1( TreeEntries.at(i).posX-TreeEntries.at(i).initX,
                          TreeEntries.at(i).posY-TreeEntries.at(i).initY,
                          TreeEntries.at(i).posZ-TreeEntries.at(i).initZ);

        G4ThreeVector v2( TreeEntries.at(i).dirX,
                          TreeEntries.at(i).dirY,
                          TreeEntries.at(i).dirZ);
        //if(v1.angle(v2)<0.0027777778*CLHEP::twopi) { // smaller than one degree
        if(v1.angle(v2)<1E-6) {
          analysisManager->FillH1(8,TreeEntries.at(i).Energy); // non-scaterred gammas
          nonScatteredGamma = i;
        } else {
          analysisManager->FillH1(7,TreeEntries.at(i).Energy); // scaterred gammas
          ScatteredGamma = i;
        }
        
      } else if(TreeEntries.at(i).PDGcode==2112) {
        G4ThreeVector v1( TreeEntries.at(i).posX-TreeEntries.at(i).initX,
                          TreeEntries.at(i).posY-TreeEntries.at(i).initY,
                          TreeEntries.at(i).posZ-TreeEntries.at(i).initZ);

        G4ThreeVector v2( TreeEntries.at(i).dirX,
                          TreeEntries.at(i).dirY,
                          TreeEntries.at(i).dirZ);
        //if(v1.angle(v2)<0.0027777778*CLHEP::twopi) { // smaller than one degree
        if(v1.angle(v2)<1E-6) {
          analysisManager->FillH1(6,TreeEntries.at(i).Energy); // non-scattered neutrons
          analysisManager->FillH2(0,TreeEntries.at(i).initX,TreeEntries.at(i).initY); // non-scattered neutrons
          nonScatteredNeutron = i;
        } else {
          analysisManager->FillH1(5,TreeEntries.at(i).Energy); // scattered neutrons
          ScatteredNeutron = i;
        }
      }
    }

    //Set intial position of decay
    AnalysisMPI::GetAnalysis()->SetInitial(TreeEntries.at(0).initX,TreeEntries.at(0).initY,TreeEntries.at(0).initZ);

    //Fill histograms
    const G4double EnergyThNeutron = 2.;
    const G4double EnergyThGamma = 2.;
    if(nonScatteredGamma>=0 && nonScatteredNeutron>=0) { // neither neutron nor gamma has been scattered
      G4double Eg =  TreeEntries.at(nonScatteredGamma).Energy;
      G4double En =  TreeEntries.at(nonScatteredNeutron).Energy;
      if(Eg>EnergyThGamma && En>EnergyThNeutron)
        analysisManager->FillH1(9,TreeEntries.at(nonScatteredNeutron).Time-TreeEntries.at(nonScatteredGamma).Time);
    }
    if(ScatteredGamma>=0 && ScatteredNeutron>=0) { // both neutron and gamma has been scattered
      G4double Eg =  TreeEntries.at(ScatteredGamma).Energy;
      G4double En =  TreeEntries.at(ScatteredNeutron).Energy;
      if(Eg>EnergyThGamma && En>EnergyThNeutron)
        analysisManager->FillH1(10,TreeEntries.at(ScatteredNeutron).Time-TreeEntries.at(ScatteredGamma).Time);
    }
    if(ScatteredGamma>=0 && nonScatteredNeutron>=0) { // only gamma has been scattered
      G4double Eg =  TreeEntries.at(ScatteredGamma).Energy;
      G4double En =  TreeEntries.at(nonScatteredNeutron).Energy;
      if(Eg>EnergyThGamma && En>EnergyThNeutron)
        analysisManager->FillH1(10,TreeEntries.at(nonScatteredNeutron).Time-TreeEntries.at(ScatteredGamma).Time);
    }
    if(nonScatteredGamma>=0 && ScatteredNeutron>=0) { // only neutron has been scattered
      G4double Eg =  TreeEntries.at(nonScatteredGamma).Energy;
      G4double En =  TreeEntries.at(ScatteredNeutron).Energy;
      if(Eg>EnergyThGamma && En>EnergyThNeutron)
        analysisManager->FillH1(10,TreeEntries.at(ScatteredNeutron).Time-TreeEntries.at(nonScatteredGamma).Time);
    }
  }
  ++fNbrGenerated;

  /*if(!EventHasGamma()) {
    //fWeight = 1;
    fNbrGenerated = fSplitFactor; //this line ensures that only events that contain at least one gamma-ray will be split

    // russian roulette on events without gammas (killing off 1/2 the events where no gamma is present)
    G4double xx = G4UniformRand();
    if(xx > 1./((G4double) fSplitFactor)) { //event killed
      GetNextEvent();
    } else { // event survives
      fWeight = fSplitFactor;
    }

  } else {
    fWeight = 1./((G4double) fSplitFactor);
  }*/

  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool ParticleSourceFromFile::EventHasGamma()
{
  for(size_t i=0;i<TreeEntries.size()-1;i++) { // the last elements are from the next event
    if(TreeEntries.at(i).PDGcode==22) {
      return true;
    }
  }
  return false;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ParticleSourceFromFile::Print(G4int ii) const
{
    G4ThreeVector position = G4ThreeVector(TreeEntries.at(ii).posX,TreeEntries.at(ii).posY,TreeEntries.at(ii).posZ);
    G4ThreeVector direction = G4ThreeVector(TreeEntries.at(ii).dirX,TreeEntries.at(ii).dirY,TreeEntries.at(ii).dirZ);

    G4cout << "-------" << G4endl;
    G4cout << "position " << position << G4endl;
    G4cout << "time " << TreeEntries.at(ii).Time << G4endl;
    G4cout << "PDG code " << TreeEntries.at(ii).Time << G4endl;
    G4cout << "direction " << direction << G4endl;
    G4cout << "weight " << TreeEntries.at(ii).weight << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ParticleSourceFromFile::SetSplitFactor(G4int aSplitFactor)
{
  fSplitFactor = aSplitFactor;
  fWeight = 1./((G4double) fSplitFactor);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParticleSourceFromFile::GeneratePrimaryVertex(G4Event* event)
{
  if(GetNextEvent()) {

    /*if(EventHasGamma()) {
      G4cout << "this event has a gamma-ray" << G4endl;
    } else {
      G4cout << "this event doesn't have a gamma-ray" << G4endl;
    }*/

    for(size_t i=0;i<TreeEntries.size()-1;i++) { // the last elements are from the next event

      //Print(i);

      G4ThreeVector position =
        G4ThreeVector(TreeEntries.at(i).posX,TreeEntries.at(i).posY,TreeEntries.at(i).posZ);

      G4PrimaryVertex *theVertex = new G4PrimaryVertex(position,TreeEntries.at(i).Time);

      G4PrimaryParticle* particle =
        new G4PrimaryParticle(TreeEntries.at(i).PDGcode,TreeEntries.at(i).dirX,TreeEntries.at(i).dirY,TreeEntries.at(i).dirZ);

      particle->SetKineticEnergy(TreeEntries.at(i).Energy);
      particle->SetWeight(fWeight*TreeEntries.at(i).weight);

      theVertex->SetPrimary(particle);
      event->AddPrimaryVertex(theVertex); 
    }
  }
  
}