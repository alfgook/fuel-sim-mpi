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
//
#include "ParallelWorldForSplitting.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

#include "G4SDManager.hh"

#include "G4PVDivision.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ParallelWorldForSplitting::ParallelWorldForSplitting(G4String worldName)
  : G4VUserParallelWorld(worldName)
{
  fMessenger = new G4GenericMessenger(this,"/ParallelWorldForSplitting/", "...doc...");

  fMessenger->DeclareProperty("InnnerRad",
                            fRad,
                            "set R; the inner radius of the parallel world cylinder");
  fMessenger->DeclareProperty("HalfHeight",
                            fHalfHeigt,
                            "set H; the half height of the parallel world cylinder");
  fMessenger->DeclareMethod("ApplyChanges",
                            &ParallelWorldForSplitting::ApplyChanges,
                            "Apply the new values (R,H) to the geometry");
  /*fMessenger->DeclareMethod("Inactive",
                            &SteppingAction::InactiveSplitting,
                            "Inactivate the event splitting");*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ParallelWorldForSplitting::~ParallelWorldForSplitting()
{;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ParallelWorldForSplitting::ApplyChanges()
{
  if(fRad>=fRadMax) {
    fRad = fRadMax - 1.;
    G4cout << "fRad >= shield outer radius: it will be set to " << fRad << G4endl;
  }
  if(fRad==fRadMin) {
    fRad = fRadMin-1;
    G4cout << "fRad == shield inner radius: it will be set to " << fRad << G4endl;
  }
  if(fHalfHeigt>=fHalfHeigtMax) {
    fHalfHeigt = fHalfHeigtMax;
    G4cout << "fHalfHeigt >= shield outer: it will be set to " << fHalfHeigt << G4endl;
  }
  if(fHalfHeigt==fHalfHeigtMin) {
    fHalfHeigt = fHalfHeigtMin-1;
    G4cout << "fHalfHeigt == shield inner: it will be set to " << fHalfHeigt << G4endl;
  }

  G4cout << "fRad = " << fRad << G4endl;
  G4cout << "fHalfHeigt = " << fHalfHeigt << G4endl;
  fInnerSolid->SetOuterRadius(fRad);
  //fInnerSolid->SetInnerRadius(fRad);
  fInnerSolid->SetZHalfLength(fHalfHeigt);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParallelWorldForSplitting::Construct()
{
  // -- Inform about construction:
  // -- (fWorldName is a protected data member of the base parallel world class)
  G4cout << "Parallel World `" << fWorldName << "' constructed." << G4endl;

  // -------------------------
  //  Build parallel geometry:
  // -------------------------
  
  // -- Obtain clone of mass geometry world from GetWorld() base class utility:
  G4VPhysicalVolume* physicalParallelWorld = GetWorld();
  G4LogicalVolume*    logicalParallelWorld = physicalParallelWorld->GetLogicalVolume();
  
  // -- 1) get back the solid used to create the concrete shield:
  //       ------------------------------------------------------
  
  // -- get back the logical volume of the shield, using its name:
  G4LogicalVolume* shieldLogical =
    G4LogicalVolumeStore::GetInstance()->GetVolume("CaskCylinderLV");
  G4Tubs* shieldSolid = (G4Tubs*) shieldLogical->GetSolid();

  G4double rOut1 = shieldSolid->GetOuterRadius();
  G4double halfZ1 = shieldSolid->GetZHalfLength();
  G4double phiStart1 = shieldSolid->GetStartPhiAngle();
  G4double phiDelta1 = shieldSolid->GetDeltaPhiAngle();

  G4LogicalVolume* cavityLogical =
    G4LogicalVolumeStore::GetInstance()->GetVolume("CaskCavityLV");
  G4Tubs* cavitySolid = (G4Tubs*) cavityLogical->GetSolid();
  G4double rIn1 = cavitySolid->GetOuterRadius();

  G4cout << "    rIn1 = "<< rIn1 << G4endl;
  G4cout << "    rOut1 = "<< rOut1 << G4endl;
  G4cout << "    halfZ1 = "<< halfZ1 << G4endl;

  fRadMin = rIn1;
  fRadMax = rOut1;
  fRad = rIn1 + 0.5*(rOut1-rIn1);

  fHalfHeigtMin = halfZ1 -1;
  fHalfHeigtMax = halfZ1+310.;
  fHalfHeigt = halfZ1 + 310./2.;

  G4cout << "    rad = "<< fRad << G4endl;
  G4cout << "    HalfHeigt = "<< fHalfHeigt << G4endl;
  G4cout << "    phiStart = "<< phiStart1 << G4endl;
  G4cout << "    phiDelta = "<< phiDelta1 << G4endl;

  fInnerSolid = new G4Tubs("Inner.solid",0.,fRad,fHalfHeigt,phiStart1,phiDelta1);
  
  G4LogicalVolume *InnerLogical =
    new G4LogicalVolume(fInnerSolid,                 // its solid
                        nullptr,                     // no material
                        "Inner.logical");  // its name
 
  
  // -- get back the translation
  G4VPhysicalVolume*
  shieldPhysical = G4PhysicalVolumeStore::GetInstance()->GetVolume("CaskCylinderPV");
  G4ThreeVector translation = shieldPhysical->GetObjectTranslation();
  
  new G4PVPlacement( nullptr,                           // no rotation
                     G4ThreeVector(0,0,0),                       // translate as for the shield
                     InnerLogical,            // its logical volume
                     "Inner.physical",        // its name
                     logicalParallelWorld,              // its mother  volume
                     false,                             // no boolean operation
                     0,                                // copy number
                     true);                             //check overlap

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParallelWorldForSplitting::ConstructSD()
{
 
}
