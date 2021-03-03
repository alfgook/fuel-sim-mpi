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
// $Id: DetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4AssemblyVolume;
class G4NistManager;
class GB01BOptrMultiParticleChangeCrossSection;
class G4GenericMessenger;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(G4String aMaterialFile);
    virtual ~DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();

    void ChangeBiasForParticle( G4String particleName, G4double biasFactor);

    // get methods
    //
	G4Material* GetDetectorMaterial() const;
    virtual void ConstructSDandField();
     
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    G4VPhysicalVolume* DefineVolumesBWR();
    void BuildDetectors(G4LogicalVolume *worldLV);
    G4AssemblyVolume* EJ309_5x5inch(G4int copyNbr, const char* name);
    G4AssemblyVolume* EJ309_3x3inch(G4int copyNbr, const char* name);
    G4AssemblyVolume* EJ309_1x2inch(G4int copyNbr, const char* name);
    G4AssemblyVolume* EJ276_detector(G4int copyNbr, const char* name, G4double size);
	
    int ReadMCNPmatCard(const char *FileName);
  
    // data members

    G4String fMaterialFile;

	G4NistManager* nistManager;
        
    G4Material*        fWorldMaterial;
    G4Material*        fVacuum;
    G4Material*        fScintilatorMat; //pointer to the Scintilator Material
    G4Material*        fAlu; // pointer to the Aluminium material
    G4Material*        fEJ309; // pointer to the LS301 material
    G4Material*        fEJ276; // pointer to the EJ276 plastic scintillator material
    G4Material*        fBoroSilicate; //pointer to the Scintilators Light Guide Material
    //G4Material*        fPb; 
    G4Material*        fCladingMat;
    G4Material*        fFuelMat;
    G4Material*        fNeutronShield;
    G4Material*        fCastIron;
    
    G4UserLimits*      fStepLimit;       // pointer to user step limits
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

    //the vector of GB01BOptrMultiParticleChangeCrossSection* is nescessary to make it possible to change the bias from UI command
    //otherwise the bias is only changed for the master
    std::vector<GB01BOptrMultiParticleChangeCrossSection*> fXSbiasVector;
    G4GenericMessenger* fMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

