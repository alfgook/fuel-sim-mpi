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
// $Id: DetectorConstruction.cc 77601 2013-11-26 17:08:44Z gcosmo $
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "Config.h"

#include "DetectorConstruction.hh"

#include "ScintilatorSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Trd.hh"
#include "G4Torus.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AutoDelete.hh"

#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4SDManager.hh"

#include "GB01BOptrMultiParticleChangeCrossSection.hh"
#include "G4GenericMessenger.hh"

#include "MaterialFile.hh"
#ifndef NOT_USING_MPI
	#include "G4MPImanager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4String aMaterialFile)
 : G4VUserDetectorConstruction(),
    fStepLimit(NULL), 
    fCheckOverlaps(false)
{
	//G4cout << "DetectorConstruction::DetectorConstruction" << G4endl;
	fMaterialFile = aMaterialFile;
	if(!fMaterialFile.size()) {
		G4String dir(INPUT_DIR);
		fMaterialFile = dir + "/SKB-TR-10-13-PWR_MCNP_matls.inp";
	}

	fMessenger = new G4GenericMessenger(this,"/Detector/XSbias/", "control the biasing of cross section");
  	auto& cmd = fMessenger->DeclareMethod("ChangeBiasForParticle",
                            &DetectorConstruction::ChangeBiasForParticle,
                            "Change the value of the bias factor (par 2) for patricle (par 1),the operator must allready exist (declared in DetectorConstruction.");
  	cmd.SetParameterName("particle",false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
	delete fStepLimit;
	delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ChangeBiasForParticle( G4String particleName, G4double biasFactor)
{
	for(auto it=fXSbiasVector.begin();it!=fXSbiasVector.end();++it) {
		GB01BOptrMultiParticleChangeCrossSection *aXSbias = *it;
		if(!aXSbias) {
		  G4ExceptionDescription ed;
	      ed << "XS biasing not initilized" << G4endl;
	      G4Exception("DetectorConstruction::ChangeBiasForParticle(...)",
	                  "DetectorConstruction.01",
	                  JustWarning,
	                  ed);
	      return;
		}
		aXSbias->ChangeBiasForParticle(particleName,biasFactor);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  
  nistManager = G4NistManager::Instance();
  
  DefineMaterials();

//	G4cout << "======================================================" << G4endl;
//	G4cout << "========== DetectorConstruction::DefineVolumes() =========" << G4endl;
//	G4cout << "======================================================" << G4endl;
  // Define volumes
  //return DefineVolumes();
  return DefineVolumesBWR();
}

int DetectorConstruction::ReadMCNPmatCard(const char *FileName)
{
	
	G4cout << "======================================================" << G4endl;
	G4cout << "== Reading material from : " << FileName << G4endl;
	G4cout << "======================================================" << G4endl;

	std::vector<std::string> ZAIDs;
	std::vector<int> ZZ;
	std::vector<int> AA;
	std::vector<double> WF;
	double TotMass = 0.;

	std::string str;
	//read AME-2016
	G4String dir(INPUT_DIR);
	G4String aFile = dir + "/atomic-mass-eval-2016.txt";
	std::ifstream textfile(aFile.data());
	if(!textfile.is_open()) {
      G4ExceptionDescription ed;
      ed << "      Could not open the file: input/input/atomic-mass-eval-2016.txt" << G4endl;
      G4Exception("DetectorConstruction::ReadMCNPmatCard(...)",
                  "DetectorConstruction.01",
                  FatalException,
                  ed);
    }
	getline(textfile,str);
	getline(textfile,str);

	std::vector<int> AME_Z;
	std::vector<int> AME_A;
	std::vector<double> AME_mass;
	while(!textfile.eof()) {
		int Z, N, A;
		double mass;
		textfile >> N >> Z >> A >> mass;
		AME_Z.push_back(Z);
		AME_A.push_back(A);
		AME_mass.push_back(mass*1.e-06);
	}
	textfile.close();

	//read MCNP material card
	textfile.open(FileName);
	if(!textfile.is_open()) {
      return 1;
    }

	
	do {
		getline(textfile,str);
		if(str.find("Zone mass (grams):")!=std::string::npos) {
			size_t pos = str.find("Zone mass (grams):");
			std::string sTotMass = str.substr(pos+18,str.size()-(pos+18));
			TotMass = atof(sTotMass.data());
		}
	} while(str[0]=='C');

	std::stringstream ss(str);
	std::string name;
	std::string ZAID;
	double WeightFrac;
	ss >> name >> ZAID >> WeightFrac;
	WeightFrac *= -1.;

	std::string sMassNumber = ZAID.substr(ZAID.size()-3,ZAID.size());
	std::string sAtomicNumber = ZAID.substr(0,ZAID.size()-3);

	int MassNumber = atoi(sMassNumber.data());
	int AtomicNumber = atoi(sAtomicNumber.data());

	ZZ.push_back(AtomicNumber);
	AA.push_back(MassNumber);
	WF.push_back(WeightFrac);
	ZAIDs.push_back(ZAID);

	G4int maxZ = 0;
	while(!textfile.eof()) {
		getline(textfile,str);
		if(str.find("nlib")!=std::string::npos) break;
		std::stringstream ss1(str);
		ss1 >> ZAID >> WeightFrac;
		WeightFrac *= -1.;

		sMassNumber = ZAID.substr(ZAID.size()-3,ZAID.size());
		sAtomicNumber = ZAID.substr(0,ZAID.size()-3);

		MassNumber = atoi(sMassNumber.data());
		AtomicNumber = atoi(sAtomicNumber.data());
		if(AtomicNumber>maxZ) maxZ = AtomicNumber;

		ZZ.push_back(AtomicNumber);
		AA.push_back(MassNumber);
		WF.push_back(WeightFrac);
		ZAIDs.push_back(ZAID);

		//G4cout << ZAID << "  " << AtomicNumber << "  " << MassNumber << "  " << WeightFrac << G4endl;
	}

	textfile.close();

//---------------------------------------------------------
	std::vector<G4Element*> theElements;
	std::vector<double> MassofElements;
	G4double TotMassMod = 0; //this is the sum of the masses after rejecting those elements with very small relative mass
	for(int Z=1;Z<=maxZ;++Z) {
		std::vector<G4Isotope*> Isotopes; //will contain all isotopes of this Z
		std::vector<double> moles; //will contain the abundacnce of this isotope
		double sum_of_mass = 0.;
		double sum_of_moles = 0.;

		G4bool found = false;
		for(size_t i=0;i<ZZ.size();i++) {
			if(Z!=ZZ[i]) continue;

			found = true;
			double atomic_mass = 0.;
			for(size_t k=0;k<AME_Z.size();k++) {
				if(ZZ[i]==AME_Z[k] && AA[i]==AME_A[k]) atomic_mass = AME_mass[k];
			}
			Isotopes.push_back(new G4Isotope(ZAIDs[i].data(),Z,AA[i],atomic_mass*g/mole));

			double mass = WF[i]*TotMass;
			moles.push_back(mass/atomic_mass);
			sum_of_mass += mass;
			sum_of_moles += mass/atomic_mass;
		}
		if(!found) continue;

		G4double massFraction = sum_of_mass/TotMass;

		if(massFraction>1.E-9) {
			TotMassMod += sum_of_mass;
			char ElName[32];
			snprintf(ElName,32,"Element-%d",Z);
			G4Element *Element = new G4Element(ElName,ElName,Isotopes.size());
			for(size_t j=0;j<Isotopes.size();j++) {
			Element->AddIsotope(Isotopes[j],moles[j]/sum_of_moles);
		}

		theElements.push_back(Element);
		MassofElements.push_back(sum_of_mass);
		}

		Isotopes.clear();
		moles.clear();
	}

	G4double PinRad = 0.5;
	G4double PinLentgh = 367.;

	double FuelVolume = 17.*17.*PinLentgh*pow(PinRad,2.)*3.141592654; //total volume of fuel in the w17x17 (17 x 17 x rod-height x pi x ror-rad²)
	double density = TotMass/FuelVolume*g/cm3;
	//double density = 9.105*g/cm3;
//	G4cout << "fuel material denisty is " << density/(g/cm3) << " g/cm3" << G4endl;
	fFuelMat = new G4Material("fuel",density,theElements.size());
	for(size_t i=0;i<theElements.size();i++) {
		G4double massFraction(MassofElements[i]/TotMassMod);
		fFuelMat->AddElement(theElements[i],massFraction);
	}

	return 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
	// Material definition

	/*
	//---------PWR-------------
	G4double PinRad = 0.5*cm;
	G4double PinLenght = 367.*cm;

	G4double FuelVolume = 17.*17.*PinLenght*pow(PinRad,2.)*3.141592654;

	FuelVolume = PinLenght*pow(22.15*cm,2.);
	MaterialFile *mFile = new MaterialFile("./input/SKB-TR-10-13-PWR-5years-nuclide-vector.plt","fuel-material",FuelVolume);
	fFuelMat = mFile->GetMaterial();
	*/

	//---------BWR-------------
	G4double FuelRodPitch = 12.4*mm;
	G4double FuelRodOuterDiameter = 9.62*mm;
	G4double FuelRodInnerDiameter = 8.36*mm;
	G4double FuelRodLength = 3712.*mm; // is acctully for SVEA-64 (I would guess it is the same)
	G4double CladdingThickness = 0.63*mm;

	G4double PinRad = 0.5*9.62*mm;
	G4double PinLenght = 3712.*mm;

	G4double FuelVolume = 96.*PinLenght*pow(PinRad,2.)*3.141592654;

	MaterialFile *mFile = new MaterialFile("./input/SKB-TR-10-13-BWR-CRAM-37years-nuclide-vector.plt","fuel-material",FuelVolume);
	fFuelMat = mFile->GetMaterial();
	#ifndef NOT_USING_MPI
	if(G4MPImanager::GetManager()->GetRank()==0) mFile->Print();
	#endif


	// Cast Iron
	G4double CastIronDensity = 7300.*kg/m3; //
	//G4double CastIronDensity = 7800.*kg/m3; // Debras number
	fCastIron = new G4Material("CastIron",CastIronDensity,9);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_Cu"),0.017*perCent);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_C"),3.41*perCent);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_Si"),2.33*perCent);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_Mn"),0.16*perCent);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_P"),0.038*perCent);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_S"),0.007*perCent);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ni"),0.52*perCent);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_Mg"),0.046*perCent);
	fCastIron->AddMaterial(nistManager->FindOrBuildMaterial("G4_Fe"),93.472*perCent);

	//fFuelMat = nistManager->FindOrBuildMaterial("G4_AIR");
	// Air defined using NIST Manager
	fWorldMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
	//fWorldMaterial = nistManager->FindOrBuildMaterial("G4_Galactic");
	fVacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
	fAlu = nistManager->FindOrBuildMaterial("G4_Al");

	// LS-309
	G4Element *fCarbon = nistManager->FindOrBuildElement("C");
	G4Element *fHydrogen = nistManager->FindOrBuildElement("H");
	fEJ309 = new G4Material("LS309",0.959*g/cm3,2);
	fEJ309->AddElement(fHydrogen,9.48*perCent);
	fEJ309->AddElement(fCarbon,90.52*perCent);

	// EJ-276
	// the electron density is slightly wrong 3.8 % too low (33.982x10^22 vs 35.33x10^22)
	// when using only hydrogen and carbon (this should lower the gamma-ray efficiency slightly)
	// the reason for this is probably the solvent which is not included in the specs.
	// I will assume the rest comes from Oxygen (Carbon oxides are mentioned in the Safety data sheet,
	// allthough no numbers are given)
	G4Element *fOxygen = nistManager->FindOrBuildElement("O");
	fEJ276 = new G4Material("EJ276",1.096*g/cm3,3, kStateSolid);
	G4double nC = 4.906; // No. of Carbon atoms per cm3 x10^22
	G4double nH = 4.546; // No. of Hydrogen atoms per cm3 x10^22
	G4double nO = 0.1685; // No. of Oxygen atoms per cm3 x10^22 (calculated from the missing electron density)
	fEJ276->AddElement(fHydrogen,nH/(nH+nC+nO));
	fEJ276->AddElement(fCarbon,nC/(nH+nC+nO));
	fEJ276->AddElement(fOxygen,nO/(nH+nC+nO));

	fCladingMat = nistManager->FindOrBuildMaterial("G4_Zr"); //assuming pure Zirkonium

	fBoroSilicate = new G4Material("BoroSilicate",2.23*g/cm3,5);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_Si"),39.209*perCent);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_B"),2.918*perCent);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_Na"),3.182*perCent);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_Al"),1.288*perCent);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_O"),53.403*perCent);

	G4Material *SWX201 = new G4Material("SWX201",0.95*g/cm3,2);
	SWX201->AddMaterial(nistManager->FindOrBuildMaterial("G4_POLYETHYLENE"),95.*perCent);
	SWX201->AddMaterial(nistManager->FindOrBuildMaterial("G4_B"),5.*perCent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    G4Material* air  = G4Material::GetMaterial("G4_AIR");
	
	//world
    G4double worldRadiusY = 250*cm;
    G4double worldRadiusX = 250*cm;
	G4double worldLength = 300*cm;
    //********************Definitions of Solids, Logical Volumes, Physical Volumes***************************
    
    // World
    
    G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

    G4Box* worldS = new G4Box("world",worldRadiusX,worldRadiusY,worldLength); 
    G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                          worldS,   //its solid
                          air,      //its material
                          "World.logical"); //its name
    
    //  Must place the World Physical volume unrotated at (0,0,0).
    // 
    G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        worldLV,         // its logical volume
                        "World.physical",         // its name
                        0,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 
    //worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
	G4VisAttributes *WorldVisAtt = new G4VisAttributes(G4Colour::Gray());
	WorldVisAtt->SetForceWireframe(true);
  	worldLV->SetVisAttributes(WorldVisAtt);
	//worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

//============= The Copper Canister ==========================================
  	//Dimensions from Fig 3-7. & Tab 3-6 copper enclosure
  	G4double CuDimA = 4835.*mm;
  	G4double CuDimB = 1050.*mm;
  	G4double CuDimC = 850.*mm;  
  	G4double CuDimE = 952.*mm;
  	G4double CuDimG = CuDimC;
  	G4double CuDimK = 35.*mm;
  	G4double CuDimL = 50.*mm;
  	G4double CuDimM = 50.*mm;
  	G4double CuDimP = 75.*mm;
  	G4double CuDimF = 821.*mm;
  	G4double InnerFreeLength = 4575.*mm;

  	//Dimensions from Fig 3-6. & Tab 3-3,3-4,3-5
  	G4double InsertDimA = 4573*mm;
  	G4double InsertDimB = 80.*mm; //PWR
  	G4double InsertDimC = 4443.*mm; //PWR
  	//G4double InsertDimB = 60.*mm; //BWR
  	//G4double InsertDimC = 4463.*mm; //BWR
  	G4double InsertDimD = 949.*mm;
  	G4double InsertDimE = 910.*mm;
  	G4double InsertDimF = 50.*mm;
  	G4double InsertDimG = 5.*deg;

  	G4Tubs *CanS = new G4Tubs("CanS",0.,CuDimB/2.,CuDimA/2.,0.,360.*deg);
  	G4Tubs *CanTopCut1S = new G4Tubs("CanTopCut1S",0.,CuDimG/2.,CuDimL/2.,0.,360.*deg);
  	G4Tubs *CanTopCut2S = new G4Tubs("CanTopCut2S",0.,CuDimF/2.,CuDimK/2.,0.,360.*deg);
  	G4Tubs *CanBottomCut1S = new G4Tubs("CanTopCut2S",0.,CuDimC/2.,CuDimP/2.,0.,360.*deg);

	G4Material* CopperMaterial = nistManager->FindOrBuildMaterial("G4_Cu");
	G4Material* ArgonMaterial = nistManager->FindOrBuildMaterial("G4_Ar");
	G4Material* StainLessSteelMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

	G4LogicalVolume *CanLV = new G4LogicalVolume(CanS,CopperMaterial,"CaskCylinderLV",0,0,0);
	G4LogicalVolume *CanTopCut1LV = new G4LogicalVolume(CanTopCut1S,air,"CanTopCut1LV",0,0,0);
	G4LogicalVolume *CanTopCut2LV = new G4LogicalVolume(CanTopCut2S,air,"CanTopCut2LV",0,0,0);
	G4LogicalVolume *CanBottomCut1LV = new G4LogicalVolume(CanBottomCut1S,air,"CanBottomCut1LV",0,0,0);

	CanTopCut1LV->SetVisAttributes(G4VisAttributes::GetInvisible());
	CanTopCut2LV->SetVisAttributes(G4VisAttributes::GetInvisible());
	CanBottomCut1LV->SetVisAttributes(G4VisAttributes::GetInvisible());

	//  interior
  	G4Tubs *CavityS = new G4Tubs("CavityS",0.,CuDimE/2.,InnerFreeLength/2.,0.,360.*deg);

  	G4double pRmax2 = InsertDimE/2.;
  	G4double pRmax1 = pRmax2 - InsertDimF*std::tan(InsertDimG);
  	G4Cons *SteelLidS = new G4Cons("SteelLidS",0.,pRmax1,0.,pRmax2,InsertDimF/2.,0.,360.*deg);

  	G4Tubs *InsertS = new G4Tubs("InsertS",0.,InsertDimD/2.,InsertDimA/2.,0.,360.*deg);


	G4LogicalVolume *CavityLV = new G4LogicalVolume(CavityS,ArgonMaterial,"CaskCavityLV",0,0,0); //
	G4LogicalVolume *InsertLV = new G4LogicalVolume(InsertS,fCastIron,"InsertLV",0,0,0); //
	G4LogicalVolume *SteelLidLV = new G4LogicalVolume(SteelLidS,StainLessSteelMaterial,"SteelLidLV",0,0,0); //this should be  iron in the end

	G4VisAttributes *GreenVisAtt = new G4VisAttributes(G4Colour::Green());
	CavityLV->SetVisAttributes(GreenVisAtt);

	G4VisAttributes *GrayVisAtt = new G4VisAttributes(G4Colour::Gray());
	InsertLV->SetVisAttributes(GrayVisAtt);

	G4VisAttributes *YellowVisAtt = new G4VisAttributes(G4Colour::Yellow());
	CanLV->SetVisAttributes(YellowVisAtt);

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        CanLV,         // its logical volume
                        "CaskCylinderPV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,CuDimA/2.- CuDimK - CuDimL/2.), // at (0,0,0)
                        CanTopCut1LV,         // its logical volume
                        "CanTopCut1PV",         // its name
                        CanLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,CuDimA/2.- CuDimK/2.), // at (0,0,0)
                        CanTopCut2LV,         // its logical volume
                        "CanTopCut2PV",         // its name
                        CanLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,-CuDimA/2. + CuDimP/2.), // at (0,0,0)
                        CanBottomCut1LV,         // its logical volume
                        "CanBottomCut1PV",         // its name
                        CanLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,+CuDimA/2.-InnerFreeLength/2.-CuDimK-CuDimL-CuDimM), // at (0,0,0)
                        CavityLV,         // its logical volume
                        "CaskCavityPV",         // its name
                        CanLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        InsertLV,         // its logical volume
                        "InsertPV",         // its name
                        CavityLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,+InsertDimA/2.-InsertDimF/2.), // at (0,0,0)
                        SteelLidLV,         // its logical volume
                        "SteelLidPV",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

//-------- Insert Channels -----------------------------------------

  	//Dimensions from Fig 3-6. & Tab 3-3,3-4,3-5
  	G4double InsertDimI = 20.*mm;
  	G4double InsertDimJ = 370.*mm;
  	G4double InsertDimK = 110.*mm;
  	G4double InsertDimL = 235.*mm;
  	G4double InsertDimM = 12.5*mm;

  	G4double ChannelLength = InsertDimA - InsertDimF - InsertDimB;

  	G4Box *ChannelOutS = new G4Box("ChannelOutS",InsertDimL/2.+InsertDimM,InsertDimL/2.+InsertDimM,InsertDimC/2.);
  	G4Box *ChannelInS = new G4Box("ChannelInS",InsertDimL/2.,InsertDimL/2.,InsertDimC/2.);

  	G4LogicalVolume *ChannelOutLV = new G4LogicalVolume(ChannelOutS,StainLessSteelMaterial,"ChannelOutLV",0,0,0);
  	G4LogicalVolume *ChannelInLV = new G4LogicalVolume(ChannelInS,ArgonMaterial,"ChannelInLV",0,0,0);

	ChannelOutLV->SetVisAttributes(GrayVisAtt);
	ChannelInLV->SetVisAttributes(GreenVisAtt);

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        ChannelInLV,         // its logical volume
                        "ChannelInPV",         // its name
                        ChannelOutLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

  	G4double deltaZ = InsertDimA/2. - InsertDimF - InsertDimC/2.;

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(InsertDimJ/2.,InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV1",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-InsertDimJ/2.,InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV2",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-InsertDimJ/2.,-InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV3",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(InsertDimJ/2.,-InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV4",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	// Create a region
	G4Region* FuelRegion = new G4Region("FuelRegion");
	FuelRegion->AddRootLogicalVolume(InsertLV);
	//FuelRegion->AddRootLogicalVolume(CaskBottomLV);
	//FuelRegion->AddRootLogicalVolume(CaskTopLV);
//******************FUEL***********************************************************
/*
	G4double PinRad = 5.*mm;
	G4double FuelRad = 3.8*mm;

	G4double HalfWidthBox = 221.5/17./2.;

	G4Box *DummyFuelBox = new G4Box("DummyFuelBox",HalfWidthBox,HalfWidthBox,3670./2.*mm);
	//G4Box *DummyFuelBox = new G4Box("DummyFuelBox",HalfWidthBox,HalfWidthBox,HalfWidthBox);
	G4LogicalVolume *DummyFuelBoxLV = new G4LogicalVolume(DummyFuelBox,ArgonMaterial,"DummyFuelBoxLV",0,0,0);
	DummyFuelBoxLV->SetVisAttributes(G4VisAttributes::GetInvisible());

	G4Box *DummyFuelBox17 = new G4Box("DummyFuelBox",17.*HalfWidthBox,HalfWidthBox,3670./2.*mm);
	//G4Box *DummyFuelBox17 = new G4Box("DummyFuelBox17",17*HalfWidthBox,HalfWidthBox,HalfWidthBox);
	G4LogicalVolume *DummyFuelBox17LV = new G4LogicalVolume(DummyFuelBox17,ArgonMaterial,"DummyFuelBox17LV",0,0,0);
	DummyFuelBox17LV->SetVisAttributes(G4VisAttributes::GetInvisible());

	G4Box *DummyFuelBox17x17 = new G4Box("DummyFuelBox",17*HalfWidthBox,17*HalfWidthBox,3670./2.*mm);
	G4LogicalVolume *DummyFuelBox17x17LV = new G4LogicalVolume(DummyFuelBox17x17,ArgonMaterial,"DummyFuelBox17x17LV",0,0,0);
	DummyFuelBox17x17LV->SetVisAttributes(G4VisAttributes::GetInvisible());

	//Since the ORIGEN calculation does not sepparate fuel and cladding I can't do it either
	G4Tubs *CladdingRodS = new G4Tubs("CladdingRodS",0.,PinRad,3670./2.*mm,0.,360.*deg);
	G4LogicalVolume *CladdingRodLV = new G4LogicalVolume(CladdingRodS,fFuelMat,"CladdingRodLV",0,0,0);
	//G4LogicalVolume *CladdingRodLV = new G4LogicalVolume(CladdingRodS,air,"CladdingRodLV",0,0,0);


	G4VisAttributes *RedVisAtt = new G4VisAttributes(G4Colour::Red());
	CladdingRodLV->SetVisAttributes(RedVisAtt);

	//place pin in dummy box
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        CladdingRodLV,         // its logical volume
                        "FuelRodPV",         // its name
                        DummyFuelBoxLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps


	//G4VPhysicalVolume *repX = 
	new G4PVReplica("FuelAssembly17",
                 DummyFuelBoxLV,
                 DummyFuelBox17LV,
                 kXAxis, 17, 2.*HalfWidthBox);
	
	//G4VPhysicalVolume *FuelAssembly =
	new G4PVReplica("FuelAssembly17x17",
                 DummyFuelBox17LV,
                 DummyFuelBox17x17LV,
                 kYAxis, 17, 2.*HalfWidthBox);

	//place 1 assembly
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        DummyFuelBox17x17LV,         // its logical volume
                        "Fuel",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps
*/
	// fuel modeled as a simple homogenous Volume

	G4double HalfWidthBox = 214./2.;
	G4Box *DummyFuelBox = new G4Box("DummyFuelBox",HalfWidthBox,HalfWidthBox,4296./2.*mm);
	G4LogicalVolume *DummyFuelBoxLV = new G4LogicalVolume(DummyFuelBox,fFuelMat,"DummyFuelBoxLV",0,0,0);
	//place 1 assembly
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        DummyFuelBoxLV,         // its logical volume
                        "Fuel",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps

//****************Detectors********************************************************

	const G4double PbWallThickness = 4.*cm;
	const G4double BoronWallThickness = 2.*mm;
	//const G4double CuWallThickness = 0.*mm;
	const G4double DummyDetectionVolume = 1.*mm;
	const G4double totWallThickness = PbWallThickness + BoronWallThickness + DummyDetectionVolume;

	//const G4double ShieldingWallDistance = CuDimB/2. + totWallThickness/2. + 5.*cm;
	//const G4double distance = CuDimB/2. + 16.3*cm; // lead wall at surface of canister

	const G4double distance = CuDimB/2. + 47.5*cm; // 1 m from center
	const G4double ShieldingWallDistance = distance - totWallThickness/2. - 5.*cm;
	//const G4double ShieldingWallDistance = CuDimB/2. + 47.5*cm/2.;

	G4double WallHalfY = InsertDimJ/2. + 90.*mm;
	G4double WallHalfZ = WallHalfY;

	G4Material *fPb = nistManager->FindOrBuildMaterial("G4_Pb");
	//G4Material *fPb = nistManager->FindOrBuildMaterial("G4_AIR");
	G4Box *ShieldingWallS = new G4Box("ShieldingWallS",totWallThickness/2.,WallHalfY,WallHalfZ);
  	G4LogicalVolume *ShieldingWallLV = new G4LogicalVolume(ShieldingWallS,fPb,"ShieldingWallLV",0,0,0);

	G4Box *B4CWallS = new G4Box("B4CWallS",BoronWallThickness/2.,WallHalfY,WallHalfZ);
	G4Material *fB4C = nistManager->FindOrBuildMaterial("G4_BORON_CARBIDE");
	//G4Material *fB4C = nistManager->FindOrBuildMaterial("G4_AIR");
  	G4LogicalVolume *B4CWallLV = new G4LogicalVolume(B4CWallS,fB4C,"B4CWallLV",0,0,0);

  	G4double PbBoxWallThicknes = 5.*cm;
  	G4Box *ShieldingBoxWallS = new G4Box("ShieldingBoxWallS",WallHalfY+PbWallThickness/2.,WallHalfY+2.*PbBoxWallThicknes,WallHalfY+2.*PbBoxWallThicknes);
  	G4Box *ShieldingBoxHoleS = new G4Box("ShieldingBoxHoleS",WallHalfY+PbWallThickness,WallHalfY,WallHalfY);
  	G4SubtractionSolid *ShieldingBoxS = new G4SubtractionSolid("ShieldingBoxS",ShieldingBoxWallS,ShieldingBoxHoleS);
  	G4LogicalVolume *ShieldingBoxLV = new G4LogicalVolume(ShieldingBoxS,fPb,"ShieldingBoxLV",0,0,0);

  	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-totWallThickness/2.+BoronWallThickness/2.,0,0), // at (0,0,0)
                        B4CWallLV,         // its logical volume
                        "B4CWallPV",         // its name
                        ShieldingWallLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps
/*
  	if(CuWallThickness>1E-3*mm) {
  		G4Box *CuWallS = new G4Box("CuWallS",CuWallThickness/2.,InsertDimJ,2.*InsertDimJ);
  		G4LogicalVolume *CuWallLV = new G4LogicalVolume(CuWallS,CopperMaterial,"CuWallLV",0,0,0);
  		CuWallLV->SetVisAttributes(YellowVisAtt);

  		new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(totWallThickness/2.-CuWallThickness/2.,0,0), // at (0,0,0)
                        CuWallLV,         // its logical volume
                        "CuWallPV",         // its name
                        ShieldingWallLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps
  	}*/

  	if(DummyDetectionVolume>1E-3*mm) {
  		G4Box *DummyS = new G4Box("DummyS",DummyDetectionVolume/2.,WallHalfY,WallHalfZ);
  		G4LogicalVolume *DummyLV = new G4LogicalVolume(DummyS,air,"DetectionVolume",0,0,0);

  		new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(totWallThickness/2.-DummyDetectionVolume/2.,0,0), // at (0,0,0)
                        DummyLV,         // its logical volume
                        "DetectionVolumePV",         // its name
                        ShieldingWallLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps
  	}

  	for(G4int i=0;i<4;i++) {
  		G4RotationMatrix *rot = new G4RotationMatrix();
  		rot->rotateZ(i*90*deg);
  		G4ThreeVector pos(ShieldingWallDistance,0,0);
  		pos *= *rot;

  		G4ThreeVector pos1(ShieldingWallDistance+WallHalfY,0,0);
  		pos1 *= *rot;

  		rot->rotateX(45.*deg);
  		new G4PVPlacement(
                        rot,               // rotation
                        pos, // at (0,0,0)
                        ShieldingWallLV,         // its logical volume
                        "ShieldingWallPV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

  		new G4PVPlacement(
                        rot,               // rotation
                        pos1, // at (0,0,0)
                        ShieldingBoxLV,         // its logical volume
                        "ShieldingBoxPV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps
  		
  	}
	/*new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(ShieldingWallDistance,0,0), // at (0,0,0)
                        ShieldingWallLV,         // its logical volume
                        "ShieldingWallPV1",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps*/


	const G4int nDetClusters = 4;
	const G4int nDetPerCluster = 4;
	const G4int nDetectors = nDetPerCluster*nDetClusters;
	G4AssemblyVolume* neutronDetectors[nDetectors];
	char DetName[32];
	for(G4int i=0;i<nDetectors;i++) {
		snprintf(DetName,32,"EJ309_%d",i+1);
		neutronDetectors[i] = EJ309_5x5inch(i+1, DetName);
	}

	//G4double phiCluster[] = {0.,180.};
	G4double phiCluster[] = {0.,90.,180.,270.};
	for(G4int cluster=0;cluster<nDetClusters;cluster++) {
		G4ThreeVector position[nDetPerCluster];
		position[0] = G4ThreeVector(0.,InsertDimJ/2.,0.);
		position[1] = G4ThreeVector(InsertDimJ/2.,0.,0.);
		position[2] = G4ThreeVector(-InsertDimJ/2.,0.,0.);
		position[3] = G4ThreeVector(0.,-InsertDimJ/2.,0.);

		G4ThreeVector posCluster = G4ThreeVector(0.,0.,distance);
		G4RotationMatrix *rotPos = new G4RotationMatrix();
		rotPos->rotateY(90.*degree);
		rotPos->rotateZ(phiCluster[cluster]*degree);

		G4RotationMatrix *rotDet = new G4RotationMatrix();
		rotDet->rotateY(90.*degree);
		rotDet->rotateZ(phiCluster[cluster]*degree);

		//posCluster *= *rotPos;
		for(G4int i=0;i<nDetPerCluster;i++) {
			position[i] += posCluster;
			position[i] *= *rotPos;
			neutronDetectors[cluster*nDetPerCluster + i]->MakeImprint(worldLV,position[i],rotDet);
		}
	}

	// Create a region
	G4LogicalVolume* logicVol = G4LogicalVolumeStore::GetInstance()->GetVolume("5inch-EJ309");
	G4Region* DetectorRegion = new G4Region("DetectorRegion");
	DetectorRegion->AddRootLogicalVolume(logicVol);
	
	// Print materials

	#ifndef NOT_USING_MPI
	if(G4MPImanager::GetManager()->GetRank()==0) G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	#else
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	#endif
/*
	// test 3x3 inch
	G4AssemblyVolume* test = EJ309_3x3inch(0,"testDet");
	G4RotationMatrix *rotTest = new G4RotationMatrix();
	G4ThreeVector posTest(0,0,2.7*m);
	test->MakeImprint(worldLV,posTest,rotTest);*/

	return worldPV;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::DefineVolumesBWR()
{
    G4Material* air  = G4Material::GetMaterial("G4_AIR");
	
	//world
    G4double worldRadiusY = 250*cm;
    G4double worldRadiusX = 250*cm;
	G4double worldLength = 300*cm;
    //********************Definitions of Solids, Logical Volumes, Physical Volumes***************************
    
    // World
    
    G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

    G4Box* worldS = new G4Box("world",worldRadiusX,worldRadiusY,worldLength); 
    G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                          worldS,   //its solid
                          air,      //its material
                          "World.logical"); //its name
    
    //  Must place the World Physical volume unrotated at (0,0,0).
    // 
    G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        worldLV,         // its logical volume
                        "World.physical",         // its name
                        0,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 
    //worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
	G4VisAttributes *WorldVisAtt = new G4VisAttributes(G4Colour::Gray());
	WorldVisAtt->SetForceWireframe(true);
  	worldLV->SetVisAttributes(WorldVisAtt);
	//worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

//============= The Copper Canister ==========================================
  	//Dimensions from Fig 3-7. & Tab 3-6 copper enclosure
  	G4double CuDimA = 4835.*mm; // total length
  	G4double CuDimB = 1050.*mm; // outer diameter
  	G4double CuDimC = 850.*mm;  
  	G4double CuDimE = 952.*mm;
  	G4double CuDimG = CuDimC;
  	G4double CuDimK = 35.*mm;
  	G4double CuDimL = 50.*mm;
  	G4double CuDimM = 50.*mm;
  	G4double CuDimP = 75.*mm;
  	G4double CuDimF = 821.*mm;
  	G4double InnerFreeLength = 4575.*mm;

  	//Dimensions from Fig 3-6. & Tab 3-3,3-4,3-5
  	G4double InsertDimA = 4573*mm;
  	//G4double InsertDimB = 80.*mm; //PWR
  	//G4double InsertDimC = 4443.*mm; //PWR
  	G4double InsertDimB = 60.*mm; //BWR
  	G4double InsertDimC = 4463.*mm; //BWR
  	G4double InsertDimD = 949.*mm;
  	G4double InsertDimE = 910.*mm;
  	G4double InsertDimF = 50.*mm;
  	G4double InsertDimG = 5.*deg;

  	G4Tubs *CanS = new G4Tubs("CanS",0.,CuDimB/2.,CuDimA/2.,0.,360.*deg);
  	G4Tubs *CanTopCut1S = new G4Tubs("CanTopCut1S",0.,CuDimG/2.,CuDimL/2.,0.,360.*deg);
  	G4Tubs *CanTopCut2S = new G4Tubs("CanTopCut2S",0.,CuDimF/2.,CuDimK/2.,0.,360.*deg);
  	G4Tubs *CanBottomCut1S = new G4Tubs("CanTopCut2S",0.,CuDimC/2.,CuDimP/2.,0.,360.*deg);

	G4Material* CopperMaterial = nistManager->FindOrBuildMaterial("G4_Cu");
	G4Material* ArgonMaterial = nistManager->FindOrBuildMaterial("G4_Ar");
	G4Material* StainLessSteelMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

	G4LogicalVolume *CanLV = new G4LogicalVolume(CanS,CopperMaterial,"CaskCylinderLV",0,0,0);
	G4LogicalVolume *CanTopCut1LV = new G4LogicalVolume(CanTopCut1S,air,"CanTopCut1LV",0,0,0);
	G4LogicalVolume *CanTopCut2LV = new G4LogicalVolume(CanTopCut2S,air,"CanTopCut2LV",0,0,0);
	G4LogicalVolume *CanBottomCut1LV = new G4LogicalVolume(CanBottomCut1S,air,"CanBottomCut1LV",0,0,0);

	CanTopCut1LV->SetVisAttributes(G4VisAttributes::GetInvisible());
	CanTopCut2LV->SetVisAttributes(G4VisAttributes::GetInvisible());
	CanBottomCut1LV->SetVisAttributes(G4VisAttributes::GetInvisible());

	//  interior
  	G4Tubs *CavityS = new G4Tubs("CavityS",0.,CuDimE/2.,InnerFreeLength/2.,0.,360.*deg);

  	G4double pRmax2 = InsertDimE/2.;
  	G4double pRmax1 = pRmax2 - InsertDimF*std::tan(InsertDimG);
  	G4Cons *SteelLidS = new G4Cons("SteelLidS",0.,pRmax1,0.,pRmax2,InsertDimF/2.,0.,360.*deg);

  	G4Tubs *InsertS = new G4Tubs("InsertS",0.,InsertDimD/2.,InsertDimA/2.,0.,360.*deg);


	G4LogicalVolume *CavityLV = new G4LogicalVolume(CavityS,ArgonMaterial,"CaskCavityLV",0,0,0); //
	G4LogicalVolume *InsertLV = new G4LogicalVolume(InsertS,fCastIron,"InsertLV",0,0,0); //
	G4LogicalVolume *SteelLidLV = new G4LogicalVolume(SteelLidS,StainLessSteelMaterial,"SteelLidLV",0,0,0); //this should be  iron in the end

	G4VisAttributes *GreenVisAtt = new G4VisAttributes(G4Colour::Green());
	CavityLV->SetVisAttributes(GreenVisAtt);

	G4VisAttributes *GrayVisAtt = new G4VisAttributes(G4Colour::Gray());
	InsertLV->SetVisAttributes(GrayVisAtt);

	G4VisAttributes *YellowVisAtt = new G4VisAttributes(G4Colour::Yellow());
	CanLV->SetVisAttributes(YellowVisAtt);

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        CanLV,         // its logical volume
                        "CaskCylinderPV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,CuDimA/2.- CuDimK - CuDimL/2.), // at (0,0,0)
                        CanTopCut1LV,         // its logical volume
                        "CanTopCut1PV",         // its name
                        CanLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,CuDimA/2.- CuDimK/2.), // at (0,0,0)
                        CanTopCut2LV,         // its logical volume
                        "CanTopCut2PV",         // its name
                        CanLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,-CuDimA/2. + CuDimP/2.), // at (0,0,0)
                        CanBottomCut1LV,         // its logical volume
                        "CanBottomCut1PV",         // its name
                        CanLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,+CuDimA/2.-InnerFreeLength/2.-CuDimK-CuDimL-CuDimM), // at (0,0,0)
                        CavityLV,         // its logical volume
                        "CaskCavityPV",         // its name
                        CanLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        InsertLV,         // its logical volume
                        "InsertPV",         // its name
                        CavityLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,+InsertDimA/2.-InsertDimF/2.), // at (0,0,0)
                        SteelLidLV,         // its logical volume
                        "SteelLidPV",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

//-------- Insert Channels -----------------------------------------

  	//Dimensions from Fig 3-6. & Tab 3-3,3-4,3-5
  	G4double InsertDimI = 20.*mm; //BWR
  	G4double InsertDimJ = 210.*mm; //BWR
  	G4double InsertDimK = 30.*mm; //BWR
  	G4double InsertDimL = 160.*mm; //BWR
  	G4double InsertDimM = 10.*mm; //BWR

  	G4double ChannelLength = InsertDimA - InsertDimF - InsertDimB;

  	G4Box *ChannelOutS = new G4Box("ChannelOutS",InsertDimL/2.+InsertDimM,InsertDimL/2.+InsertDimM,InsertDimC/2.);
  	G4Box *ChannelInS = new G4Box("ChannelInS",InsertDimL/2.,InsertDimL/2.,InsertDimC/2.);

  	G4LogicalVolume *ChannelOutLV = new G4LogicalVolume(ChannelOutS,StainLessSteelMaterial,"ChannelOutLV",0,0,0);
  	G4LogicalVolume *ChannelInLV = new G4LogicalVolume(ChannelInS,ArgonMaterial,"ChannelInLV",0,0,0);

	ChannelOutLV->SetVisAttributes(GrayVisAtt);
	ChannelInLV->SetVisAttributes(GreenVisAtt);

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        ChannelInLV,         // its logical volume
                        "ChannelInPV",         // its name
                        ChannelOutLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

  	G4double deltaZ = InsertDimA/2. - InsertDimF - InsertDimC/2.;

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-InsertDimJ/2.,InsertDimJ*3./2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV1",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(InsertDimJ/2.,InsertDimJ*3./2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV2",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-InsertDimJ*3./2.,InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV3",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-InsertDimJ/2.,InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV4",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(InsertDimJ/2.,InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV5",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(InsertDimJ*3./2.,InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV6",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-InsertDimJ*3./2.,-InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV7",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-InsertDimJ/2.,-InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV8",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(InsertDimJ/2.,-InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV9",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(InsertDimJ*3./2.,-InsertDimJ/2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV10",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-InsertDimJ/2.,-InsertDimJ*3./2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV11",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(InsertDimJ/2.,-InsertDimJ*3./2.,deltaZ), // at (0,0,0)
                        ChannelOutLV,         // its logical volume
                        "ChannelOutPV12",         // its name
                        InsertLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps

	// Create a region
	G4Region* FuelRegion = new G4Region("FuelRegion");
	FuelRegion->AddRootLogicalVolume(InsertLV);

	//FuelRegion->AddRootLogicalVolume(CaskBottomLV);
	//FuelRegion->AddRootLogicalVolume(CaskTopLV);
//******************FUEL***********************************************************
	// THe BWR assemblies are 10x10 we assume the SVEA-96
	// parameters taken from SKB TR10-13 Table A-3. /SKBdoc 1193244/
	G4double FuelRodPitch = 12.4*mm;
	G4double FuelRodOuterDiameter = 9.62*mm;
	G4double FuelRodInnerDiameter = 8.36*mm;
	G4double FuelRodLength = 3712.*mm; // is acctully for SVEA-64 (I would guess it is the same)
	G4double CladdingThickness = 0.63*mm;



	G4double PinRad = FuelRodOuterDiameter*0.5;

	G4double HalfWidthBox = FuelRodPitch*0.5;

	G4Box *DummyFuelBox = new G4Box("DummyFuelBox",HalfWidthBox,HalfWidthBox,0.5*FuelRodLength);
	//G4Box *DummyFuelBox = new G4Box("DummyFuelBox",HalfWidthBox,HalfWidthBox,HalfWidthBox.);
	G4LogicalVolume *DummyFuelBoxLV = new G4LogicalVolume(DummyFuelBox,ArgonMaterial,"DummyFuelBoxLV",0,0,0);
	DummyFuelBoxLV->SetVisAttributes(G4VisAttributes::GetInvisible());

	G4Box *DummyFuelBox5 = new G4Box("DummyFuelBox5",5.*HalfWidthBox,HalfWidthBox,0.5*FuelRodLength);
	G4LogicalVolume *DummyFuelBox5LV = new G4LogicalVolume(DummyFuelBox5,ArgonMaterial,"DummyFuelBox5LV",0,0,0);
	DummyFuelBox5LV->SetVisAttributes(G4VisAttributes::GetInvisible());

	G4Box *DummyFuelBox4 = new G4Box("DummyFuelBox4",4.*HalfWidthBox,HalfWidthBox,0.5*FuelRodLength);
	G4LogicalVolume *DummyFuelBox4LV = new G4LogicalVolume(DummyFuelBox4,ArgonMaterial,"DummyFuelBox4LV",0,0,0);
	DummyFuelBox4LV->SetVisAttributes(G4VisAttributes::GetInvisible());

	G4Box *DummyFuelBox5x4 = new G4Box("DummyFuelBox5x4",5*HalfWidthBox,4*HalfWidthBox,0.5*FuelRodLength);
	G4LogicalVolume *DummyFuelBox5x4LV = new G4LogicalVolume(DummyFuelBox5x4,ArgonMaterial,"DummyFuelBox5x4LV",0,0,0);
	DummyFuelBox5x4LV->SetVisAttributes(G4VisAttributes::GetInvisible());

	//Since the ORIGEN calculation does not sepparate fuel and cladding I can't do it either
	G4Tubs *CladdingRodS = new G4Tubs("CladdingRodS",0.,PinRad,0.5*FuelRodLength,0.,360.*deg);
	G4LogicalVolume *CladdingRodLV = new G4LogicalVolume(CladdingRodS,fFuelMat,"CladdingRodLV",0,0,0);
	//G4LogicalVolume *CladdingRodLV = new G4LogicalVolume(CladdingRodS,air,"CladdingRodLV",0,0,0);


	G4VisAttributes *RedVisAtt = new G4VisAttributes(G4Colour::Red());
	CladdingRodLV->SetVisAttributes(RedVisAtt);

	//place pin in dummy box
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        CladdingRodLV,         // its logical volume
                        "FuelRodPV",         // its name
                        DummyFuelBoxLV,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps


	//G4VPhysicalVolume *repX = 
	new G4PVReplica("FuelAssembly4",
                 DummyFuelBoxLV,
                 DummyFuelBox4LV,
                 kXAxis, 4, 2.*HalfWidthBox);

	new G4PVReplica("FuelAssembly5",
                 DummyFuelBoxLV,
                 DummyFuelBox5LV,
                 kXAxis, 5, 2.*HalfWidthBox);
	
	//G4VPhysicalVolume *FuelAssembly =
	new G4PVReplica("FuelAssembly5x4",
                 DummyFuelBox5LV,
                 DummyFuelBox5x4LV,
                 kYAxis, 4, 2.*HalfWidthBox);

	//place 1 assembly
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(5.*HalfWidthBox,6*HalfWidthBox,0), // at (0,0,0)
                        DummyFuelBox5x4LV,         // its logical volume
                        "Fuel5x4_A",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps

	//place 1 assembly
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-5.*HalfWidthBox,6*HalfWidthBox,0), // at (0,0,0)
                        DummyFuelBox5x4LV,         // its logical volume
                        "Fuel5x4_B",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps

	//place 1 assembly
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-5.*HalfWidthBox,-6*HalfWidthBox,0), // at (0,0,0)
                        DummyFuelBox5x4LV,         // its logical volume
                        "Fuel5x4_B",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps

	//place 1 assembly
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(5.*HalfWidthBox,-6*HalfWidthBox,0), // at (0,0,0)
                        DummyFuelBox5x4LV,         // its logical volume
                        "Fuel5x4_D",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps


	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(6.*HalfWidthBox,HalfWidthBox,0), // at (0,0,0)
                        DummyFuelBox4LV,         // its logical volume
                        "Fuel4x1_A",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-6.*HalfWidthBox,HalfWidthBox,0), // at (0,0,0)
                        DummyFuelBox4LV,         // its logical volume
                        "Fuel4x1_B",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(-6.*HalfWidthBox,-HalfWidthBox,0), // at (0,0,0)
                        DummyFuelBox4LV,         // its logical volume
                        "Fuel4x1_C",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps
	new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(6.*HalfWidthBox,-HalfWidthBox,0), // at (0,0,0)
                        DummyFuelBox4LV,         // its logical volume
                        "Fuel4x1_D",         // its name
                        ChannelInLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps

	G4double maxTime = 10000.*ns;
	InsertLV->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,maxTime));
	// fuel modeled as a simple homogenous Volume

  	/*const G4double distance = 600.*mm;
  	const G4int nRings = 16;
  	const G4int DetectorsPerRing = 18;
  	const G4int nDetectors = DetectorsPerRing*nRings;
  	G4AssemblyVolume* neutronDetectors[nDetectors];
	char DetName[32];

	G4double x_offset = twopi/DetectorsPerRing*distance - 2.*deltaZ;
	G4double x_0 = -0.5*x_offset*nRings;

	for(G4int ring=0;ring<nRings;++ring) {
		for(G4int i=ring*DetectorsPerRing;i<(ring+1)*DetectorsPerRing;i++) {
			snprintf(DetName,32,"EJ309_%d",i+1);
			neutronDetectors[i] = EJ309_1x2inch(i+1, DetName);

			G4double phiDet = (G4double(i-DetectorsPerRing)/G4double(DetectorsPerRing))*360.*degree;

			G4RotationMatrix *rotDet = new G4RotationMatrix();
			rotDet->rotateY(90.*degree);
			rotDet->rotateZ(phiDet);

			G4ThreeVector detPos(x_0 + ring*x_offset,0.,distance);
			detPos *= *rotDet;
			neutronDetectors[i]->MakeImprint(worldLV,detPos,rotDet);

			G4cout << "detectorPos[" << i << "] = TVector3("
			       << detPos.x() << ","
			       << detPos.y() << ","
			       << detPos.z() << ");" << G4endl;
		}
	}

	// Create a region
	G4LogicalVolume* logicVol = G4LogicalVolumeStore::GetInstance()->GetVolume("1inch-EJ309");
	G4Region* DetectorRegion = new G4Region("DetectorRegion");
	DetectorRegion->AddRootLogicalVolume(logicVol);
	*/
	//G4double CuDimA = 4835.*mm; // total length
  	//G4double CuDimB = 1050.*mm; // outer diameter

	/*G4Tubs *testShieldS = new G4Tubs("testShieldS",0.5*CuDimB + 5*mm,0.5*CuDimB + 65*mm,0.5*CuDimA,0.,360.*deg);
	G4LogicalVolume *testShieldLV = new G4LogicalVolume(testShieldS,nistManager->FindOrBuildMaterial("G4_Pb"),"CladdingRodLV",0,0,0);
	new G4PVPlacement(
                    0,               // no rotation
                    G4ThreeVector(0,0,0), // at (0,0,0)
                    testShieldLV,         // its logical volume
                    "testShieldPV",         // its name
                    worldLV,               // its mother  volume
                    false,           // no boolean operations
                    1,               // copy number
                    fCheckOverlaps); // checking overlaps
	*/

	//BuildDetectors(worldLV);
	
	// Print materials
	#ifndef NOT_USING_MPI
	if(G4MPImanager::GetManager()->GetRank()==0) G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	#else
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	#endif

	return worldPV;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::BuildDetectors(G4LogicalVolume *worldLV)
{
	const G4double PbWallThickness = 8.*cm;
	const G4double BoronWallThickness = 2.*mm;
	const G4double WallHalfX = 110.*mm;
	G4Material *fPb = nistManager->FindOrBuildMaterial("G4_Pb");
	//G4Box *ShieldingWallS = new G4Box("ShieldingWallS",PbWallThickness/2.,WallHalfX,WallHalfX);

	//G4Box *ShieldingWallS = new G4Box("ShieldingWallS",WallHalfX,WallHalfX,PbWallThickness/2.);
	G4Trd *ShieldingWallS = new G4Trd("ShieldingWallS",WallHalfX*0.9,WallHalfX,WallHalfX*0.9,WallHalfX,PbWallThickness/2.);
  	G4LogicalVolume *ShieldingWallLV = new G4LogicalVolume(ShieldingWallS,fPb,"PbShieldingWallLV",0,0,0);

  	G4double sideWallThickn = 2.*cm;
  	G4Box *SideWallS1 = new G4Box("ShieldingWallS1",sideWallThickn,WallHalfX,PbWallThickness/2.);
  	G4LogicalVolume *SideWallLV1 = new G4LogicalVolume(SideWallS1,fPb,"SideWallLV1",0,0,0);

  	G4Box *SideWallS2 = new G4Box("ShieldingWallS2",WallHalfX-2.*sideWallThickn,sideWallThickn,PbWallThickness/2.);
  	G4LogicalVolume *SideWallLV2 = new G4LogicalVolume(SideWallS2,fPb,"SideWallLV2",0,0,0);

	G4Box *B4CWallS = new G4Box("B4CWallS",WallHalfX,WallHalfX,BoronWallThickness/2.);
	G4Material *fB4C = nistManager->FindOrBuildMaterial("G4_BORON_CARBIDE");
  	G4LogicalVolume *B4CWallLV = new G4LogicalVolume(B4CWallS,fB4C,"B4CWallLV",0,0,0);

  	const G4int nRings = 9;
  	const G4double ring_dZ = 45.*cm;

	const G4double distance = 630.*mm;
	const G4int nDetectorClusters = 12;
	const G4int nDetPerCluster = 4;
	const G4int nDetectors = nDetectorClusters*nDetPerCluster*nRings;
	G4AssemblyVolume* neutronDetectors[nDetectors];
	char DetName[32];

	for(G4int i=0;i<nDetectors;i++) {
		snprintf(DetName,32,"EJ276_%d",i+1);
		neutronDetectors[i] = EJ276_detector(i+1, DetName, 6.*cm);
	}

	G4double phiCluster[] = {23.89,45.,66.11,90.+23.89,90.+45.,90.+66.11,180.+23.89,180.+45.,180.+66.11,270.+23.89,270.+45.,270.+66.11};
	for(G4int ring=0;ring<nRings;++ring) {
		G4double detZ0 = (ring-4)*ring_dZ;
		for(G4int cluster=0;cluster<nDetectorClusters;cluster++) {
			G4double dx = 33.*mm; // 3 mm gap between detectors

			G4ThreeVector position[nDetPerCluster];
			position[0] = G4ThreeVector(dx,dx,0.);
			position[1] = G4ThreeVector(dx,-dx,0.);
			position[2] = G4ThreeVector(-dx,dx,0.);
			position[3] = G4ThreeVector(-dx,-dx,0.);

			G4ThreeVector posPbWall(0.,0.,-(14.04+0.5*PbWallThickness));
			G4ThreeVector posB4CWall(0.,0.,-(14.04+PbWallThickness+0.5*BoronWallThickness));

			G4ThreeVector posCluster = G4ThreeVector(detZ0,0.,distance);
			G4RotationMatrix *rotPos = new G4RotationMatrix();
			rotPos->rotateY(90.*degree);
			rotPos->rotateZ(phiCluster[cluster]*degree);

			G4RotationMatrix *rotDet = new G4RotationMatrix();
			rotDet->rotateY(90.*degree);
			rotDet->rotateZ(phiCluster[cluster]*degree);

			for(G4int i=0;i<nDetPerCluster;i++) {
				position[i] += posCluster;
				position[i] *= *rotPos;
				G4int detNbr = ring*nDetectorClusters*nDetPerCluster + cluster*nDetPerCluster + i;
				neutronDetectors[detNbr]->MakeImprint(worldLV,position[i],rotDet,detNbr + 1);
			}

			posPbWall += posCluster;


			G4ThreeVector sideWallPos[4];
			sideWallPos[0] = posPbWall + G4ThreeVector((WallHalfX-sideWallThickn),0,PbWallThickness);
			sideWallPos[1] = posPbWall + G4ThreeVector(-(WallHalfX-sideWallThickn),0,PbWallThickness);
			sideWallPos[2] = posPbWall + G4ThreeVector(0,(WallHalfX-sideWallThickn),PbWallThickness);
			sideWallPos[3] = posPbWall + G4ThreeVector(0,-(WallHalfX-sideWallThickn),PbWallThickness);

			posPbWall *= *rotPos;
			G4Transform3D transformPb(*rotPos,posPbWall);
			new G4PVPlacement(transformPb, ShieldingWallLV, "ShieldingWallPV", worldLV, false, cluster, fCheckOverlaps);

			for(G4int i=0;i<4;i++) sideWallPos[i] *= *rotPos;
			G4RotationMatrix *rotSideWall = new G4RotationMatrix();
			rotSideWall->rotateY(90.*degree);
			rotSideWall->rotateZ(phiCluster[cluster]*degree);

			G4Transform3D transformSideWall1(*rotSideWall,sideWallPos[0]);
			new G4PVPlacement(transformSideWall1, SideWallLV1, "SideWallPV1", worldLV, false, cluster, fCheckOverlaps);

			G4Transform3D transformSideWall2(*rotSideWall,sideWallPos[1]);
			new G4PVPlacement(transformSideWall2, SideWallLV1, "SideWallPV2", worldLV, false, cluster, fCheckOverlaps);

			G4Transform3D transformSideWall3(*rotSideWall,sideWallPos[2]);
			new G4PVPlacement(transformSideWall3, SideWallLV2, "SideWallPV3", worldLV, false, cluster, fCheckOverlaps);

			G4Transform3D transformSideWall4(*rotSideWall,sideWallPos[3]);
			new G4PVPlacement(transformSideWall4, SideWallLV2, "SideWallPV4", worldLV, false, cluster, fCheckOverlaps);

			posB4CWall += posCluster;
			posB4CWall *= *rotPos;
			G4Transform3D transformB4C(*rotPos,posB4CWall);
			//new G4PVPlacement(transformB4C, B4CWallLV, "B4CWallPV", worldLV, false, cluster, fCheckOverlaps);

			 
		}
	}

	// Create a region
	G4LogicalVolume* logicVol = G4LogicalVolumeStore::GetInstance()->GetVolume("PlasticDetectorHouse");
	G4Region* DetectorRegion = new G4Region("DetectorRegion");
	DetectorRegion->AddRootLogicalVolume(logicVol);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
//================== Sensitive detectors =================================

  G4String ScintillatorSDname = "ScintillatorSD";
  G4String ScintillatorHCname = "ScintillatorHC";
  ScintilatorSD* sd = new ScintilatorSD(ScintillatorSDname,ScintillatorHCname);
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  //SetSensitiveDetector("5x5-EJ-309-LV", sd, true);
  //SetSensitiveDetector("1x2-EJ-309-LV", sd, true);
  if(G4LogicalVolumeStore::GetInstance()->GetVolume("EJ276_LV")) {
  	SetSensitiveDetector("EJ276_LV", sd, true);
  }

  G4double parsLinear[] = {1.0,0.0,0.,0.};
  sd->SetLightFunc(11,0,parsLinear); //electrons
  sd->SetLightFunc(-11,0,parsLinear); //positrons
  sd->SetLightFunc(22,0,parsLinear); //photons
  
  G4double parsProtons[] = {0.74787,2.4077,0.29866,1.0}; //parameters according to Enqvist et al. (DOI: 10.1016/j.nima.2013.03.032)
  sd->SetLightFunc(2212,1,parsProtons);


//================== cross section biasing for (alpha,n)-reactions =========
  /*GB01BOptrMultiParticleChangeCrossSection *aXSbias = new GB01BOptrMultiParticleChangeCrossSection("XSbias");

  G4double BiasFactor=1.0e+07;
  //G4double BiasFactor=1.0;
  G4String biasedParticleName="alpha";
  aXSbias->AddParticle(biasedParticleName,BiasFactor);

  //G4LogicalVolume* logicVol = G4LogicalVolumeStore::GetInstance()->GetVolume("FuelRodLV");
  //aXSbias->AttachTo(logicVol);
  //G4cout << " Attaching biasing operator to logical volume " << logicVol->GetName()<<" with biasFactor = "<<BiasFactor<<" for "<<biasedParticleName<<G4endl;

  G4LogicalVolume* logicVol = G4LogicalVolumeStore::GetInstance()->GetVolume("CladdingRodLV");
  aXSbias->AttachTo(logicVol); 
  G4cout << " Attaching biasing operator to logical volume " << logicVol->GetName()<<" with biasFactor = "<<BiasFactor<<" for "<<biasedParticleName<<G4endl;
*/
}

G4AssemblyVolume* DetectorConstruction::EJ309_5x5inch(G4int copyNbr, const char* name)
{
	const G4bool CheckOverlaps = true;
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();
	//-------------------Scintillator CuDimEnsions-------------------------------------------
	G4double ScintRad = 127./2.*mm;
	G4double ScintHeight = 127.*mm;
	G4double ScintHouseWall = 1.52*mm;
	G4double LightGuideRad = 5.5*cm;
	G4double LightGuideHeight = 31.97*mm;

	//Scintillator Housing
	G4int Nzplanes = 14;
	const G4double zPlane[14] = {0.,0.,120.9,120.9,159.,159.,165.35,165.35,175.35,175.35,236.35,278.4,386.4,386.4};
	const G4double rInner[14] = {0.,0.,0.,0.,0.,0.,0.};
	const G4double rOuter[14] = {0.,65.0,65.0,69.85,69.85,88.,88.,80.,80.,71.,71.,29.,29.,0.};
	G4Polycone *DetHouseS = new G4Polycone("DetHouseS",0.,360.*deg,Nzplanes,zPlane,rInner,rOuter);
	G4LogicalVolume *fLogicScintDetector = new G4LogicalVolume(DetHouseS, fAlu, "5inch-EJ309",0,0,0);

	//Inside Scintillator Housing

	G4int NzplanesIn = 8;
	const G4double zPlaneIn[8] = {ScintHouseWall+ScintHeight,ScintHouseWall+ScintHeight,ScintHouseWall+ScintHeight+LightGuideHeight,ScintHouseWall+ScintHeight+LightGuideHeight,236.35,278.4,384.88,384.88};
	const G4double rInnerIn[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
	const G4double rOuterIn[8] = {0.,ScintRad,ScintRad,71.-0.64,71.-0.64,29.-0.64,29.-0.64,0.};
	G4Polycone *PMT5inchInS = new G4Polycone("5inchPMTIns",0.,360.*deg,NzplanesIn,zPlaneIn,rInnerIn,rOuterIn);
	G4LogicalVolume* InsidePMT5inchLV = new G4LogicalVolume(PMT5inchInS, fVacuum, "InsidePMT",0,0,0);
	new G4PVPlacement(0,              // no rotation
	                  G4ThreeVector(0,0,0), // at (x,y,z) relative to the house
	                  InsidePMT5inchLV,   // its logical volume
	                  "PMT",       // its name
	                  fLogicScintDetector,        // its mother volume
	                  false,          // no boolean operations
	                  copyNbr,              // copy number
	                  CheckOverlaps); // checking overlaps
	G4VisAttributes *emptyVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
	InsidePMT5inchLV->SetVisAttributes(emptyVisAtt);

	char nameLV[32];
	snprintf(nameLV,32,"%s_LV",name);
  
	G4Tubs *Scintillator5inchS = new G4Tubs("Scintillator5inch",0.,ScintRad,ScintHeight/2.,0.*deg,360.*deg);
	G4LogicalVolume* ScintillatorLV = new G4LogicalVolume(Scintillator5inchS, fEJ309, "5x5-EJ-309-LV",0,0,0);
	/*G4VPhysicalVolume* ScintillatorPV =*/ new G4PVPlacement(0,              // no rotation
	                  G4ThreeVector(0,0,ScintHouseWall+ScintHeight/2.), // at (x,y,z) relative to the house
	                  ScintillatorLV,   // its logical volume
	                  name,       // its name
	                  fLogicScintDetector,        // its mother volume
	                  false,          // no boolean operations
	                  copyNbr,              // copy number
	                  CheckOverlaps); // checking overlaps

	G4VisAttributes *scintVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
	ScintillatorLV->SetVisAttributes(scintVisAtt);  

	G4Tubs *LightGuide5inchS = new G4Tubs("LightGuide5inchS",0.,LightGuideRad,LightGuideHeight/2.,0.*deg,360.*deg);
	G4LogicalVolume* fLogicLightGuide5inch = new G4LogicalVolume(LightGuide5inchS, fBoroSilicate, "LightGuide",0,0,0);
	new G4PVPlacement(0,              // no rotation
	                  G4ThreeVector(0,0,LightGuideHeight/2.+ScintHouseWall+ScintHeight), // at (x,y,z) relative to the house
	                  fLogicLightGuide5inch,   // its logical volume
	                  "LightGuide",       // its name
	                  InsidePMT5inchLV,        // its mother volume
	                  false,          // no boolean operations
	                  copyNbr,              // copy number
	                  CheckOverlaps); // checking overlaps
	G4VisAttributes *lightguideVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
	fLogicLightGuide5inch->SetVisAttributes(lightguideVisAtt);

	G4VisAttributes *VisAtt = new G4VisAttributes(G4Colour::Gray());
	//VisAtt->SetForceWireframe(true);
	fLogicScintDetector->SetVisAttributes(VisAtt);

	G4RotationMatrix rot0;
	G4ThreeVector pos0(0.,0.,0.);
	detectorAssembly->AddPlacedVolume(fLogicScintDetector,pos0, &rot0);

  //---- Detector shielding------------------------------------------
	G4double ThicknessPb = 8.*mm;
	G4double ThicknessCu = 2.*mm;
	G4double ThicknessTot = ThicknessPb + ThicknessCu;
	G4double CylinderLength = 120.9*mm;
	G4double CylinderInnerRad = 65.*mm;

	G4Material *fPb = nistManager->FindOrBuildMaterial("G4_Pb");
	//G4Material *fCu = nistManager->FindOrBuildMaterial("G4_Cu");
	G4Material *fCu = G4Material::GetMaterial("G4_AIR"); //test to remove the copper

	G4Tubs *shieldOuterCylS = new G4Tubs("shieldOuterCyl.solid",CylinderInnerRad,CylinderInnerRad+ThicknessTot,CylinderLength/2.,0.*deg,360.*deg);
	//G4Tubs *shieldInnerCylS = new G4Tubs("shieldInnerCyl.solid",CylinderInnerRad,CylinderInnerRad+ThicknessCu,CylinderLength/2.,0.*deg,360.*deg);
	//G4Tubs *shieldOuterLidS = new G4Tubs("shieldOuterLid.solid",0.,CylinderInnerRad+ThicknessTot,ThicknessTot/2.,0.*deg,360.*deg);
	//G4Tubs *shieldInnerLidS = new G4Tubs("shieldInnerLid.solid",0.,CylinderInnerRad+ThicknessCu,ThicknessCu/2.,0.*deg,360.*deg);

	G4LogicalVolume *shieldOuterCylLV = new G4LogicalVolume(shieldOuterCylS, fPb, "shieldOuterCyl.logic",0,0,0);
	//G4LogicalVolume *shieldInnerCylLV = new G4LogicalVolume(shieldInnerCylS, fCu, "shieldInnerCyl.logic",0,0,0);
	//G4LogicalVolume *shieldOuterLidLV = new G4LogicalVolume(shieldOuterLidS, fPb, "shieldOuterLid.logic",0,0,0);
	//G4LogicalVolume *shieldInnerLidLV = new G4LogicalVolume(shieldInnerLidS, fCu, "shieldInnerLid.logic",0,0,0);

	G4VisAttributes *shieldOuterVis = new G4VisAttributes(G4Colour::Gray());
	//G4VisAttributes *shieldInnerVis = new G4VisAttributes(G4Colour::Yellow());
	shieldOuterCylLV->SetVisAttributes(shieldOuterVis);
	//shieldInnerCylLV->SetVisAttributes(shieldInnerVis);
	//shieldOuterLidLV->SetVisAttributes(shieldOuterVis);
	//shieldInnerLidLV->SetVisAttributes(shieldInnerVis);

	// place the copper in the lead
	//new G4PVPlacement(0,              // no rotation
	//                  G4ThreeVector(0,0,0), // at (x,y,z) relative to the house
	//                  shieldInnerCylLV,   // its logical volume
	//                  "shieldOuterLid.phys",       // its name
	//                  shieldOuterCylLV,        // its mother volume
	//                  false,          // no boolean operations
	//                  copyNbr,              // copy number
	//                  CheckOverlaps); // checking overlaps

	//new G4PVPlacement(0,              // no rotation
	//                  G4ThreeVector(0,0,ThicknessTot/2.-ThicknessCu/2.), // at (x,y,z) relative to the house
	//                  shieldInnerLidLV,   // its logical volume
	//                  "shieldOuterLid.phys",       // its name
	//                  shieldOuterLidLV,        // its mother volume
	//                  false,          // no boolean operations
	//                  copyNbr,              // copy number
	//                  CheckOverlaps); // checking overlaps

	//G4ThreeVector posLid(0.,0.,-ThicknessTot/2.);
	//detectorAssembly->AddPlacedVolume(shieldOuterLidLV,posLid,&rot0);

	G4ThreeVector posCyl(0.,0.,+CylinderLength/2.);
	detectorAssembly->AddPlacedVolume(shieldOuterCylLV,posCyl,&rot0);

  	G4double InsertDimJ = 370.*mm;
	G4Box* PolyShield1 = new G4Box("PolyShield1",InsertDimJ/sqrt(8.),InsertDimJ/sqrt(8.),159./2.+ThicknessTot/2.);
	G4Tubs* PolyHole = new G4Tubs("PolyHole",0.,CylinderInnerRad+ThicknessTot,159./2.+ThicknessTot/2.+1.,0.*deg,360.*deg);
	G4SubtractionSolid *PolyShield = new G4SubtractionSolid("PolyShield",PolyShield1,PolyHole);

    G4LogicalVolume* PolyShieldLV
    = new G4LogicalVolume(
                          PolyShield,   //its solid
                          G4Material::GetMaterial("SWX201"),      //its material
                          "PolyShieldLV"); //its name

	G4RotationMatrix rotPoly;
	rotPoly.rotateZ(45.*deg);
	G4ThreeVector posPoly(0.,0.,159./2.-ThicknessTot/2.);
	detectorAssembly->AddPlacedVolume(PolyShieldLV,posPoly, &rotPoly);
    

  return detectorAssembly;
}

G4AssemblyVolume* DetectorConstruction::EJ309_3x3inch(G4int copyNbr, const char* name)
{
	const G4bool CheckOverlaps = false;
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();
	//-------------------Scintillator CuDimEnsions-------------------------------------------
	G4double ScintRad = 76./2.*mm;
	G4double ScintHeight = 76.*mm;
	G4double ScintHouseWall = 1.52*mm;
	G4double LightGuideRad = 70./2.*mm;
	G4double LightGuideHeight = 20*mm;

	//Scintillator Housing
	const G4int Nzplanes = 13;
	const G4double zPlane[Nzplanes] = {0.,0.,   70.9, 70.9,92.2,92.2,107.1,107.1,157.1,157.1,187.1,324.,324.};
	const G4double rInner[Nzplanes] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	const G4double rOuter[Nzplanes] = {0.,39.52,39.52,44.4,44.4,50.8,50.8, 41.6, 41.6, 41.6, 29.4, 29.4, 0.};
	G4Polycone *DetHouseS = new G4Polycone("DetHouseS",0.,360.*deg,Nzplanes,zPlane,rInner,rOuter);
	G4LogicalVolume *fLogicScintDetector = new G4LogicalVolume(DetHouseS, fAlu, "3inch-EJ309",0,0,0);

	//Inside Scintillator Housing

	G4int NzplanesIn = 8;
	const G4double zPlaneIn[8] = {ScintHouseWall+ScintHeight,ScintHouseWall+ScintHeight,ScintHouseWall+ScintHeight+LightGuideHeight,ScintHouseWall+ScintHeight+LightGuideHeight,157.1,187.1,322.48,322.48};
	const G4double rInnerIn[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
	const G4double rOuterIn[8] = {0.,ScintRad,ScintRad,41.6-0.64,41.6-0.64,29.4-0.64,29.4-0.64,0.};
	G4Polycone *PMT5inchInS = new G4Polycone("3inchPMTIns",0.,360.*deg,NzplanesIn,zPlaneIn,rInnerIn,rOuterIn);
	G4LogicalVolume* InsidePMT5inchLV = new G4LogicalVolume(PMT5inchInS, fVacuum, "InsidePMT",0,0,0);
	new G4PVPlacement(0,              // no rotation
	                  G4ThreeVector(0,0,0), // at (x,y,z) relative to the house
	                  InsidePMT5inchLV,   // its logical volume
	                  "PMT",       // its name
	                  fLogicScintDetector,        // its mother volume
	                  false,          // no boolean operations
	                  copyNbr,              // copy number
	                  CheckOverlaps); // checking overlaps
	G4VisAttributes *emptyVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
	InsidePMT5inchLV->SetVisAttributes(emptyVisAtt);

	char nameLV[32];
	snprintf(nameLV,32,"%s_LV",name);
  
	G4Tubs *Scintillator5inchS = new G4Tubs("Scintillator3inch",0.,ScintRad,ScintHeight/2.,0.*deg,360.*deg);
	G4LogicalVolume* ScintillatorLV = new G4LogicalVolume(Scintillator5inchS, fEJ309, "3x3-EJ-309-LV",0,0,0);
	/*G4VPhysicalVolume* ScintillatorPV =*/ new G4PVPlacement(0,              // no rotation
	                  G4ThreeVector(0,0,ScintHouseWall+ScintHeight/2.), // at (x,y,z) relative to the house
	                  ScintillatorLV,   // its logical volume
	                  name,       // its name
	                  fLogicScintDetector,        // its mother volume
	                  false,          // no boolean operations
	                  copyNbr,              // copy number
	                  CheckOverlaps); // checking overlaps

	G4VisAttributes *scintVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
	ScintillatorLV->SetVisAttributes(scintVisAtt);  

	G4Tubs *LightGuide5inchS = new G4Tubs("LightGuide5inchS",0.,LightGuideRad,LightGuideHeight/2.,0.*deg,360.*deg);
	G4LogicalVolume* fLogicLightGuide5inch = new G4LogicalVolume(LightGuide5inchS, fBoroSilicate, "LightGuide",0,0,0);
	new G4PVPlacement(0,              // no rotation
	                  G4ThreeVector(0,0,LightGuideHeight/2.+ScintHouseWall+ScintHeight), // at (x,y,z) relative to the house
	                  fLogicLightGuide5inch,   // its logical volume
	                  "LightGuide",       // its name
	                  InsidePMT5inchLV,        // its mother volume
	                  false,          // no boolean operations
	                  copyNbr,              // copy number
	                  CheckOverlaps); // checking overlaps
	G4VisAttributes *lightguideVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
	fLogicLightGuide5inch->SetVisAttributes(lightguideVisAtt);

	G4VisAttributes *VisAtt = new G4VisAttributes(G4Colour::Gray());
	//VisAtt->SetForceWireframe(true);
	fLogicScintDetector->SetVisAttributes(VisAtt);

	G4RotationMatrix rot0;
	G4ThreeVector pos0(0.,0.,0.);
	detectorAssembly->AddPlacedVolume(fLogicScintDetector,pos0, &rot0);
/*
  //---- Detector shielding------------------------------------------
	G4double ThicknessPb = 8.*mm;
	G4double ThicknessCu = 2.*mm;
	G4double ThicknessTot = ThicknessPb + ThicknessCu;
	G4double CylinderLength = 120.9*mm;
	G4double CylinderInnerRad = 65.*mm;

	G4Material *fPb = nistManager->FindOrBuildMaterial("G4_Pb");
	//G4Material *fCu = nistManager->FindOrBuildMaterial("G4_Cu");
	G4Material *fCu = G4Material::GetMaterial("G4_AIR"); //test to remove the copper

	G4Tubs *shieldOuterCylS = new G4Tubs("shieldOuterCyl.solid",CylinderInnerRad,CylinderInnerRad+ThicknessTot,CylinderLength/2.,0.*deg,360.*deg);
	//G4Tubs *shieldInnerCylS = new G4Tubs("shieldInnerCyl.solid",CylinderInnerRad,CylinderInnerRad+ThicknessCu,CylinderLength/2.,0.*deg,360.*deg);
	//G4Tubs *shieldOuterLidS = new G4Tubs("shieldOuterLid.solid",0.,CylinderInnerRad+ThicknessTot,ThicknessTot/2.,0.*deg,360.*deg);
	//G4Tubs *shieldInnerLidS = new G4Tubs("shieldInnerLid.solid",0.,CylinderInnerRad+ThicknessCu,ThicknessCu/2.,0.*deg,360.*deg);

	G4LogicalVolume *shieldOuterCylLV = new G4LogicalVolume(shieldOuterCylS, fPb, "shieldOuterCyl.logic",0,0,0);
	//G4LogicalVolume *shieldInnerCylLV = new G4LogicalVolume(shieldInnerCylS, fCu, "shieldInnerCyl.logic",0,0,0);
	//G4LogicalVolume *shieldOuterLidLV = new G4LogicalVolume(shieldOuterLidS, fPb, "shieldOuterLid.logic",0,0,0);
	//G4LogicalVolume *shieldInnerLidLV = new G4LogicalVolume(shieldInnerLidS, fCu, "shieldInnerLid.logic",0,0,0);

	G4VisAttributes *shieldOuterVis = new G4VisAttributes(G4Colour::Gray());
	//G4VisAttributes *shieldInnerVis = new G4VisAttributes(G4Colour::Yellow());
	shieldOuterCylLV->SetVisAttributes(shieldOuterVis);
	//shieldInnerCylLV->SetVisAttributes(shieldInnerVis);
	//shieldOuterLidLV->SetVisAttributes(shieldOuterVis);
	//shieldInnerLidLV->SetVisAttributes(shieldInnerVis);

	// place the copper in the lead
	//new G4PVPlacement(0,              // no rotation
	//                  G4ThreeVector(0,0,0), // at (x,y,z) relative to the house
	//                  shieldInnerCylLV,   // its logical volume
	//                  "shieldOuterLid.phys",       // its name
	//                  shieldOuterCylLV,        // its mother volume
	//                  false,          // no boolean operations
	//                  copyNbr,              // copy number
	//                  CheckOverlaps); // checking overlaps

	//new G4PVPlacement(0,              // no rotation
	//                  G4ThreeVector(0,0,ThicknessTot/2.-ThicknessCu/2.), // at (x,y,z) relative to the house
	//                  shieldInnerLidLV,   // its logical volume
	//                  "shieldOuterLid.phys",       // its name
	//                  shieldOuterLidLV,        // its mother volume
	//                  false,          // no boolean operations
	//                  copyNbr,              // copy number
	//                  CheckOverlaps); // checking overlaps

	//G4ThreeVector posLid(0.,0.,-ThicknessTot/2.);
	//detectorAssembly->AddPlacedVolume(shieldOuterLidLV,posLid,&rot0);

	G4ThreeVector posCyl(0.,0.,+CylinderLength/2.);
	detectorAssembly->AddPlacedVolume(shieldOuterCylLV,posCyl,&rot0);

  	G4double InsertDimJ = 370.*mm;
	G4Box* PolyShield1 = new G4Box("PolyShield1",InsertDimJ/sqrt(8.),InsertDimJ/sqrt(8.),159./2.+ThicknessTot/2.);
	G4Tubs* PolyHole = new G4Tubs("PolyHole",0.,CylinderInnerRad+ThicknessTot,159./2.+ThicknessTot/2.+1.,0.*deg,360.*deg);
	G4SubtractionSolid *PolyShield = new G4SubtractionSolid("PolyShield",PolyShield1,PolyHole);

    G4LogicalVolume* PolyShieldLV
    = new G4LogicalVolume(
                          PolyShield,   //its solid
                          G4Material::GetMaterial("SWX201"),      //its material
                          "PolyShieldLV"); //its name

	G4RotationMatrix rotPoly;
	rotPoly.rotateZ(45.*deg);
	G4ThreeVector posPoly(0.,0.,159./2.-ThicknessTot/2.);
	detectorAssembly->AddPlacedVolume(PolyShieldLV,posPoly, &rotPoly);
*/

  return detectorAssembly;
}

G4AssemblyVolume* DetectorConstruction::EJ309_1x2inch(G4int copyNbr, const char* name)
{
	//for now this will just be a cylinder inside an aluminium casing
	const G4bool CheckOverlaps = false;
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();
	//-------------------Scintillator CuDimEnsions-------------------------------------------
	G4double ScintRad = 25.4/2.*mm;
	G4double ScintHeight = 50.8*mm;
	G4double ScintHouseWall = 1.52*mm;

	//Scintillator Housing
	G4Tubs *SolidScintDetector = new G4Tubs("SolidScintDetector",0.,ScintRad+ScintHouseWall,ScintHeight/2.+ScintHouseWall/2.,0.*deg,360.*deg);
	G4LogicalVolume *LogicScintDetector = new G4LogicalVolume(SolidScintDetector, fAlu, "1inch-EJ309",0,0,0);

	//Inside Scintillator Housing

	char nameLV[32];
	snprintf(nameLV,32,"%s_LV",name);
  
	G4Tubs *Scintillator1inchS = new G4Tubs("Scintillator1inch",0.,ScintRad,ScintHeight/2.,0.*deg,360.*deg);
	G4LogicalVolume* ScintillatorLV = new G4LogicalVolume(Scintillator1inchS, fEJ309, "1x2-EJ-309-LV",0,0,0);
	/*G4VPhysicalVolume* ScintillatorPV =*/ new G4PVPlacement(0,              // no rotation
	                  G4ThreeVector(0,0,ScintHouseWall*0.5), // at (x,y,z) relative to the house
	                  ScintillatorLV,   // its logical volume
	                  name,       // its name
	                  LogicScintDetector,        // its mother volume
	                  false,          // no boolean operations
	                  copyNbr,              // copy number
	                  CheckOverlaps); // checking overlaps

	G4VisAttributes *scintVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
	ScintillatorLV->SetVisAttributes(scintVisAtt); 

	G4RotationMatrix rot0;
	G4ThreeVector pos0(0.,0.,0.);
	detectorAssembly->AddPlacedVolume(LogicScintDetector,pos0, &rot0);

  //---- Detector shielding------------------------------------------
	/*G4double ThicknessPb = 8.*mm;
	G4double CylinderLength = ScintHeight + ScintHouseWall+0.2*mm;
	G4double CylinderInnerRad = ScintRad+ScintHouseWall+0.1*mm;
	G4double CylinderOuterRad = CylinderInnerRad+ThicknessPb+0.1*mm;
	G4double LidOuterRad = CylinderOuterRad;

	G4Material *fPb = nistManager->FindOrBuildMaterial("G4_Pb");

	G4Tubs *shieldOuterCylS = new G4Tubs("shieldOuterCyl.solid",CylinderInnerRad,CylinderOuterRad,CylinderLength/2.,0.*deg,360.*deg);
	G4Tubs *shieldLidS = new G4Tubs("shieldLid.solid",0,LidOuterRad,ThicknessPb/2.,0.*deg,360.*deg);
	
	G4LogicalVolume *shieldOuterCylLV = new G4LogicalVolume(shieldOuterCylS, fPb, "shieldOuterCyl.logic",0,0,0);
	G4LogicalVolume *shieldLidLV = new G4LogicalVolume(shieldLidS, fPb, "shieldLid.logic",0,0,0);
	
	G4VisAttributes *shieldOuterVis = new G4VisAttributes(G4Colour::Gray());
	shieldOuterCylLV->SetVisAttributes(shieldOuterVis);
	shieldLidLV->SetVisAttributes(shieldOuterVis);

	G4ThreeVector posCyl(0.,0.,0.);
	detectorAssembly->AddPlacedVolume(shieldOuterCylLV,posCyl,&rot0);
	G4ThreeVector posLid(0.,0.,-0.5*(CylinderLength+ThicknessPb));
	detectorAssembly->AddPlacedVolume(shieldLidLV,posLid,&rot0);
    */

  return detectorAssembly;
}

G4AssemblyVolume* DetectorConstruction::EJ276_detector(G4int copyNbr, const char* name, G4double size)
{
	//for now this will just be a cylinder inside an aluminium casing
	const G4bool CheckOverlaps = true;
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();
	//-------------------Scintillator DimEnsions-------------------------------------------
	G4double ScintHeight = 25.*mm;
	G4double ScintHouseWall = 1.52*mm;

	//Scintillator Housing
	G4Box *SolidScintDetector = new G4Box("SolidScintDetector",0.5*size+ScintHouseWall,0.5*size+ScintHouseWall,0.5*(ScintHeight+ScintHouseWall));
	G4LogicalVolume *LogicScintDetector = new G4LogicalVolume(SolidScintDetector, fAlu, "PlasticDetectorHouse",0,0,0);

	//Inside Scintillator Housing
  
	G4Box *ScintillatorS = new G4Box("ScintillatorS",0.5*size,0.5*size,0.5*ScintHeight);
	G4LogicalVolume* ScintillatorLV = new G4LogicalVolume(ScintillatorS, fEJ276, "EJ276_LV",0,0,0);
	/*G4VPhysicalVolume* ScintillatorPV =*/ new G4PVPlacement(0,              // no rotation
	                  G4ThreeVector(0,0,ScintHouseWall*0.5), // at (x,y,z) relative to the house
	                  ScintillatorLV,   // its logical volume
	                  name,       // its name
	                  LogicScintDetector,        // its mother volume
	                  false,          // no boolean operations
	                  copyNbr,              // copy number
	                  CheckOverlaps); // checking overlaps

	G4VisAttributes *scintVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
	ScintillatorLV->SetVisAttributes(scintVisAtt); 

	G4RotationMatrix rot0;
	G4ThreeVector pos0(0.,0.,0.);
	detectorAssembly->AddPlacedVolume(LogicScintDetector,pos0, &rot0);

	return detectorAssembly;
}