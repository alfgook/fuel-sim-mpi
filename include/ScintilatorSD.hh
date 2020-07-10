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
// $Id: ScintilatorSD.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file ScintilatorSD.hh
/// \brief Definition of the ScintilatorSD class

#ifndef ScintilatorSD_h
#define ScintilatorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4Trajectory.hh"

#include "ScintilatorHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Scintilator sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

class ScintilatorSD : public G4VSensitiveDetector
{
  public:
    ScintilatorSD(const G4String& name, const G4String& HCname);
    virtual ~ScintilatorSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

	// other methods
	G4int 		   SetLightFunc(G4int,G4int,G4double*);

  private:
    ScintilatorHitsCollection* fHitsCollection;

	std::vector<G4int>		LightFuncPDG;
	std::vector<G4int>		LightFuncType;
	std::vector<G4double>	LightFuncPar0;
	std::vector<G4double>	LightFuncPar1;
	std::vector<G4double>	LightFuncPar2;
	std::vector<G4double>	LightFuncPar3;

	inline G4double Light(G4double T, G4int PDGcode);
	inline void InitLightFuncs();
};

inline G4double ScintilatorSD::Light(G4double T, G4int PDGcode) {

	if(T<0.) return 0.;

	G4double L =0.;
	G4int index = -1;
	for(size_t i=0;i<LightFuncPDG.size();++i) {
		//find the index of the current particle type in the list of light output functions
		if(PDGcode==LightFuncPDG.at(i)) {
			index = i;
			break;
		}
	}
	if(index == -1) return 0.; //no light out-put set for this particle type

	switch(LightFuncType.at(index))
	{
		case 0:
		L = LightFuncPar0.at(index)*T + LightFuncPar1.at(index);
		if ( verboseLevel>2 ) {
			G4cout <<G4endl<< "PDGcode " << PDGcode << " case 0 : T = " << T << " : L = " << L << G4endl;
			G4cout << "LightFuncPar0 = " << LightFuncPar0.at(index) << G4endl;
			G4cout << "LightFuncPar1 = " << LightFuncPar1.at(index) << G4endl;
		}
		break;

		case 1:
		L = LightFuncPar0.at(index)*T - LightFuncPar1.at(index)*(1.0-exp(-LightFuncPar2.at(index)*pow(T,LightFuncPar3.at(index))));
		if ( verboseLevel>2 ) {
			G4cout <<G4endl<< "PDGcode " << PDGcode << "case 1 : T = " << T << " L = " << L << G4endl;
			G4cout << "LightFuncPar0 = " << LightFuncPar0.at(index) << G4endl;
			G4cout << "LightFuncPar1 = " << LightFuncPar1.at(index) << G4endl;
			G4cout << "LightFuncPar2 = " << LightFuncPar2.at(index) << G4endl;
			G4cout << "LightFuncPar3 = " << LightFuncPar3.at(index) << G4endl;
		}
		break;

		case 2:
		L = LightFuncPar0.at(index)*T*T/(LightFuncPar1.at(index)+T);
		if ( verboseLevel>2 ) {
			G4cout <<G4endl<< "PDGcode " << PDGcode << "case 2 : T = " << T << " L = " << L << G4endl;
			G4cout << "LightFuncPar0 = " << LightFuncPar0.at(index) << G4endl;
			G4cout << "LightFuncPar1 = " << LightFuncPar1.at(index) << G4endl;
		}
		break;

		case 3:
		L = LightFuncPar0.at(index)*pow(T,LightFuncPar1.at(index));
		if ( verboseLevel>2 ) {
			G4cout <<G4endl<< "PDGcode " << PDGcode << "case 3 : T = " << T << " L = " << L << G4endl;
			G4cout << "LightFuncPar0 = " << LightFuncPar0.at(index) << G4endl;
			G4cout << "LightFuncPar1 = " << LightFuncPar1.at(index) << G4endl;
		}
		break;

		default: //should never be reached
		L = 0;
		if ( verboseLevel>2 ) {
			G4cout <<G4endl<< "PDGcode " << PDGcode << "case default : T = " << T << " L = " << L << G4endl;
		}
		break;
	}

	if(L<0.) L = 0.;

	return L;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int ScintilatorSD::SetLightFunc(G4int PDGcode, G4int type, G4double *pars)
{
	//First find out if a light-out function was already set for the requested PDGcode
	G4int index = -1;
	for(size_t i=0;i<LightFuncPDG.size();++i) {
		if(PDGcode==LightFuncPDG.at(i)) {
			index = i;
			break;
		}
	}

	if(index == -1) { //function is not set for the requested particle type. Add new one to the list
		LightFuncType.push_back(type);
	} else { //function was already set for the requested particle type. modify the list
		LightFuncType.at(index) = type;
	}

	switch(type)
	{
		case 0:
		if(index == -1) { //function is not set for the requested particle type. Add new one to the list
			LightFuncPDG.push_back(PDGcode);
			LightFuncPar0.push_back(pars[0]);
			LightFuncPar1.push_back(pars[1]);
			LightFuncPar2.push_back(0.);
			LightFuncPar3.push_back(0.);
		} else { //function was already set for the requested particle type. modify the list
			LightFuncPar0.at(index) = pars[0];
			LightFuncPar1.at(index) = pars[1];
			LightFuncPar2.at(index) = 0.;
			LightFuncPar3.at(index) = 0.;
		}
		break;

		case 1:
		if(index == -1) { //function is not set for the requested particle type. Add new one to the list
			LightFuncPDG.push_back(PDGcode);
			LightFuncPar0.push_back(pars[0]);
			LightFuncPar1.push_back(pars[1]);
			LightFuncPar2.push_back(pars[2]);
			LightFuncPar3.push_back(pars[3]);
		} else { //function was already set for the requested particle type. modify the list
			LightFuncPar0.at(index) = pars[0];
			LightFuncPar1.at(index) = pars[1];
			LightFuncPar2.at(index) = pars[2];
			LightFuncPar3.at(index) = pars[3];
		}
		break;

		case 2:
		if(index == -1) { //function is not set for the requested particle type. Add new one to the list
			LightFuncPar0.push_back(pars[0]);
			LightFuncPar1.push_back(pars[1]);
			LightFuncPar2.push_back(0.);
			LightFuncPar3.push_back(0.);
		} else { //function was already set for the requested particle type. modify the list
			LightFuncPar0.at(index) = pars[0];
			LightFuncPar1.at(index) = pars[1];
			LightFuncPar2.push_back(0.);
			LightFuncPar3.push_back(0.);
		}
		break;

		case 3:
		if(index == -1) { //function is not set for the requested particle type. Add new one to the list
			LightFuncPar0.push_back(pars[0]);
			LightFuncPar1.push_back(pars[1]);
			LightFuncPar2.push_back(0.);
			LightFuncPar3.push_back(0.);
		} else { //function was already set for the requested particle type. modify the list
			LightFuncPar0.at(index) = pars[0];
			LightFuncPar1.at(index) = pars[1];
			LightFuncPar2.push_back(0.);
			LightFuncPar3.push_back(0.);
		}
		break;

		default:
		G4cout << "ScintilatorSD::SetLightFunc() : syntax error : type " << type << "is not supported" << G4endl;
		G4cout << "*****************************************************************" << G4endl;
		G4cout << "Supported Light-output functions types" << G4endl;
		G4cout << "0: L = pars[0]*T + pars[1]" << G4endl;
		G4cout << "1: L = pars[0]*T - pars[1]*(1.0-exp(-pars[2]*pow(T,pars[3])))" << G4endl;
		G4cout << "2: L = pars[0]*T*T/(pars[1]+T);" << G4endl;
		G4cout << "3: L = pars[0]*pow(T,pars[1]);" << G4endl;
		G4cout << "*****************************************************************" << G4endl;
		return 1; //return 1 upon failure!
		
	}

	return 0; //return 0 upon success!

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void ScintilatorSD::InitLightFuncs()
{ //set default values for the most common light out-put functions (EJ-301 liquid scintillator)

	G4double pars1[] = {0.83,2.82,0.25,0.93}; //NIMA704 (2013) 104-110
	SetLightFunc(2212,1,pars1); //proton

	G4double pars2[] = {1.0,0.,0.,0.};
	SetLightFunc(11,0,pars2); //electrons
	SetLightFunc(-11,0,pars2); //positrons
	SetLightFunc(22,0,pars2); //photons

	G4double pars3[] = {0.535*0.138,0.535*75.786,0.001745,1.007};
	SetLightFunc(1000060120,1,pars3); //carbon-12
	SetLightFunc(1000060130,1,pars3); //carbon-13
	
	G4double pars4[] = {0.535*0.391,0.535*4.159,0.085,1.1};
	SetLightFunc(1000020040,1,pars4); //alpha-particles

	G4double pars5[] = {0.535*0.138,0.535*75.786/3.,0.001745,1.007};
	SetLightFunc(1000040090,1,pars5); //Be-9

	G4double pars6[] = {0.83,5.64,0.131,0.93}; //deuterons: L(T,deuteron) = 2.*L(T/2.,proton)
	SetLightFunc(1000010020,1,pars6); //deuterons

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
