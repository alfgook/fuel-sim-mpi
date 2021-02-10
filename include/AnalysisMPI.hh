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
/// @file AnalysisMPI.hh
/// @brief Define histograms

#ifndef ANALYSIS_MANAGER_H
#define ANALYSIS_MANAGER_H

#include "G4ThreeVector.hh"
#include <tools/histo/h1d>
#include <tools/histo/h2d>


#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

class AnalysisMPI {
public:
  ~AnalysisMPI();

  static AnalysisMPI* GetAnalysis();

  void Book();
  void EndOfRun();

  //void OpenFile(const G4String& fname);
  void OpenFile();
  void Write();
  void CloseFile(G4bool reset = true);

  void FillScintillatorHit(G4int eventID, G4int copyNbr, G4int PDGcode, G4double time, G4double light, G4double weight);
  void FillSplitEvent(G4int eventID, G4int PDGcode, G4double time, G4double KinE, G4double posX, G4double posY, G4double posZ, G4double dirX, G4double dirY, G4double dirZ, G4double weight);
  void SetInitial(G4double x, G4double y, G4double z) { initX = x; initY = y; initZ = z;}

  void ClearIncidentFlag();

private:
  AnalysisMPI();
  DISALLOW_COPY_AND_ASSIGN(AnalysisMPI);

  static G4ThreadLocal G4int fincidentFlag;

  G4double initX, initY, initZ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void AnalysisMPI::ClearIncidentFlag()
{
  fincidentFlag = false;
}

#endif
