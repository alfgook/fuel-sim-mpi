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
// $Id: ScintilatorHit.cc 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file ScintilatorHit.cc
/// \brief Implementation of the ScintilatorHit class

#include "ScintilatorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<ScintilatorHit>* ScintilatorHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintilatorHit::ScintilatorHit() : G4VHit(), fVolumeCopyNo(0), fParentID(-1), fPDGcode(0), fLight(0.), fTime(0.), fWeight(1.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintilatorHit::~ScintilatorHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintilatorHit::ScintilatorHit(const ScintilatorHit& right) : G4VHit()
{
  fVolumeCopyNo = right.fVolumeCopyNo;
  fParentID   = right.fParentID;
  fPDGcode   = right.fPDGcode;
  fLight      = right.fLight;
  fTime      = right.fTime;
  fWeight      = right.fWeight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintilatorHit::ScintilatorHit(G4int VolumeCopyNo,G4int ParentID, G4int PDGcode, G4double weight = 1.) : G4VHit(), fVolumeCopyNo(VolumeCopyNo), fParentID(ParentID), fPDGcode(PDGcode), fLight(0.), fTime(0.), fWeight(weight)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const ScintilatorHit& ScintilatorHit::operator=(const ScintilatorHit& right)
{
  fVolumeCopyNo = right.fVolumeCopyNo;
  fParentID   = right.fParentID;
  fPDGcode   = right.fPDGcode;
  fLight      = right.fLight;
  fTime      = right.fTime;
  fWeight      = right.fWeight;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ScintilatorHit::operator==(const ScintilatorHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintilatorHit::Draw()
{
/*  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }*/
}

void ScintilatorHit::AddLight(G4double dL)
{
	fLight += dL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintilatorHit::Print()
{
  G4cout << "  PDGcode: " << fPDGcode << " Light: " << G4BestUnit(fLight,"Energy") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
