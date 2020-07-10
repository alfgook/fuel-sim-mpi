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
//  File:   MySFDecay.hh                                                      //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   27 March 2019                                                     //
//  Description: performs emission of neutrons and gammas from spontaneous    //
//               fission of the parent nucleus.  The LLNL model for neutron   //
//               and gamma spectra is used.  Final state nuclear fragments    //
//               are not currently generated.                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef MySFDecay_h
#define MySFDecay_h 1

#include "G4NuclearDecay.hh"
#include "fissionEvent.h"
#include "FissionGenerator.hh"


class MySFDecay : public G4NuclearDecay
{
  public:
    MySFDecay(const G4ParticleDefinition* theParentNucleus,
              const G4double& theBR, const G4double& Qvalue,
              const G4double& excitation, const G4Ions::G4FloatLevelBase& flb);

    virtual ~MySFDecay();

    virtual G4DecayProducts* DecayIt(G4double);

    virtual void DumpNuclearInfo();

    void SetARM(G4bool onoff) {applyARM = onoff;}
 
  private:
    const G4double transitionQ;
    G4int parentZ;
    G4int parentA;
    G4bool applyARM;

    std::vector<G4int> ZAIDvalid;
    std::vector<G4int> ZAIDvalidFREYA;

    FissionGenerator *theFissionGenerator;
};

#endif

