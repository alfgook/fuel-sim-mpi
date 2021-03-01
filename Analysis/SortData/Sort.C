#define Sort_cxx
#include "Sort.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream> 

#define NBR_DETECTORS 8
#define MAX_NBR_HITS 512

using namespace std;

void Sort::Loop(const char *OutPutFileName)
{
//   In a ROOT session, you can do:
//      root> .L Sort.C
//      root> Sort t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;



   Long64_t nentries = fChain->GetEntries();
   cout << "nentries = " << nentries << endl;


//-------leafs-------------------------
//   Int_t           EventID;
//   Int_t           DetectorNbr;
//   Int_t           PDGcode;
//   Double_t        TimeOfEvent;
//   Double_t        TimeInEvent;
//   Double_t        Light;
//   Double_t        Weight;

   Double_t  fWeight; // is the same for all hits etc in an event
   UInt_t  NbrOfHits;
   UShort_t  fDetectorNbr[MAX_NBR_HITS];
   Int_t     fPDGcode[MAX_NBR_HITS];
   Double_t  fLight[MAX_NBR_HITS];
   Double_t  fTime[MAX_NBR_HITS];
   Double_t  fInitX;
   Double_t  fInitY;
   Double_t  fInitZ;

   TFile *fOut = new TFile(OutPutFileName,"recreate");
   TTree *tree = new TTree("Sorted","Sorted");

   tree->Branch("eventWeight",&fWeight,"eventWeight/D");
   tree->Branch("NbrOfHits",&NbrOfHits,"NbrOfHits/i");
   tree->Branch("DetectorNbr",fDetectorNbr,"DetectorNbr[NbrOfHits]/s");
   tree->Branch("PDGcodes",fPDGcode,"PDGcodes[NbrOfHits]/I");
   tree->Branch("Light",fLight,"Light[NbrOfHits]/D");
   tree->Branch("Time",fTime,"Time[NbrOfHits]/D");
   tree->Branch("InitX",&fInitX,"InitX/D");
   tree->Branch("InitY",&fInitY,"InitY/D");
   tree->Branch("InitZ",&fInitZ,"InitZ/D");

   Int_t FormerEventID = -1;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      DetectorNbr--; // in geant the numbers go from 1 -- max, but here i need 0 -- (max-1)

      if(EventID!=FormerEventID) {
         // push back the former event to the output tree
         if(FormerEventID!=-1) tree->Fill();

         // reset variables for the new event
         fWeight = Weight;
         NbrOfHits = 0;
      }

      if(fWeight != Weight) cout << "weight changed within an event: should not happen!!!" << endl;
         
      if(NbrOfHits<MAX_NBR_HITS) {
         fDetectorNbr[NbrOfHits] = DetectorNbr;
         fPDGcode[NbrOfHits]     = PDGcode;
         fTime[NbrOfHits]        = TimeInEvent;
         fLight[NbrOfHits]       = Light;
         fInitX                  = InitX;
         fInitY                  = InitY;
         fInitZ                  = InitZ;
         ++NbrOfHits;
      } else {
         cout << "too many hits for event " << EventID << endl;
      };
      
      //cout << "NbrOfHits[" << DetectorNbr << "] = " << NbrOfHits[DetectorNbr] << endl;

     // if(jentry>100) break;
      if(!(jentry%10000)) cout << jentry << endl;

      FormerEventID = EventID;
   }

  fOut->Write();
  fOut->Close();
}
