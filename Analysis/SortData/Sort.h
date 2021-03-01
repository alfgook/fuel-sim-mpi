//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar 21 11:08:44 2020 by ROOT version 5.34/38
// from TTree FUEL-PIN/ResponseNtuple
// found on file: PWR_run1_actinides_t0.root
//////////////////////////////////////////////////////////

#ifndef Sort_h
#define Sort_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Sort {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EventID;
   Int_t           DetectorNbr;
   Int_t           PDGcode;
   Float_t         TimeOfEvent;
   Float_t         TimeInEvent;
   Float_t         Light;
   Float_t         Weight;
   Float_t         InitX;
   Float_t         InitY;
   Float_t         InitZ;

   // List of branches
   TBranch        *b_EventID;   //!
   TBranch        *b_DetectorNbr;   //!
   TBranch        *b_PDGcode;   //!
   TBranch        *b_TimeOfEvent;   //!
   TBranch        *b_TimeInEvent;   //!
   TBranch        *b_Light;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_InitX;   //!
   TBranch        *b_InitY;   //!
   TBranch        *b_InitZ;   //!

   Sort(TTree *tree=0);
   virtual ~Sort();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const char*);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Sort_cxx
Sort::Sort(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PWR_run1_actinides_t0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("PWR_run1_actinides_t0.root");
      }
      f->GetObject("FUEL-PIN",tree);

   }
   Init(tree);
}

Sort::~Sort()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Sort::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Sort::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Sort::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("DetectorNbr", &DetectorNbr, &b_DetectorNbr);
   fChain->SetBranchAddress("PDGcode", &PDGcode, &b_PDGcode);
   fChain->SetBranchAddress("TimeOfEvent", &TimeOfEvent, &b_TimeOfEvent);
   fChain->SetBranchAddress("TimeInEvent", &TimeInEvent, &b_TimeInEvent);
   fChain->SetBranchAddress("Light", &Light, &b_Light);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("InitX", &InitX, &b_InitX);
   fChain->SetBranchAddress("InitY", &InitY, &b_InitY);
   fChain->SetBranchAddress("InitZ", &InitZ, &b_InitZ);
   Notify();
}

Bool_t Sort::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Sort::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Sort::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Sort_cxx
