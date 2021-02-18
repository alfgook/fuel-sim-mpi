//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 26 14:28:15 2020 by ROOT version 6.18/04
// from TTree Sorted/Sorted
// found on file: PWR_run1_alphaN_sorted.root
//////////////////////////////////////////////////////////

#ifndef FuelAnalysis_h
#define FuelAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <TH1D.h>
#include <TH2D.h>

#define NBR_DETECTORS 72


class FuelAnalysis : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Double_t> rv_eventWeight = {fReader, "eventWeight"};
   TTreeReaderValue<UInt_t> rv_NbrOfHits = {fReader, "NbrOfHits"};
   TTreeReaderArray<UShort_t> DetectorNbr = {fReader, "DetectorNbr"};
   TTreeReaderArray<Int_t> PDGcodes = {fReader, "PDGcodes"};
   TTreeReaderArray<Double_t> Light = {fReader, "Light"};
   TTreeReaderArray<Double_t> Time = {fReader, "Time"};


   FuelAnalysis(TTree * /*tree*/ =0) { }
   virtual ~FuelAnalysis() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(FuelAnalysis,0);

private :

   TH1D *hLightAll;
   TH1D *hLightAll_neutrons;
   TH1D *hLightAll_gammas;
   TH1D **hLight;
   TH1D ***hToF;
   TH1D *hToFany;
   TH1D *hToFanyGN;
   TH1D *hToFanyNN;

   Double_t scalingFactor;

};

#endif

#ifdef FuelAnalysis_cxx
void FuelAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t FuelAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef FuelAnalysis_cxx
