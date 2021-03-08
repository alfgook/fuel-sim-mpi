#define FuelAnalysis_cxx
// The class definition in FuelAnalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("FuelAnalysis.C")
// root> T->Process("FuelAnalysis.C","some options")
// root> T->Process("FuelAnalysis.C+")
//


#include "FuelAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TParameter.h>

void FuelAnalysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void FuelAnalysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   TParameter<Double_t> *pScalingFactor = dynamic_cast<TParameter<Double_t> *>(fInput->FindObject("pScalingFactor"));
   if(!pScalingFactor) {
      scalingFactor = 1.;
      cout << "didn't find the parameter" << endl;
   } else {
      scalingFactor = pScalingFactor->GetVal();
   }

   Char_t name[64];
   Char_t title[128];

   hLightAll = new TH1D("hLightAll","Light (any detector);Light (MeVee);counts/sec",2000,0.,20.);
   hLightAll->Sumw2();
   GetOutputList()->Add(hLightAll);

   hLightAll_neutrons = new TH1D("hLightAll_neutrons","Light (any detector,only neutrons);Light (MeVee);counts/sec",2000,0.,20.);
   hLightAll_neutrons->Sumw2();
   GetOutputList()->Add(hLightAll_neutrons);

   hLightAll_gammas = new TH1D("hLightAll_gammas","Light (any detector,only gammas);Light (MeVee);counts/sec",2000,0.,20.);
   hLightAll_gammas->Sumw2();
   GetOutputList()->Add(hLightAll_gammas);

   hSingleNeutrons = new TH1D("hSingleNeutrons","hSingleNeutrons",NBR_DETECTORS,-0.5,NBR_DETECTORS-0.5);
   hSingleNeutrons->Sumw2();
   GetOutputList()->Add(hSingleNeutrons);
   hSingleGammas = new TH1D("hSingleGammas","hSingleGammas",NBR_DETECTORS,-0.5,NBR_DETECTORS-0.5);
   hSingleGammas->Sumw2();
   GetOutputList()->Add(hSingleGammas);
   hGNcoincs = new TH2D("hGNcoincs","hGNcoincs",NBR_DETECTORS,-0.5,NBR_DETECTORS-0.5,NBR_DETECTORS,-0.5,NBR_DETECTORS-0.5);
   hGNcoincs->Sumw2();
   GetOutputList()->Add(hGNcoincs);
   hNNcoincs  = new TH2D("hNNcoincs","hNNcoincs",NBR_DETECTORS,-0.5,NBR_DETECTORS-0.5,NBR_DETECTORS,-0.5,NBR_DETECTORS-0.5);
   hNNcoincs->Sumw2();
   GetOutputList()->Add(hNNcoincs);

   hToFany= new TH1D("hToFany","time-of-flight any-coinc. (time(1)-time(2s));time (ns);counts/sec",500,0.,500.);
   hToFany->Sumw2();
   GetOutputList()->Add(hToFany);

   hToFanyGN = new TH1D("hToFanyGN","time-of-flight #gamma-n-coinc. (time(n)-time(#gamma));time (ns);counts/sec",1200,-100.,500.);
   hToFanyGN->Sumw2();
   GetOutputList()->Add(hToFanyGN);

   hToFanyNN = new TH1D("hToFanyNN","time-of-flight n-n-coinc. (time(n)-time(n));time (ns);counts/sec",1000,0.,500.);
   hToFanyNN->Sumw2();
   GetOutputList()->Add(hToFanyNN);


   for(Int_t i=0;i<nDetClusters+1;i++) {
      snprintf(name,64,"hInitPosSN%d",i);
      snprintf(title,128,"hInitPosSN");
      hInitPosSN[i] = new TH2D(name,title,800,-400.,400.,800.,-400.,400.);
      hInitPosSN[i]->Sumw2();
      GetOutputList()->Add(hInitPosSN[i]);
   }

   for(Int_t i=0;i<nDetClusters+1;i++) {
      snprintf(name,64,"hInitPosSG%d",i);
      snprintf(title,128,"hInitPosSG");
      hInitPosSG[i] = new TH2D(name,title,800,-400.,400.,800.,-400.,400.);
      hInitPosSG[i]->Sumw2();
      GetOutputList()->Add(hInitPosSG[i]);
   }

   for(Int_t i=0;i<nDetClusters+1;i++) {
      snprintf(name,64,"hInitPosGN%d",i);
      snprintf(title,128,"hInitPosGN");
      hInitPosGN[i] = new TH2D(name,title,800,-400.,400.,800.,-400.,400.);
      hInitPosGN[i]->Sumw2();
      GetOutputList()->Add(hInitPosGN[i]);
   }

   for(Int_t i=0;i<nDetClusters+1;i++) {
      snprintf(name,64,"hInitPosNN%d",i);
      snprintf(title,128,"hInitPosNN");
      hInitPosNN[i] = new TH2D(name,title,800,-400.,400.,800.,-400.,400.);
      hInitPosNN[i]->Sumw2();
      GetOutputList()->Add(hInitPosNN[i]);
   }

   hLight = new TH1D*[NBR_DETECTORS];
   for(Int_t i=0;i<NBR_DETECTORS;i++) {
      snprintf(name,64,"hLight%d",i);
      snprintf(title,128,"Light;Light (MeVee);counts/sec");
      hLight[i] = new TH1D(name,title,2000,0.,20.);
      hLight[i]->Sumw2();
      GetOutputList()->Add(hLight[i]);
   }

   hLight_neutrons = new TH1D*[NBR_DETECTORS];
   for(Int_t i=0;i<NBR_DETECTORS;i++) {
      snprintf(name,64,"hLight_neutrons%d",i);
      snprintf(title,128,"Light_neutrons;Light (MeVee);counts/sec");
      hLight_neutrons[i] = new TH1D(name,title,2000,0.,20.);
      hLight_neutrons[i]->Sumw2();
      GetOutputList()->Add(hLight_neutrons[i]);
   }

   hLight_gammas = new TH1D*[NBR_DETECTORS];
   for(Int_t i=0;i<NBR_DETECTORS;i++) {
      snprintf(name,64,"hLight_gammas%d",i);
      snprintf(title,128,"Light_gammas;Light (MeVee);counts/sec");
      hLight_gammas[i] = new TH1D(name,title,2000,0.,20.);
      hLight_gammas[i]->Sumw2();
      GetOutputList()->Add(hLight_gammas[i]);
   }
   
   hToF = new TH1D**[NBR_DETECTORS];
   for(Int_t i=0;i<NBR_DETECTORS;i++) {
      hToF[i] = new TH1D*[NBR_DETECTORS];
      for(Int_t j=0;j<NBR_DETECTORS;j++) {
         if(i==j) continue;
         snprintf(name,64,"hToF_%d_%d",i,j);
         snprintf(title,128,"time-of-flight #gamma-n-coinc. (time(det%d) - time(det%d));time (ns);counts/sec",i,j);
         hToF[i][j] = new TH1D(name,title,400,-100.,100.);
         hToF[i][j]->Sumw2();
         GetOutputList()->Add(hToF[i][j]);
      }
   }
}

Bool_t FuelAnalysis::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);
   Double_t eventWeight = *rv_eventWeight;
   eventWeight *= scalingFactor; // to get the histograms in counts/sec
   UShort_t NbrOfHits = *rv_NbrOfHits;

   //cout << "event weigt =" << eventWeight << endl;
   //cout << "scalingFactor =" << scalingFactor << endl;
   /*if(eventWeight==0.) {
      cout << "aaa" << endl;
   } else {
      cout << "bbb" << endl;
   }*/

   const Double_t t_min_GN = 4.;
   const Double_t t_max_GN = 120.;

   // is the main peak due to cross-talk? most probably
   const Double_t t_min_NN = 5.;
   const Double_t t_max_NN = 150.;

   const Double_t A2 = pow(0.102,2.);
   const Double_t B2 = pow(0.102,2.);
   const Double_t C2 = pow(0.036,2.);

   const Double_t Lmin = 0.05; // MeVee

   for(UShort_t hit1=0;hit1<NbrOfHits;++hit1) {
      UShort_t det1 = DetectorNbr[hit1];
      Int_t cluster = det1/nDetPerCluster + 1;

      Int_t gamma1 = 0; // 0 means it is a neutron
      if(PDGcodes[hit1]==22 || PDGcodes[hit1]==11 || PDGcodes[hit1]==-11) {
         gamma1 = 1; // 1 means it is a gamma
      }

      //weight due to light output threshold
      Double_t RL1 = sqrt(A2 + B2/Light[hit1] + C2/pow(Light[hit1],2.));
      Double_t sigmaL1 = RL1/2.3548*Light[hit1];
      //Double_t wL1 = 0.5*( 1. - TMath::Erf( (Lmin - Light[hit1])/(sqrt(2.)*sigmaL1) ) );
      Double_t wL1 = 1.;

      //eventWeight *= wL1;
      if(Light[hit1]<Lmin) continue;

      hLight[det1]->Fill(Light[hit1],eventWeight*wL1);
      hLightAll->Fill(Light[hit1],eventWeight*wL1);

      if(gamma1) {
         hLightAll_gammas->Fill(Light[hit1],eventWeight*wL1);
         hLight_gammas[DetectorNbr[hit1]]->Fill(Light[hit1],eventWeight*wL1);
         hSingleGammas->Fill(DetectorNbr[hit1],eventWeight*wL1);

         hInitPosSG[cluster]->Fill(*InitX,*InitY,eventWeight*wL1);
         hInitPosSG[0]->Fill(*InitX,*InitY,eventWeight*wL1);
      } else {
         hLightAll_neutrons->Fill(Light[hit1],eventWeight*wL1);
         hLight_neutrons[DetectorNbr[hit1]]->Fill(Light[hit1],eventWeight*wL1);
         hSingleNeutrons->Fill(DetectorNbr[hit1],eventWeight*wL1);

         hInitPosSN[cluster]->Fill(*InitX,*InitY,eventWeight*wL1);
         hInitPosSN[0]->Fill(*InitX,*InitY,eventWeight*wL1);
      }

      if(NbrOfHits>1) {
         for(UShort_t hit2=hit1+1;hit2<NbrOfHits;++hit2) {
            UShort_t det2 = DetectorNbr[hit2];

            UShort_t d_max = std::max(det1,det2);
            UShort_t d_min = std::min(det1,det2);
            
            if(det1>=NBR_DETECTORS || det2>=NBR_DETECTORS) continue;

            //weight due to light output threshold
            Double_t RL2 = sqrt(A2 + B2/Light[hit2] + C2/pow(Light[hit1],2.));
            Double_t sigmaL2 = RL2/2.3548*Light[hit2];
            //Double_t wL2 = 0.5*( 1. - TMath::Erf( (Lmin - Light[hit2])/(sqrt(2.)*sigmaL2) ) );
            Double_t wL2 = 1.;

            //eventWeight *= wL2;
            if(Light[hit2]<Lmin) continue;
            //PDGcode = 2212 = proton
            //PDGcode = 22 = gamma
            //PDGcode = 11 = electron
            //PDGcode = -11 = positron
            
            Int_t gamma2 = 0; // 0 means it is a neutron
            if(PDGcodes[hit2]==22 || PDGcodes[hit2]==11 || PDGcodes[hit2]==-11) {
               gamma2 = 1; // 1 means it is a gamma
            }

            // gamma-neutron coincidences
            if((gamma2 && !gamma1) || (!gamma2 && gamma1)) {
               Double_t tof = Time[hit1] - Time[hit2];
               Int_t detG = det2;
               Int_t detN = det1;
               if(!gamma2 && gamma1) {
                  tof = Time[hit2] - Time[hit1];
                  detG = det1;
                  detN = det2;
               }
               //hToF[det1][det2]->Fill(tof,eventWeight);
               if(det1!=det2) {
                  hToFanyGN->Fill(tof,eventWeight*wL1*wL2);
                  hToF[detG][detN]->Fill(tof,eventWeight*wL1*wL2);
                  if(tof<t_max_GN && tof>=t_min_GN) hGNcoincs->Fill(det2,det1,eventWeight*wL1*wL2);
               }

               Int_t clusterG = detG/nDetPerCluster;
               Int_t clusterN = detN/nDetPerCluster;
               Int_t clusterGN = 0; // cluster 0 is if the coincidence comes from detectorsof different clusters
               if(clusterG==clusterN) {
                  clusterGN = clusterG + 1;
               }

               if(tof<t_max_GN && tof>=t_min_GN) hInitPosGN[clusterGN]->Fill(*InitX,*InitY,eventWeight*wL1*wL2);

            }

            // neutron-neutron coincidences
            if(!gamma1 && !gamma2) {
               Double_t tof = fabs(Time[hit2] - Time[hit1]);
               if(det1!=det2) {
                  hToFanyNN->Fill(tof,eventWeight*wL1*wL2); //(d2-d1)>=4 if and only if the detectors dont belong to the same cluster_max
                  if(tof<t_max_NN && tof>=t_min_NN) hNNcoincs->Fill(det2,det1,eventWeight*wL1*wL2);

                  Int_t cluster1 = det1/nDetPerCluster;
                  Int_t cluster2 = det2/nDetPerCluster;
                  Int_t clusterNN = 0; // cluster 0 is if the coincidence comes from detectorsof different clusters
                  if(cluster1==cluster2) {
                     clusterNN = cluster1 + 1;
                  }

                  if(tof<t_max_NN && tof>=t_min_NN) hInitPosNN[clusterNN]->Fill(*InitX,*InitY,eventWeight*wL1*wL2);
               }
            }

            Double_t tof = fabs(Time[hit2] - Time[hit1]);
            if(det1!=det2) hToFany->Fill(tof,eventWeight*wL1*wL2);

         }
      }
      
   }

   return kTRUE;
}

void FuelAnalysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void FuelAnalysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TString option = GetOption();
   if(!option.Length()) option = "out.root";

   TFile fOut(option.Data(),"recreate");
   GetOutputList()->Write();
   fOut.Close();
}