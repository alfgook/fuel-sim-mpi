{
	TChain *chain = new TChain("Sorted");

	const Double_t totalActivity = 4.30871e+15*4.; // activity for the background runs; 7 assemblies a 4,3E15 Bq
	//const Double_t initialEvents = 1.6e7; // number of initial events in Geant4
	const Double_t initialEvents = 19.E10; // number of initial events in Geant4
	Double_t scalingFactor = totalActivity/initialEvents;

	TProof *p = TProof::Open("workers=1");
	p->SetParameter("pScalingFactor",(Double_t) scalingFactor);
	chain->SetProof();

	//chain->Add("../../build/PWR_CopperCan_BG_19E10events_testPb_sorted.root"); // 4.e6 events

	chain->Add("../../build/testPb_final_step_sorted.root"); // 4.e6 events

	chain->Process("FuelAnalysis.C++");
}