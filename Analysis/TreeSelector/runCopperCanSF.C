{
	TChain *chain = new TChain("Sorted");

	const Double_t totalActivity = 4.*3.89763e+07; // activity for the background runs; 7 assemblies a 4,3E15 Bq
	//const Double_t initialEvents = 1.6e7; // number of initial events in Geant4
	const Double_t initialEvents = 55.*1.e7; // number of initial events in Geant4
	Double_t scalingFactor = totalActivity/initialEvents;

	TProof *p = TProof::Open("workers=10");
	p->SetParameter("pScalingFactor",(Double_t) scalingFactor);
	chain->SetProof();

	chain->Add("/Data1/alf/test/testSF/step2/PWR_CopperCan_SF_asterix_sorted.root"); // 4.e6 events

	chain->Process("FuelAnalysis.C++");
}