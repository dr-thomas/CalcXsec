#include "./xsecBinsHelper.hxx"
#include "./xsecBinsHelper.cxx"

Float_t* GetMCEventRateFromFitIPS(TString inFMCStr){

	TFile* inFMC = new TFile(inFMCStr,"OPEN");
	TTree* selTreeMC = (TTree*)inFMC->Get("selectedEvents");

	Int_t cutBranch=-1;
	Int_t mectopology=-1;
	Float_t D1True=-1;//cos
	Float_t D1Rec=-1;
	Float_t D2True=-1;//p
	Float_t D2Rec=-1;
	Float_t weight=-1;
	Int_t IsOnWater=-1;

	selTreeMC->SetBranchAddress("cutBranch",&cutBranch);
	selTreeMC->SetBranchAddress("mectopology",&mectopology);
	selTreeMC->SetBranchAddress("D1True",&D1True);
	selTreeMC->SetBranchAddress("D1Rec",&D1Rec);
	selTreeMC->SetBranchAddress("D2True",&D2True);
	selTreeMC->SetBranchAddress("D2Rec",&D2Rec);
	selTreeMC->SetBranchAddress("weight",&weight);
	selTreeMC->SetBranchAddress("IsOnWater",&IsOnWater);

	Int_t nEntryMC = selTreeMC->GetEntriesFast();

	Float_t* out = new Float_t[19];
	for(int ii=0; ii<19; ii++){
		out[ii]=0.;
	}

	xsecBinsHelper* binHelper = new xsecBinsHelper();

	for(int iEvent=0; iEvent<nEntryMC; iEvent++){
		selTreeMC->GetEntry(iEvent);
		int binIndex = binHelper->GetBinIPS(D2True,D1True);
		if(binIndex<0) continue;
		//signal cut
		if(cutBranch==1 && mectopology<3 && IsOnWater==1) out[binIndex]+=weight;
	}
	inFMC->Close();
	return out;
}
