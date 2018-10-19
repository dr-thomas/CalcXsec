#include <TString.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TTree.h>
#include "./util/p0dCCEventClass.C"
#include "./util/GetMCEventRateFromFitIPS.C"
//to clone: git clone tcampbell@ens-hpc.engr.colostate.edu:/home/other/tcampbell/git/CalcXsec.git

void CalcXsec(){
	//define files
	TString inFMCGenieStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/GenWithFlagGENIE.root";
	TString inFMCEffStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/GenWithFlag.root";
	TString inFFitStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/DataFits/NEUT/fitBaseOutNoReg.root";

	TFile* inFMCEff = new TFile(inFMCEffStr,"OPEN");
	TFile* inFMCGenie = new TFile(inFMCGenieStr,"OPEN");
	TFile* inFFit = new TFile(inFFitStr,"OPEN");

	//TODO: fit should have been done against SF->RFG+RPA tuning... investigate!
	//acutlally not really that important as signal parameters are free,
	//but definately TODO TODO need to calc data event rate with used MC 
	//from the fit
	
	//define trees
	TTree* truthT = (TTree*)inFMCEff->Get("outTtruth");
	p0dCCEvent* genEvt = new p0dCCEvent();
	genEvt->SetBranchAddresses(truthT, true, false);
	int nEntriesTruth = truthT->GetEntries();

	TTree* defaultT = (TTree*)inFMCEff->Get("outTdefault");
	p0dCCEvent* selEvt = new p0dCCEvent();
	selEvt->SetBranchAddresses(defaultT, false, false);
	int nEntriesD = defaultT->GetEntries();

	//Genie
	//TODO: add weight to Genie stuff, re-produce those comparisions...
	TTree* truthTG = (TTree*)inFMCGenie->Get("outTtruth");
	p0dCCEvent* genEvtGENIE = new p0dCCEvent();
	genEvtGENIE->SetBranchAddresses(truthTG, true, true);
	int nEntriesTruthG = truthTG->GetEntries();

	//set up cov matrices 
	//intilize results (ngen, nsel, etc) 
	//loop over events, loop over toys
	//calc xsec
	Float_t* nomEvtRate = GetMCEventRateFromFitIPS();
	//draw results
}
