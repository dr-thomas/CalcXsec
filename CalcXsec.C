#include <TString.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TTree.h>
#include "./util/p0dCCEventClass.C"
//to clone: git clone tcampbell@ens-hpc.engr.colostate.edu:/home/other/tcampbell/git/CalcXsec.git

void CalcXsec(){
	//define files
	TString inFMCStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/WaterSubFullNom.root";
	TString inFMCGenieStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/GenWithFlagGENIE.root";
	TString inFMCEffStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/GenWithFlag.root";
	TString inFFitStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/DataFits/NEUT/fitBaseOutNoReg.root";
	TString inFDataStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/WaterSubData.root";

	TFile* inFMC = new TFile(inFMCStr,"OPEN");
	TFile* inFMCEff = new TFile(inFMCEffStr,"OPEN");
	TFile* inFMCGenie = new TFile(inFMCGenieStr,"OPEN");
	TFile* inFFit = new TFile(inFFitStr,"OPEN");
	TFile* inFData = new TFile(inFDataStr,"OPEN");

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
	//draw results
}
