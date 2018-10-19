#include <TString.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TAxis.h>

#include "./util/p0dCCEventClass.C"
#include "./util/GetMCEventRateFromFitIPS.C"
#include "./covmat/covmat.hxx"
#include "./covmat/covmat.cxx"
//to clone: git clone tcampbell@ens-hpc.engr.colostate.edu:/home/other/tcampbell/git/CalcXsec.git

int GetFluxBinIndex(Float_t);

void CalcXsec(){
	//define files
	TString inFMCGenieStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/GenWithFlagGENIE.root";
	TString inFMCEffStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/GenWithFlag.root";
	TString inFFitStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/DataFits/NEUT/fitBaseOutNoReg.root";

	TFile* inFMCEff = new TFile(inFMCEffStr,"OPEN");
	TFile* inFMCGenie = new TFile(inFMCGenieStr,"OPEN");
	TFile* inFFit = new TFile(inFFitStr,"OPEN");

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
	
	//flux
	TFile *finFluxCov = new TFile("/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/flux_covariance_banff_13av1.1.root","OPEN");
	TAxis* numubBins = (TAxis*)finFluxCov->Get("nd5_anumode_numub_bins");
	TMatrixDSym* covInFlux   = (TMatrixDSym*)finFluxCov->Get("total_flux_cov");

	covMatD* fluxCov = new covMatD(numubBins->GetNbins());
	for(int ii=0;ii<numubBins->GetNbins();ii++){
		for(int jj=0;jj<numubBins->GetNbins();jj++){
			fluxCov->SetMat(ii,jj,(*covInFlux)(ii+30,jj+30));
      }
    }
    finFluxCov->Close();

	//fit
	TMatrixD* covInFit = (TMatrixD*)inFFit->Get("res_cov_matrix");
	covMatD* fitCov = new covMatD(numubBins->GetNbins());
	for(int ii=0; ii<covInFit->GetNcols(); ii++){
		for(int jj=0; jj<covInFit->GetNcols(); jj++){
			fitCov->SetMat(ii,jj,(*covInFit)(ii,jj));
		}
	}

	//intialize results (ngen, nsel, etc) 
	Float_t** nData = new Float_t*[400];
	Float_t** nSel = new Float_t*[400];
	Float_t** nGen = new Float_t*[400];
	for(int ii=0; ii<400; ii++){
		nData[ii] = new Float_t[19];
		nSel[ii] = new Float_t[19];
		nGen[ii] = new Float_t[19];
	}

	xsecBinsHelper* binHelper = new xsecBinsHelper();
	//loop over events, loop over toys
	//generated loop
	int nToys=10;//max 400
	for(int iEntry=0; iEntry<nEntriesTruth; iEntry++){
		if((iEntry%10000)==0) cout << iEntry*1.0/nEntriesTruth << endl;
		truthT->GetEntry(iEntry);
		int bin=binHelper->GetBinIPS(genEvt->truelepton_mom, genEvt->truelepton_costheta);
		if((!genEvt->IsOnWater) || bin<0) continue;

		fluxCov->SetSeed(134987);
		for(int iToy=0; iToy<nToys; iToy++){
			fluxCov->Throw();

			Float_t weight=1.0;
			for(int ii=0; ii<15; ii++){
				weight*=(*(genEvt->WeightsMatrix))(ii,iToy);
			}
			weight*=genEvt->weightHL;
			int fluxBin=GetFluxBinIndex(genEvt->nu_trueE);
			if(fluxBin<0) continue;
			weight*=(fluxCov->varVec[fluxBin]);
			if(weight>-1e-6 && weight<10.) nGen[iToy][bin]+=weight;
		}
		//selected loop
		if(iEntry<nEntriesD){
			defaultT->GetEntry(iEntry);
			int binSel=binHelper->GetBinIPS(selEvt->truelepton_mom, selEvt->truelepton_costheta);
			if((selEvt->IsOnWater)==1 && binSel>0 && (selEvt->topology)==0 && (selEvt->nu_pdg)==-14) {
				//TODO: FV cut here!?!?!?
				fluxCov->SetSeed(134987);
				for(int iToy=0; iToy<nToys; iToy++){
					fluxCov->Throw();
					Float_t weight=1.0;
					for(int ii=0; ii<15; ii++){
						weight*=(*(selEvt->WeightsMatrix))(ii,iToy);
					}
					weight*=selEvt->weightHL;
					int fluxBin=GetFluxBinIndex(selEvt->nu_trueE);
					if(fluxBin<0) continue;
					weight*=(fluxCov->varVec[fluxBin]);
					if(weight>-1e-6 && weight<10.) nSel[iToy][binSel]+=weight;
				}
			}
		} // selected events
	} // iEntry

	//get data event rates
	fitCov->SetSeed(12334987);
	for(int iToy=0; iToy<400; iToy++){
		fitCov->Throw();
	}


	//calc xsec
	Float_t* nomEvtRate = GetMCEventRateFromFitIPS();
	//draw results
}

int GetFluxBinIndex(Float_t nuE){
	Float_t FluxBins[12]={0.0,0.4,0.5,0.6,0.7,1.0,1.5,2.5,3.5,5.0,7.0,30.0};
	int out=0;
	if(!(nuE<=1 || nuE>=29999)){
		while(nuE>=(FluxBins[out]*1000)){
			out++;
		}
	}
	if (out<1||out>11) return -1;
	return (out-1);
}

