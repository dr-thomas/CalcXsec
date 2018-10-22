#include <TString.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TRandom3.h>

#include "./util/p0dCCEventClass.C"
#include "./util/GetMCEventRateFromFitIPS.C"
#include "./covmat/covmat.hxx"
#include "./covmat/covmat.cxx"
#include "./util/suffstat.hxx"
#include "./util/suffstat.cxx"
#include "./draw/DrawXsec.C"
//to clone: git clone tcampbell@ens-hpc.engr.colostate.edu:/home/other/tcampbell/git/CalcXsec.git

int GetFluxBinIndex(Float_t);
bool isInWTFV(Float_t* pos);

void CalcXsec(){
	//define files
	TString inFMCGenieStr = "/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/XsecDrawFiles/GenWithFlagGENIE.root";
	//TString inFMCEffStr = "/Users/thomascampbell/Desktop/GenCheck/GenWithFlag.root";
	TString inFMCEffStr = "./GenWithFlag.root";
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
	//(is there even flux weighting for GENIE MC?!?!)
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
	fluxCov->Decompose();
    finFluxCov->Close();

	//fit
	TMatrixD* covInFit = (TMatrixD*)inFFit->Get("res_cov_matrix");
	TVectorD* priorVec = (TVectorD*)inFFit->Get("res_vector");
	covMatD* fitCov = new covMatD(covInFit->GetNcols());
	for(int ii=0; ii<covInFit->GetNcols(); ii++){
		for(int jj=0; jj<covInFit->GetNcols(); jj++){
			fitCov->SetMat(ii,jj,(*covInFit)(ii,jj));
		}
	}
	fitCov->Decompose();

	//intialize results (ngen, nsel, etc) 
	Float_t** nData = new Float_t*[400];
	Float_t** nSel = new Float_t*[400];
	Float_t** nGen = new Float_t*[400];
	for(int ii=0; ii<400; ii++){
		nData[ii] = new Float_t[19];
		nSel[ii] = new Float_t[19];
		nGen[ii] = new Float_t[19];
	}
	for(int iToy=0; iToy<400; iToy++){
		for(int ii=0; ii<19; ii++){
			nData[iToy][ii]=0.;
			nSel[iToy][ii]=0.;
			nGen[iToy][ii]=0.;
		}
	}

	xsecBinsHelper* binHelper = new xsecBinsHelper();
	//loop over events, loop over toys
	//generated loop
	int nToys=400;//max 400

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
				weight*=((*(genEvt->WeightsMatrix))(ii,iToy));
			}
			weight*=genEvt->weightHL;
			weight*=genEvt->weightSF2RFG;
			int fluxBin=GetFluxBinIndex(genEvt->nu_trueE);
			if(fluxBin>-1) weight*=(1.+(fluxCov->varVec[fluxBin]));
			if(weight>-1e-6 && weight<10.) nGen[iToy][bin]+=weight;
			else nGen[iToy][bin]+=1.;
		}
	} // iEntry
	//selected loop
	for(int iEntry=0; iEntry<nEntriesD; iEntry++){
		defaultT->GetEntry(iEntry);
		int binSel=binHelper->GetBinIPS(selEvt->truelepton_mom, selEvt->truelepton_costheta);
		if((selEvt->IsOnWater)==1 && binSel>-1 && (selEvt->topology)==0 && (selEvt->nu_pdg)==-14 && isInWTFV(selEvt->vtx_truepos)) {
			fluxCov->SetSeed(134987);
			for(int iToy=0; iToy<nToys; iToy++){
				fluxCov->Throw();
				Float_t weight=1.0;
				for(int ii=0; ii<15; ii++){
					weight*=(*(selEvt->WeightsMatrix))(ii,iToy);
				}
				weight*=selEvt->weightHL;
				weight*=selEvt->weightSF2RFG;
				int fluxBin=GetFluxBinIndex(selEvt->nu_trueE);
				if(fluxBin>-1) weight*=(1.+(fluxCov->varVec[fluxBin]));
				if(weight>-1e-6 && weight<10.) nSel[iToy][binSel]+=weight;
				else nSel[iToy][binSel]+=1.;
			}
		}
	} // selected events

	//get data event rates
	Float_t* nomEvtRate = GetMCEventRateFromFitIPS();
	fitCov->SetSeed(12334987);
	for(int iToy=0; iToy<nToys; iToy++){
		fitCov->Throw();
		for(int ii=0; ii<19; ii++){
			nData[iToy][ii]=((*priorVec)(binHelper->ips2full[ii])*nomEvtRate[ii]*(1.+(fitCov->varVec[binHelper->ips2full[ii]])));
		}
	}

	//integrated flux
	//correleted with efficiency using prefit covariacne, no correlations to data
	TFile* inFflux = new TFile("/Users/thomascampbell/p0dCCAnalysis/FitResults/Macros/tuned13av1.1/run5c/nd5_tuned13av1.1_13anom_run5c_antinumode_fine.root");
	Double_t FluxBinsPass[12]={0.0,0.4,0.5,0.6,0.7,1.0,1.5,2.5,3.5,5.0,7.0,30.0};
	TH1F* tempFluxHist = (TH1F*)inFflux->Get("enu_nd5_tuned13a_numub");
	TH1F* FluxHist = (TH1F*)tempFluxHist->Rebin(11,"FluxHist",FluxBinsPass);
	Float_t* nomFlux = new Float_t[11];
	for(int i=0; i<11; i++) nomFlux[i]=(Float_t)FluxHist->GetBinContent(i+1);

	Float_t* intFlux = new Float_t[400];
	for(int ii=0; ii<400; ii++){
		intFlux[ii] = 0.;
	}

	fluxCov->SetSeed(134987);
	for(int iToy=0; iToy<nToys; iToy++){
		fluxCov->Throw();
		for(int ii=0; ii<11; ii++){
			intFlux[iToy]+=(nomFlux[ii]*(1.+(fluxCov->varVec[ii]))*2.08);
		}
	}

	//number of water molequles
	//1902 kg from TN-82, 0.8% uncertainty
	Float_t nTargetsNom=6.36e+28;
	//1930 kg from TN-72
	Float_t nTargetsNomMC=6.46e+28; 
	TRandom3 randN;
	randN.SetSeed(13857993);
	Float_t* nTargets = new Float_t[400];
	for(int ii=0; ii<nToys; ii++){
		nTargets[ii] = (1.+0.008*randN.Gaus())*nTargetsNom;
	}

	Float_t* binWidth = binHelper->GetBinWidths();

	//draw results
	DrawXsec(nData,nSel,nGen,binWidth,intFlux,nTargets,nTargetsNomMC,nToys);
	//for drawing, copy and paste gross code into new macro and pass around 
	//the nesseary calculated stuff 
	//-> or just take style stuff and write better/neater 
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

bool isInWTFV(Float_t* pos){
	return (pos[0]>-836. && pos[0]<764. && pos[1]>-871. && pos[1]<869. && pos[2]>-2968 && pos[2]<-1264);
}
