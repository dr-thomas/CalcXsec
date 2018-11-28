#include <iostream>
#include <iomanip>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector.h>

#include "TColor.h"
#include "/Volumes/ThomasDrive/p0dCCAnalysis/FitResults/Macros/Phil_root-style-files-master/palettes.C"
#include "TPaletteAxis.h"

#include "../util/suffstat.hxx"
#include "../util/suffstat.cxx"


void DrawXsec(Float_t** nData, Float_t** nSel, Float_t** nGen, Float_t** nGenSF, Float_t* binWidth, Float_t* intFlux, Float_t* nTargets, Float_t nTargetsNomMC, int nToys){
	/*
	//temp print stuff
	suffStat** tempDataStat = new suffStat*[19];
	suffStat** tempEffStat = new suffStat*[19];
	for(int ii=0; ii<19; ii++){
		tempDataStat[ii] = new suffStat();
		tempEffStat[ii] = new suffStat();
	}
	suffStat* tempFluxStat = new suffStat();
	suffStat* tempTargetsStat = new suffStat();
	for(int ii=0; ii<nToys; ii++){
		tempFluxStat->Fill(intFlux[ii]);
		tempTargetsStat->Fill(nTargets[ii]);
		for(int jj=0; jj<19; jj++){
			tempDataStat[jj]->Fill(nData[ii][jj]/2.08*0.287);
			tempEffStat[jj]->Fill(nSel[ii][jj]/nGen[ii][jj]);
		}
	}
	for(int ii=0; ii<19; ii++){
		cout << ii << "&";
		cout << std::setprecision(4) << tempDataStat[ii]->GetMean();
		cout << "&";
		cout << std::setprecision(4) << tempEffStat[ii]->GetMean();
		cout << "&";
		cout << std::setprecision(4) << binWidth[ii];
		cout << "\\\\" << endl;
	}
	cout << (tempFluxStat->GetMean()*0.287) << endl;
	cout << tempTargetsStat->GetMean() << endl;
	return;
	*/

	//calc xsec
	suffStat** xsecStat = new suffStat*[19];
	suffStat** xsecStatNEUT = new suffStat*[19];
	suffStat** xsecStatNEUTSF = new suffStat*[19];
	for(int ii=0; ii<19; ii++){
		xsecStat[ii] = new suffStat(1e-39);
		xsecStatNEUT[ii] = new suffStat(1e-39);
		xsecStatNEUTSF[ii] = new suffStat(1e-39);
	}
	suffStat** xsecMomStat = new suffStat*[7];
	suffStat** xsecMomStatNEUT = new suffStat*[7];
	suffStat** xsecMomStatNEUTSF = new suffStat*[7];
	for(int ii=0; ii<7; ii++){
		xsecMomStat[ii] = new suffStat(1e-39);
		xsecMomStatNEUT[ii] = new suffStat(1e-39);
		xsecMomStatNEUTSF[ii] = new suffStat(1e-39);
	}
	Float_t* singleBin = new Float_t[400];
	Float_t* singleBinN = new Float_t[400];
	Float_t* singleBinNSF = new Float_t[400];
	for(int ii=0; ii<nToys; ii++){
		singleBin[ii]=0.;
		singleBinN[ii]=0.;
		singleBinNSF[ii]=0.;
	}

	//TODO: single diff shit is sloppy af, should re-write
	//
	
	//xsecs for covariance calc
	Float_t** xsec4cov = new Float_t*[400];
	Float_t** xsec4covMom = new Float_t*[400];
	for(int ii=0; ii<400; ii++){
		xsec4cov[ii] = new Float_t[19];
		xsec4covMom[ii] = new Float_t[7];
	}
	for(int ii=0; ii<400; ii++){
		for(int jj=0; jj<19; jj++){
			xsec4cov[ii][jj] = 0.;
		}
		for(int jj=0; jj<7; jj++){
			xsec4covMom[ii][jj] = 0.;
		}
	}

	for(int iToy=0; iToy<nToys; iToy++){
		Float_t xsecMom=0.;
		Float_t xsecMomN=0.;
		Float_t xsecMomNSF=0.;
		int nBinsMom=2;
		int nDrawnP=0;
		int pIndex=0;
		for(int ii=0; ii<19; ii++){
			Float_t xsec = nData[iToy][ii]/(nSel[iToy][ii]/nGen[iToy][ii]);
			xsec*=1.0/binWidth[ii];
			xsec*=1.0/intFlux[iToy];
			xsec*=1.0/nTargets[iToy];
			xsecStat[ii]->Fill(xsec);
			xsec4cov[iToy][ii] = xsec*(1e40);
			singleBin[iToy]+=xsec;
			Float_t xsecN = nGen[iToy][ii]/(binWidth[ii]*intFlux[iToy]);
			Float_t xsecNSF = nGenSF[iToy][ii]/(binWidth[ii]*intFlux[iToy]);
			xsecN*=1.0/nTargetsNomMC;
			xsecNSF*=1.0/nTargetsNomMC;
			xsecStatNEUT[ii]->Fill(xsecN);
			xsecStatNEUTSF[ii]->Fill(xsecNSF);
			singleBinN[iToy]+=xsecN;
			singleBinNSF[iToy]+=xsecNSF;
			//Momentum single (Note: not normalized by bin width)
			xsec*=binWidth[ii];
			xsecN*=binWidth[ii];
			xsecNSF*=binWidth[ii];
			if(nDrawnP<nBinsMom){
				xsecMom+=xsec;
				xsecMomN+=xsecN;
				xsecMomNSF+=xsecNSF;
				nDrawnP++;
				if(ii==18) {
					xsecMomStat[pIndex]->Fill(xsecMom);
					xsec4covMom[iToy][pIndex] = xsecMom*(1e40);
					xsecMomStatNEUT[pIndex]->Fill(xsecMomN);
					xsecMomStatNEUTSF[pIndex]->Fill(xsecMomNSF);
				}
			} else {
				nDrawnP=0;
				xsecMomStat[pIndex]->Fill(xsecMom);
				xsec4covMom[iToy][pIndex] = xsecMom*(1e40);
				xsecMomStatNEUT[pIndex]->Fill(xsecMomN);
				xsecMomStatNEUTSF[pIndex]->Fill(xsecMomNSF);
				xsecMom=0.;
				xsecMomN=0.;
				xsecMomNSF=0.;
				pIndex++;
				if(ii==16) nBinsMom=2;
				else nBinsMom=3;
				//add again
				xsecMom+=xsec;
				xsecMomN+=xsecN;
				xsecMomNSF+=xsecNSF;
				nDrawnP++;
			}
		}
	}

	//calc final cov matrix 
	Float_t** xsecCov = new Float_t*[19];
	Float_t* xsecCovMeans = new Float_t[19];
	for(int ii=0; ii<19; ii++){
		xsecCov[ii] = new Float_t[19];
		xsecCovMeans[ii]=xsecStat[ii]->GetMean()*(1e40);
	}
	for(int ii=0; ii<19; ii++){
		for(int jj=0; jj<19; jj++){
			xsecCov[ii][jj] = 0.;
		}
	}
	for(int iToy=0; iToy<nToys; iToy++){
		for(int ii=0; ii<19; ii++){
			for(int jj=0; jj<19; jj++){
				xsecCov[ii][jj]+=(xsec4cov[iToy][ii]-xsecCovMeans[ii])*(xsec4cov[iToy][jj]-xsecCovMeans[jj])/(xsecCovMeans[ii]*xsecCovMeans[jj]);
			}
		}
	}

	for(int ii=0; ii<19; ii++){
		for(int jj=0; jj<19; jj++){
			xsecCov[ii][jj]*=(1.0/((Float_t)nToys));
		}
	}


	Float_t** xsecCovMom = new Float_t*[7];
	Float_t* xsecCovMeansMom = new Float_t[7];
	for(int ii=0; ii<7; ii++){
		xsecCovMom[ii] = new Float_t[7];
		xsecCovMeansMom[ii]=xsecMomStat[ii]->GetMean()*(1e40);
	}
	for(int ii=0; ii<7; ii++){
		for(int jj=0; jj<7; jj++){
			xsecCovMom[ii][jj]=0.;
		}
	}
	for(int iToy=0; iToy<nToys; iToy++){
		for(int ii=0; ii<7; ii++){
			for(int jj=0; jj<7; jj++){
				xsecCovMom[ii][jj]+=(xsec4covMom[iToy][ii]-xsecCovMeansMom[ii])*(xsec4covMom[iToy][jj]-xsecCovMeansMom[jj])/(xsecCovMeansMom[ii]*xsecCovMeansMom[jj]);
			}
		}
	}
	for(int ii=0; ii<7; ii++){
		for(int jj=0; jj<7; jj++){
			xsecCovMom[ii][jj]*=(1.0/((Float_t)nToys));
		}
	}

	cout << "xsec results in data (x+-error | NEUT): " << endl;
	for(int ii=0; ii<19; ii++){
		cout << xsecStat[ii]->GetMean() << " +- " << xsecStat[ii]->GetRMS();
		cout << " | " << xsecStatNEUT[ii]->GetMean() << endl;
	}


	cout << "Drawing Cross Section Plots" << endl;

	Float_t CosBinDraw[9][4]={
		{1,1,1,1},//only one bin here
		{0.84,0.94,1,1},//only three bins here
		{0.85,0.92,0.96,1},//Col. Forbins...
		{0.88,0.93,0.97,1},
		{0.90,0.94,0.97,1},
		{0.91,0.95,0.97,1},
		{0.94,0.97,0.98,1},
		{0.95,0.98,1,1},//only three bins here
		{1,1,1,1}//only one bin here
	};

	TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

	//NUISANCE
	TFile* inFnuisance = new TFile("rootFiles/nuiscompAnuCiroScaled.root", "OPEN");
	TH1F* nuisNeutRes = (TH1F*) inFnuisance->Get("T2K_CC0pi_XSec_2DPcos_anu_P0D_MC");


	TH1F** CosHists = new TH1F*[7];
	TH1F** CosHistsNEUT = new TH1F*[7];
	TH1F** CosHistsNEUTSF = new TH1F*[7];
	TH1F** CosHistsNuisNEUT = new TH1F*[7];

	for(int ii=0; ii<7; ii++){
		int tempNbins=0;
		if(ii==0 || ii==6){
			tempNbins=2;
		} else {
			tempNbins=3;
		}
		TString labelStr = "CosHists";
		labelStr += ii;
		CosHists[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);
		labelStr = "CosHistsNEUT";
		labelStr += ii;
		CosHistsNEUT[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);
		labelStr = "CosHistsNEUTSF";
		labelStr += ii;
		CosHistsNEUTSF[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);
		labelStr = "CosHistsNuisNEUT";
		labelStr += ii;
		CosHistsNuisNEUT[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);
	}

	int nFilled=0;
	for(int ii=0; ii<7; ii++){
		int tempNbins=0;
		if(ii==0 || ii==6){
			tempNbins=2;
		} else {
			tempNbins=3;
		}
		for(int jj=0; jj<tempNbins; jj++){
			CosHists[ii]->SetBinContent(jj+1,xsecStat[nFilled]->GetMean());
			CosHists[ii]->SetBinError(jj+1,xsecStat[nFilled]->GetRMS());
			CosHistsNEUT[ii]->SetBinContent(jj+1,xsecStatNEUT[nFilled]->GetMean());
			CosHistsNEUTSF[ii]->SetBinContent(jj+1,xsecStatNEUTSF[nFilled]->GetMean());
			CosHistsNuisNEUT[ii]->SetBinContent(jj+1, nuisNeutRes->GetBinContent(nFilled+1));
			nFilled++;
		}
	}

	// Print single bin results
	suffStat* sbStat = new suffStat(1e-39);
	suffStat* sbStatN = new suffStat(1e-39);
	suffStat* sbStatNSF = new suffStat(1e-39);
	for(int ii=0; ii<nToys; ii++){
		sbStat->Fill(singleBin[ii]);
		sbStatN->Fill(singleBinN[ii]);
		sbStatNSF->Fill(singleBinNSF[ii]);
	}
	Float_t singleBinNuisN = 0.;
	for(int ii=0; ii<19; ii++){
		singleBinNuisN += nuisNeutRes->GetBinContent(ii+1);
	}
	cout << "Single Bin Results:" << endl;
	cout << "Data: " << sbStat->GetMean() << " +- " << sbStat->GetRMS() << endl;
	cout << "NEUT: " << sbStatN->GetMean() << endl;
	cout << "NEUT SF: " << sbStatNSF->GetMean() << endl;
	cout << "NUISANCE NEUT: " << singleBinNuisN << endl;


	//Cos by momentum xsec slice 
	TString PBinsStr[9]={"0","400","530","670","800","1000","1380","2010","3410"};
	TString CosTitleStr;

	gStyle->SetOptStat(0);

	TCanvas* c;

	for(int ii=0; ii<7; ii++){
		c = new TCanvas;
		gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
		CosHists[ii]->GetYaxis()->SetTitleOffset(1.4);

		CosTitleStr="";
		CosTitleStr+=PBinsStr[ii+1];
		CosTitleStr+=" < P_{#mu} < ";
		CosTitleStr+=PBinsStr[ii+2];

		CosHists[ii]->GetXaxis()->SetTitle("sel. #mu^{+} cos#theta");
		CosHists[ii]->GetYaxis()->SetTitle("#frac{d#sigma}{dpd(cos#theta)} (cm^{2}/MeV/H_{2}O molecule)");
		CosHists[ii]->SetMarkerStyle(8);
		CosHists[ii]->SetMarkerSize(1);
		CosHists[ii]->SetMinimum(0);
		CosHists[ii]->SetTitle(CosTitleStr);

		CosHistsNEUT[ii]->SetLineColor(kRed+2);
		CosHistsNEUT[ii]->SetLineStyle(7);

		CosHistsNuisNEUT[ii]->SetLineColor(kBlack);
		CosHistsNuisNEUT[ii]->SetLineStyle(7);

		CosHistsNEUTSF[ii]->SetLineColor(kGreen+2);
		CosHistsNEUTSF[ii]->SetLineStyle(7);

		leg = new TLegend(0.1,0.1,0.9,0.9);
		leg->AddEntry(CosHists[ii],"Data","pel");
		leg->AddEntry(CosHistsNEUT[ii],"NEUT RFG+RPA","l");
		leg->AddEntry(CosHistsNuisNEUT[ii],"NUISANCE NEUT","l");
		leg->AddEntry(CosHistsNEUTSF[ii],"NEUT SF","l");

		CosHists[ii]->Draw("PE0");
		CosHistsNEUT[ii]->Draw("same");
		CosHistsNuisNEUT[ii]->Draw("same");
		CosHistsNEUTSF[ii]->Draw("same");

		TString printStr = "./plots/XsecSlice";
		printStr += ii;
		printStr += ".pdf";
		c->Print(printStr);
	}

	c = new TCanvas;
	leg->Draw();
	c->Print("./plots/XsecSliceLegend.pdf");
	
	// single diff
	Float_t PBinsDraw[8]={400,530,670,800,1000,1380,2010,3410};
	Float_t* PBinWidths = new Float_t[7];
	for(int ii=0; ii<7; ii++){
		PBinWidths[ii]=(PBinsDraw[ii+1]-PBinsDraw[ii]);
	}
	TH1F* MomFitResult = new TH1F("MomFitResult","",7,PBinsDraw);
	TH1F* MomFitResultNEUT = new TH1F("MomFitResultNEUT","",7,PBinsDraw);
	TH1F* MomFitResultNEUTSF = new TH1F("MomFitResultNEUTSF","",7,PBinsDraw);
	for(int ii=0; ii<7; ii++){
		MomFitResult->SetBinContent(ii+1, xsecMomStat[ii]->GetMean()/PBinWidths[ii]);
		MomFitResult->SetBinError(ii+1, xsecMomStat[ii]->GetRMS()/PBinWidths[ii]);
		MomFitResultNEUT->SetBinContent(ii+1, xsecMomStatNEUT[ii]->GetMean()/PBinWidths[ii]);
		MomFitResultNEUTSF->SetBinContent(ii+1, xsecMomStatNEUTSF[ii]->GetMean()/PBinWidths[ii]);
	}

	c = new TCanvas;
	gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
	MomFitResult->GetYaxis()->SetTitleOffset(1.4);

	MomFitResult->GetXaxis()->SetTitle("sel. #mu^{+} momentum (MeV/c)");
	MomFitResult->GetYaxis()->SetTitle("#frac{d#sigma}{dp} (cm^{2}/MeV/H_{2}O molecule)");
	MomFitResult->SetMarkerStyle(8);
	MomFitResult->SetMarkerSize(1);
	MomFitResult->SetMarkerColor(kBlue+3);
	MomFitResultNEUT->SetLineColor(kRed+2);
	MomFitResultNEUT->SetLineStyle(7);
	MomFitResultNEUTSF->SetLineColor(kGreen+2);
	MomFitResultNEUTSF->SetLineStyle(7);

	leg = new TLegend(0.7,0.7,0.9,0.9);
	leg->AddEntry(MomFitResult,"Data","pel");
	leg->AddEntry(MomFitResultNEUT,"NEUT RFG+RPA","l");
	leg->AddEntry(MomFitResultNEUTSF,"NEUT SF","l");

	MomFitResult->Draw("PE0");
	MomFitResultNEUT->Draw("same");
	MomFitResultNEUTSF->Draw("same");
	leg->Draw("same");
	c->Print("./plots/MomOverlay.pdf");

	// Covariance matrices
	
	TMatrixD FitCov(19,19);
	for(int ii=0; ii<19; ii++){
		for(int jj=0; jj<19; jj++){
			FitCov(ii,jj)=xsecCov[ii][jj];
		}
	}
	SetUserPalette(2);
	gStyle->SetOptStat(0);

	c = new TCanvas("c","",700,600);
	TH2D* tempFitCovHist = new TH2D(FitCov);
	tempFitCovHist->GetZaxis()->SetRangeUser(-0.3,0.3);
	tempFitCovHist->Draw("colz");

	gPad->Update();
	gPad->SetRightMargin(0.15);
	TPaletteAxis *palette = (TPaletteAxis*)tempFitCovHist->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.85);
	palette->SetX2NDC(0.9);
	palette->SetLabelSize(0.01);
	palette->SetLabelColor(kGreen);
	gPad->Modified();
	gPad->Update();

	c->Print("./plots/xsecCovMatrix.pdf");
	delete c;
	delete palette;
	delete tempFitCovHist;

	// Momentum single
	TMatrixD FitCovMom(7,7);
	for(int ii=0; ii<7; ii++){
		for(int jj=0; jj<7; jj++){
			FitCovMom(ii,jj)=xsecCovMom[ii][jj];
		}
	}
	SetUserPalette(2);
	gStyle->SetOptStat(0);

	c = new TCanvas("c","",700,600);
	tempFitCovHist = new TH2D(FitCovMom);
	tempFitCovHist->GetZaxis()->SetRangeUser(-0.3,0.3);
	tempFitCovHist->Draw("colz");

	gPad->Update();
	gPad->SetRightMargin(0.15);
	palette = (TPaletteAxis*)tempFitCovHist->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.85);
	palette->SetX2NDC(0.9);
	palette->SetLabelSize(0.01);
	palette->SetLabelColor(kGreen);
	gPad->Modified();
	gPad->Update();

	c->Print("./plots/xsecCovMatrixMom.pdf");
	delete c;

	//Chi2s
	FitCov.Invert();
	TVectorD chiVec1(19);
	TVectorD chiVec2(19);

	for(int ii=0; ii<19; ii++){
		chiVec1(ii) = (xsecStat[ii]->GetMean() - xsecStatNEUT[ii]->GetMean())/xsecStat[ii]->GetMean();
	}
	chiVec2 = FitCov*chiVec1;
	Float_t tempChi2=0;
	for(int ii=0; ii<19; ii++){
		tempChi2 += chiVec1(ii)*chiVec2(ii);
	}
	cout << "Chi2 data to NEUT: " << tempChi2/19. << endl;

	for(int ii=0; ii<19; ii++){
		chiVec1(ii) = (xsecStat[ii]->GetMean() - xsecStatNEUTSF[ii]->GetMean())/xsecStat[ii]->GetMean();
	}
	chiVec2 = FitCov*chiVec1;
	tempChi2=0;
	for(int ii=0; ii<19; ii++){
		tempChi2 += chiVec1(ii)*chiVec2(ii);
	}
	cout << "Chi2 data to NEUT SF: " << tempChi2/19. << endl;

	for(int ii=0; ii<19; ii++){
		chiVec1(ii) = (xsecStat[ii]->GetMean() - nuisNeutRes->GetBinContent(ii+1))/xsecStat[ii]->GetMean();
	}
	chiVec2 = FitCov*chiVec1;
	tempChi2=0;
	for(int ii=0; ii<19; ii++){
		tempChi2 += chiVec1(ii)*chiVec2(ii);
	}
	cout << "Chi2 data to NUISANCE NEUT: " << tempChi2/19. << endl;

	//TODO: print results from DrawXsecFast for Latex
	
	//CovMat:
	cout << "======Cov Mat===1==" << endl;
		cout << "&";
		for(int i=0; i<6; i++){
			cout << i+1;
			if(i<5) cout << "&";
			else cout << "\\\\" << endl;
		}
		for(int i=0; i<19; i++){
			cout << i+1 << "&";
			for(int j=0; j<6; j++){
				cout << std::setprecision(4) << FitCov(i,j);
				if(j<5) cout << "&";
				else cout << "\\\\" << endl;
			}
		}

		cout << "======Cov Mat===2==" << endl;
		cout << "&";
		for(int i=6; i<12; i++){
			cout << i+1;
			if(i<11) cout << "&";
			else cout << "\\\\" << endl;
		}
		for(int i=0; i<19; i++){
			cout << i+1 << "&";
			for(int j=6; j<12; j++){
				cout << std::setprecision(4) << FitCov(i,j);
				if(j<11) cout << "&";
				else cout << "\\\\" << endl;
			}
		}

		cout << "======Cov Mat===3==" << endl;
		cout << "&";
		for(int i=12; i<19; i++){
			cout << i+1;
			if(i<18) cout << "&";
			else cout << "\\\\" << endl;
		}
		for(int i=0; i<19; i++){
			cout << i+1 << "&";
			for(int j=12; j<19; j++){
				cout << std::setprecision(4) << FitCov(i,j);
				if(j<18) cout << "&";
				else cout << "\\\\" << endl;
			}
		}



}
