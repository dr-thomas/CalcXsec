#include <iostream>
#include <iomanip>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>

#include "TColor.h"
#include "/Volumes/ThomasDrive/p0dCCAnalysis/FitResults/Macros/Phil_root-style-files-master/palettes.C"
#include "TPaletteAxis.h"


void DrawNuisComp(Float_t** nData, Float_t** nSel, Float_t** nGen, Float_t** nGenSF, Float_t* binWidth, Float_t* intFlux, Float_t* nTargets, Float_t nTargetsNomMC, int nToys){
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

	suffStat* sbStat = new suffStat(1e-39);
	suffStat* sbStatN = new suffStat(1e-39);
	suffStat* sbStatNSF = new suffStat(1e-39);
	for(int ii=0; ii<nToys; ii++){
		sbStat->Fill(singleBin[ii]);
		sbStatN->Fill(singleBinN[ii]);
		sbStatNSF->Fill(singleBinNSF[ii]);
	}
	cout << "Single Bin Results:" << endl;
	cout << "Data: " << sbStat->GetMean() << " +- " << sbStat->GetRMS() << endl;
	cout << "NEUT: " << sbStatN->GetMean() << endl;
	cout << "NEUT SF: " << sbStatNSF->GetMean() << endl;

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


	TH1F* CosHist1 = new TH1F("CosHist1","",2,CosBinDraw[1]);
	TH1F* CosHist2 = new TH1F("CosHist2","",3,CosBinDraw[2]);
	TH1F* CosHist3 = new TH1F("CosHist3","",3,CosBinDraw[3]);
	TH1F* CosHist4 = new TH1F("CosHist4","",3,CosBinDraw[4]);
	TH1F* CosHist5 = new TH1F("CosHist5","",3,CosBinDraw[5]);
	TH1F* CosHist6 = new TH1F("CosHist6","",3,CosBinDraw[6]);
	TH1F* CosHist7 = new TH1F("CosHist7","",2,CosBinDraw[7]);

	TH1F* CosHistNEUT1 = new TH1F("CosHistNEUT1","",2,CosBinDraw[1]);
	TH1F* CosHistNEUT2 = new TH1F("CosHistNEUT2","",3,CosBinDraw[2]);
	TH1F* CosHistNEUT3 = new TH1F("CosHistNEUT3","",3,CosBinDraw[3]);
	TH1F* CosHistNEUT4 = new TH1F("CosHistNEUT4","",3,CosBinDraw[4]);
	TH1F* CosHistNEUT5 = new TH1F("CosHistNEUT5","",3,CosBinDraw[5]);
	TH1F* CosHistNEUT6 = new TH1F("CosHistNEUT6","",3,CosBinDraw[6]);
	TH1F* CosHistNEUT7 = new TH1F("CosHistNEUT7","",2,CosBinDraw[7]);

	TH1F* CosHistNEUTSF1 = new TH1F("CosHistNEUTSF1","",2,CosBinDraw[1]);
	TH1F* CosHistNEUTSF2 = new TH1F("CosHistNEUTSF2","",3,CosBinDraw[2]);
	TH1F* CosHistNEUTSF3 = new TH1F("CosHistNEUTSF3","",3,CosBinDraw[3]);
	TH1F* CosHistNEUTSF4 = new TH1F("CosHistNEUTSF4","",3,CosBinDraw[4]);
	TH1F* CosHistNEUTSF5 = new TH1F("CosHistNEUTSF5","",3,CosBinDraw[5]);
	TH1F* CosHistNEUTSF6 = new TH1F("CosHistNEUTSF6","",3,CosBinDraw[6]);
	TH1F* CosHistNEUTSF7 = new TH1F("CosHistNEUTSF7","",2,CosBinDraw[7]);


	TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

	CosHist1->SetBinContent(1,xsecStat[0]->GetMean());
	CosHist1->SetBinError(1,xsecStat[0]->GetRMS());
	CosHist1->SetBinContent(2,xsecStat[1]->GetMean());
	CosHist1->SetBinError(2,xsecStat[1]->GetRMS());

	CosHist2->SetBinContent(1,xsecStat[2]->GetMean());
	CosHist2->SetBinError(1,xsecStat[2]->GetRMS());
	CosHist2->SetBinContent(2,xsecStat[3]->GetMean());
	CosHist2->SetBinError(2,xsecStat[3]->GetRMS());
	CosHist2->SetBinContent(3,xsecStat[4]->GetMean());
	CosHist2->SetBinError(3,xsecStat[4]->GetRMS());

	CosHist3->SetBinContent(1,xsecStat[5]->GetMean());
	CosHist3->SetBinError(1,xsecStat[5]->GetRMS());
	CosHist3->SetBinContent(2,xsecStat[6]->GetMean());
	CosHist3->SetBinError(2,xsecStat[6]->GetRMS());
	CosHist3->SetBinContent(3,xsecStat[7]->GetMean());
	CosHist3->SetBinError(3,xsecStat[7]->GetRMS());

	CosHist4->SetBinContent(1,xsecStat[8]->GetMean());
	CosHist4->SetBinError(1,xsecStat[8]->GetRMS());
	CosHist4->SetBinContent(2,xsecStat[9]->GetMean());
	CosHist4->SetBinError(2,xsecStat[9]->GetRMS());
	CosHist4->SetBinContent(3,xsecStat[10]->GetMean());
	CosHist4->SetBinError(3,xsecStat[10]->GetRMS());

	CosHist5->SetBinContent(1,xsecStat[11]->GetMean());
	CosHist5->SetBinError(1,xsecStat[11]->GetRMS());
	CosHist5->SetBinContent(2,xsecStat[12]->GetMean());
	CosHist5->SetBinError(2,xsecStat[12]->GetRMS());
	CosHist5->SetBinContent(3,xsecStat[13]->GetMean());
	CosHist5->SetBinError(3,xsecStat[13]->GetRMS());

	CosHist6->SetBinContent(1,xsecStat[14]->GetMean());
	CosHist6->SetBinError(1,xsecStat[14]->GetRMS());
	CosHist6->SetBinContent(2,xsecStat[15]->GetMean());
	CosHist6->SetBinError(2,xsecStat[15]->GetRMS());
	CosHist6->SetBinContent(3,xsecStat[16]->GetMean());
	CosHist6->SetBinError(3,xsecStat[16]->GetRMS());

	CosHist7->SetBinContent(1,xsecStat[17]->GetMean());
	CosHist7->SetBinError(1,xsecStat[17]->GetRMS());
	CosHist7->SetBinContent(2,xsecStat[18]->GetMean());
	CosHist7->SetBinError(2,xsecStat[18]->GetRMS());

	//NEUT
	CosHistNEUT1->SetBinContent(1,xsecStatNEUT[0]->GetMean());
	CosHistNEUT1->SetBinContent(2,xsecStatNEUT[1]->GetMean());

	CosHistNEUT2->SetBinContent(1,xsecStatNEUT[2]->GetMean());
	CosHistNEUT2->SetBinContent(2,xsecStatNEUT[3]->GetMean());
	CosHistNEUT2->SetBinContent(3,xsecStatNEUT[4]->GetMean());

	CosHistNEUT3->SetBinContent(1,xsecStatNEUT[5]->GetMean());
	CosHistNEUT3->SetBinContent(2,xsecStatNEUT[6]->GetMean());
	CosHistNEUT3->SetBinContent(3,xsecStatNEUT[7]->GetMean());

	CosHistNEUT4->SetBinContent(1,xsecStatNEUT[8]->GetMean());
	CosHistNEUT4->SetBinContent(2,xsecStatNEUT[9]->GetMean());
	CosHistNEUT4->SetBinContent(3,xsecStatNEUT[10]->GetMean());

	CosHistNEUT5->SetBinContent(1,xsecStatNEUT[11]->GetMean());
	CosHistNEUT5->SetBinContent(2,xsecStatNEUT[12]->GetMean());
	CosHistNEUT5->SetBinContent(3,xsecStatNEUT[13]->GetMean());

	CosHistNEUT6->SetBinContent(1,xsecStatNEUT[14]->GetMean());
	CosHistNEUT6->SetBinContent(2,xsecStatNEUT[15]->GetMean());
	CosHistNEUT6->SetBinContent(3,xsecStatNEUT[16]->GetMean());

	CosHistNEUT7->SetBinContent(1,xsecStatNEUT[17]->GetMean());
	CosHistNEUT7->SetBinContent(2,xsecStatNEUT[18]->GetMean());

	//NEUTSF
	CosHistNEUTSF1->SetBinContent(1,xsecStatNEUTSF[0]->GetMean());
	CosHistNEUTSF1->SetBinContent(2,xsecStatNEUTSF[1]->GetMean());

	CosHistNEUTSF2->SetBinContent(1,xsecStatNEUTSF[2]->GetMean());
	CosHistNEUTSF2->SetBinContent(2,xsecStatNEUTSF[3]->GetMean());
	CosHistNEUTSF2->SetBinContent(3,xsecStatNEUTSF[4]->GetMean());

	CosHistNEUTSF3->SetBinContent(1,xsecStatNEUTSF[5]->GetMean());
	CosHistNEUTSF3->SetBinContent(2,xsecStatNEUTSF[6]->GetMean());
	CosHistNEUTSF3->SetBinContent(3,xsecStatNEUTSF[7]->GetMean());

	CosHistNEUTSF4->SetBinContent(1,xsecStatNEUTSF[8]->GetMean());
	CosHistNEUTSF4->SetBinContent(2,xsecStatNEUTSF[9]->GetMean());
	CosHistNEUTSF4->SetBinContent(3,xsecStatNEUTSF[10]->GetMean());

	CosHistNEUTSF5->SetBinContent(1,xsecStatNEUTSF[11]->GetMean());
	CosHistNEUTSF5->SetBinContent(2,xsecStatNEUTSF[12]->GetMean());
	CosHistNEUTSF5->SetBinContent(3,xsecStatNEUTSF[13]->GetMean());

	CosHistNEUTSF6->SetBinContent(1,xsecStatNEUTSF[14]->GetMean());
	CosHistNEUTSF6->SetBinContent(2,xsecStatNEUTSF[15]->GetMean());
	CosHistNEUTSF6->SetBinContent(3,xsecStatNEUTSF[16]->GetMean());

	CosHistNEUTSF7->SetBinContent(1,xsecStatNEUTSF[17]->GetMean());
	CosHistNEUTSF7->SetBinContent(2,xsecStatNEUTSF[18]->GetMean());

	//NUISANCE
	TFile* inFnuisance = new TFile("rootFiles/nuiscompAnuCiroScaled.root", "OPEN");
	TH1F* nuisNeutRes = (TH1F*) inFnuisance->Get("T2K_CC0pi_XSec_2DPcos_anu_P0D_MC");

	TH1F** CosHistsNuisNEUT = new TH1F*[7];
	for(int ii=0; ii<7; ii++){
		int tempNbins=0;
		if(ii==0 || ii==6){
			tempNbins=2;
		} else {
			tempNbins=3;
		}
		TString labelStr = "CosHistsNuisNEUT";
		labelStr += ii;
		CosHistsNuisNEUT[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii]);
	}


	// ------------------  going to re write all this -------------------//
	//Cos by momentum slice 
	TString PBinsStr[9]={"0","400","530","670","800","1000","1380","2010","3410"};
	TString CosTitleStr;

	gStyle->SetOptStat(0);

	//1
	TCanvas* c = new TCanvas;
	gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
	CosHist1->GetYaxis()->SetTitleOffset(1.4);


	CosTitleStr="";
	CosTitleStr+=PBinsStr[1];
	CosTitleStr+=" < P_{#mu} < ";
	CosTitleStr+=PBinsStr[2];

	CosHist1->GetXaxis()->SetTitle("sel. #mu^{+} cos#theta");
	CosHist1->GetYaxis()->SetTitle("#frac{d#sigma}{dpd(cos#theta)} (cm^{2}/MeV/H_{2}O molecule)");
	CosHist1->SetMarkerStyle(8);
	CosHist1->SetMarkerSize(1);
	CosHist1->SetMinimum(0);
	CosHist1->SetTitle(CosTitleStr);

	CosHistNEUT1->SetLineColor(kRed+2);
	CosHistNEUT1->SetLineStyle(7);

	CosHistNEUTSF1->SetLineColor(kGreen+2);
	CosHistNEUTSF1->SetLineStyle(7);


	leg = new TLegend(0.1,0.1,0.9,0.9);
	leg->AddEntry(CosHist1,"Data","pel");
	leg->AddEntry(CosHistNEUT1,"NEUT RFG+RPA","l");
	leg->AddEntry(CosHistNEUTSF1,"NEUT SF","l");

	CosHist1->Draw("PE0");
	CosHistNEUT1->Draw("same");
	CosHistNEUTSF1->Draw("same");
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/CosByMom1.pdf");

	c = new TCanvas;
	leg->Draw();
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/CosByMomLegend.pdf");

	//2
	c = new TCanvas;
	gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
	CosHist2->GetYaxis()->SetTitleOffset(1.4);


	CosTitleStr="";
	CosTitleStr+=PBinsStr[2];
	CosTitleStr+=" < P_{#mu} < ";
	CosTitleStr+=PBinsStr[3];

	CosHist2->GetXaxis()->SetTitle("sel. #mu^{+} cos#theta");
	CosHist2->GetYaxis()->SetTitle("#frac{d#sigma}{dpd(cos#theta)} (cm^{2}/MeV/H_{2}O molecule)");
	CosHist2->SetMarkerStyle(8);
	CosHist2->SetMarkerSize(1);
	CosHist2->SetMinimum(0);
	CosHist2->SetTitle(CosTitleStr);

	CosHistNEUT2->SetLineColor(kRed+2);
	CosHistNEUT2->SetLineStyle(7);

	CosHistNEUTSF2->SetLineColor(kGreen+2);
	CosHistNEUTSF2->SetLineStyle(7);


	CosHist2->Draw("PE0");
	CosHistNEUT2->Draw("same");
	CosHistNEUTSF2->Draw("same");
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/CosByMom2.pdf");

	//3
	c = new TCanvas;
	gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
	CosHist3->GetYaxis()->SetTitleOffset(1.4);


	CosTitleStr="";
	CosTitleStr+=PBinsStr[3];
	CosTitleStr+=" < P_{#mu} < ";
	CosTitleStr+=PBinsStr[4];

	CosHist3->GetXaxis()->SetTitle("sel. #mu^{+} cos#theta");
	CosHist3->GetYaxis()->SetTitle("#frac{d#sigma}{dpd(cos#theta)} (cm^{2}/MeV/H_{2}O molecule)");
	CosHist3->SetMarkerStyle(8);
	CosHist3->SetMarkerSize(1);
	CosHist3->SetMinimum(0);
	CosHist3->SetTitle(CosTitleStr);

	CosHistNEUT3->SetLineColor(kRed+2);
	CosHistNEUT3->SetLineStyle(7);

	CosHistNEUTSF3->SetLineColor(kGreen+2);
	CosHistNEUTSF3->SetLineStyle(7);

	CosHist3->Draw("PE0");
	CosHistNEUT3->Draw("same");
	CosHistNEUTSF3->Draw("same");
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/CosByMom3.pdf");

	//4
	c = new TCanvas;
	gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
	CosHist4->GetYaxis()->SetTitleOffset(1.4);


	CosTitleStr="";
	CosTitleStr+=PBinsStr[4];
	CosTitleStr+=" < P_{#mu} < ";
	CosTitleStr+=PBinsStr[5];

	CosHist4->GetXaxis()->SetTitle("sel. #mu^{+} cos#theta");
	CosHist4->GetYaxis()->SetTitle("#frac{d#sigma}{dpd(cos#theta)} (cm^{2}/MeV/H_{2}O molecule)");
	CosHist4->SetMarkerStyle(8);
	CosHist4->SetMarkerSize(1);
	CosHist4->SetMinimum(0);
	CosHist4->SetTitle(CosTitleStr);

	CosHistNEUT4->SetLineColor(kRed+2);
	CosHistNEUT4->SetLineStyle(7);

	CosHistNEUTSF4->SetLineColor(kGreen+2);
	CosHistNEUTSF4->SetLineStyle(7);

	CosHist4->Draw("PE0");
	CosHistNEUT4->Draw("same");
	CosHistNEUTSF4->Draw("same");
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/CosByMom4.pdf");

	//5
	c = new TCanvas;
	gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
	CosHist5->GetYaxis()->SetTitleOffset(1.4);


	CosTitleStr="";
	CosTitleStr+=PBinsStr[5];
	CosTitleStr+=" < P_{#mu} < ";
	CosTitleStr+=PBinsStr[6];

	CosHist5->GetXaxis()->SetTitle("sel. #mu^{+} cos#theta");
	CosHist5->GetYaxis()->SetTitle("#frac{d#sigma}{dpd(cos#theta)} (cm^{2}/MeV/H_{2}O molecule)");
	CosHist5->SetMarkerStyle(8);
	CosHist5->SetMarkerSize(1);
	CosHist5->SetMinimum(0);
	CosHist5->SetTitle(CosTitleStr);

	CosHistNEUT5->SetLineColor(kRed+2);
	CosHistNEUT5->SetLineStyle(7);

	CosHistNEUTSF5->SetLineColor(kGreen+2);
	CosHistNEUTSF5->SetLineStyle(7);

	CosHist5->Draw("PE0");
	CosHistNEUT5->Draw("same");
	CosHistNEUTSF5->Draw("same");
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/CosByMom5.pdf");

	//6
	c = new TCanvas;
	gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
	CosHist6->GetYaxis()->SetTitleOffset(1.4);


	CosTitleStr="";
	CosTitleStr+=PBinsStr[6];
	CosTitleStr+=" < P_{#mu} < ";
	CosTitleStr+=PBinsStr[7];

	CosHist6->GetXaxis()->SetTitle("sel. #mu^{+} cos#theta");
	CosHist6->GetYaxis()->SetTitle("#frac{d#sigma}{dpd(cos#theta)} (cm^{2}/MeV/H_{2}O molecule)");
	CosHist6->SetMarkerStyle(8);
	CosHist6->SetMarkerSize(1);
	CosHist6->SetMinimum(0);
	CosHist6->SetMaximum(1.8*CosHist6->GetMaximum());
	CosHist6->SetTitle(CosTitleStr);

	CosHistNEUT6->SetLineColor(kRed+2);
	CosHistNEUT6->SetLineStyle(7);

	CosHistNEUTSF6->SetLineColor(kGreen+2);
	CosHistNEUTSF6->SetLineStyle(7);

	CosHist6->Draw("PE0");
	CosHistNEUT6->Draw("same");
	CosHistNEUTSF6->Draw("same");
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/CosByMom6.pdf");

	//7
	c = new TCanvas;
	gPad->SetLeftMargin(0.15); //left margin is 15 per cent of the pad width
	CosHist7->GetYaxis()->SetTitleOffset(1.4);


	CosTitleStr="";
	CosTitleStr+=PBinsStr[7];
	CosTitleStr+=" < P_{#mu} < ";
	CosTitleStr+=PBinsStr[8];

	CosHist7->GetXaxis()->SetTitle("sel. #mu^{+} cos#theta");
	CosHist7->GetYaxis()->SetTitle("#frac{d#sigma}{dpd(cos#theta)} (cm^{2}/MeV/H_{2}O molecule)");
	CosHist7->SetMarkerStyle(8);
	CosHist7->SetMarkerSize(1);
	CosHist7->SetMinimum(0);
	CosHist7->SetTitle(CosTitleStr);

	CosHistNEUT7->SetLineColor(kRed+2);
	CosHistNEUT7->SetLineStyle(7);

	CosHistNEUTSF7->SetLineColor(kGreen+2);
	CosHistNEUTSF7->SetLineStyle(7);

	CosHist7->Draw("PE0");
	CosHistNEUT7->Draw("same");
	CosHistNEUTSF7->Draw("same");
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/CosByMom7.pdf");

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
	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/MomOverlay.pdf");

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

	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/xsecCovMatrix.pdf");
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

	c->Print("/Users/thomascampbell/Documents/Research/NewCC0piAnalysisPlots/Xsec/xsecCovMatrixMom.pdf");
	delete c;

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
