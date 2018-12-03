#include <iostream>
#include <iomanip>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector.h>

#include "TLatex.h"
#include "TPad.h"

#include "TPaletteAxis.h"

#include "../util/suffstat.hxx"
#include "../util/suffstat.cxx"


void DrawCallumXsec(Float_t** nData, Float_t** nSel, Float_t** nGen, Float_t** nGenSF, Float_t* binWidth, Float_t* intFlux, Float_t* nTargets, Float_t nTargetsNomMC, int nToys){

	//calc xsec
	suffStat** xsecStat = new suffStat*[19];
	suffStat** xsecStatNEUT = new suffStat*[19];
	suffStat** xsecStatNEUTSF = new suffStat*[19];
	for(int ii=0; ii<19; ii++){
		xsecStat[ii] = new suffStat(1e-39);
		xsecStatNEUT[ii] = new suffStat(1e-39);
		xsecStatNEUTSF[ii] = new suffStat(1e-39);
	}

	for(int iToy=0; iToy<nToys; iToy++){
		Float_t xsecMom=0.;
		Float_t xsecMomN=0.;
		Float_t xsecMomNSF=0.;
		for(int ii=0; ii<19; ii++){
			Float_t xsec = nData[iToy][ii]/(nSel[iToy][ii]/nGen[iToy][ii]);
			Float_t xsecN = nGen[iToy][ii]/(binWidth[ii]*intFlux[iToy]);
			Float_t xsecNSF = nGenSF[iToy][ii]/(binWidth[ii]*intFlux[iToy]);
			xsec*=1.0/binWidth[ii];
			xsec*=1.0/intFlux[iToy];
			xsec*=1.0/nTargets[iToy];
			xsecN*=1.0/nTargetsNomMC;
			xsecNSF*=1.0/nTargetsNomMC;
			xsecStat[ii]->Fill(xsec);
			xsecStatNEUT[ii]->Fill(xsecN);
			xsecStatNEUTSF[ii]->Fill(xsecNSF);
		}
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

	//Cos by momentum xsec slice 
	TString PBinsStr[9]={"0","400","530","670","800","1000","1380","2010","3410"};

	TCanvas* c;

	gStyle->SetOptStat(0);

	gStyle->SetLineWidth(1);
	gStyle->SetLabelSize(0.09, "xyzt");
	gStyle->SetTitleSize(0.10, "xyzt");
	gStyle->SetTitleSize(0.15, "t");
	gStyle->SetTitleX(0.55);
	gStyle->SetTitleY(0.24);
	gStyle->SetTitleW(0.8);
	gStyle->SetTitleH(0.11);
	gStyle->SetOptTitle(1);
	gStyle->SetLineWidth(2);

	gStyle->SetTextFont(42);

	TGaxis::SetMaxDigits(3);

	c = new TCanvas("can", "can", 2800, 1500);

	TLatex xText(0.48,0.04,"cos#theta_{#mu}");
	TLatex yText(0.058,0.2,"#frac{#it{#it{d}}^{2}#sigma}{dp_{#mu}d(cos#theta_{#mu})} (#times 10^{-41} cm^{2}/H_{2}O)");
	yText.SetTextAngle(90);
	xText.Draw();
	yText.Draw();

	gPad->SetRightMargin(0);
	gPad->SetLeftMargin(0);
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0);

    TPad* subPad = new TPad("pad", "pad", 0.1, 0.1, 1.0, 1.0);
    subPad->Divide(4, 2, 0.0001, 0.0001);
    subPad->Draw();

	TString CosTitleStr;

	for(int ii=0; ii<7; ii++){

        subPad->cd(ii+1);

		CosTitleStr="";
		CosTitleStr+=PBinsStr[ii+1];
		CosTitleStr+=" < P_{#mu} < ";
		CosTitleStr+=PBinsStr[ii+2];

		CosHists[ii]->SetTitle(CosTitleStr);

		CosHists[ii]->Scale(1e41);
		CosHistsNEUT[ii]->Scale(1e41);

		float maxval = CosHists[ii]->GetMaximum();
		maxval *= 1.5;
		CosHists[ii]->SetMaximum(maxval);
		CosHists[ii]->SetMinimum(0);

		CosHists[ii]->GetYaxis()->SetLabelOffset(0.015);
		CosHists[ii]->GetXaxis()->SetLabelOffset(0.015);
		CosHists[ii]->GetXaxis()->SetNdivisions(505);
		CosHists[ii]->GetYaxis()->SetNdivisions(505);
		CosHists[ii]->GetXaxis()->SetLabelSize(0.07);
		CosHists[ii]->GetYaxis()->SetLabelSize(0.07);

		CosHists[ii]->SetLineColor(kBlack);
		CosHists[ii]->SetLineWidth(3);

		CosHistsNEUT[ii]->SetLineColor(kRed);
		CosHistsNEUT[ii]->SetLineWidth(3);
		CosHistsNEUT[ii]->SetLineStyle(7);

		CosHists[ii]->Draw("PE0");
		CosHistsNEUT[ii]->Draw("same hist");

        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.10);  
        gPad->SetRightMargin(0.05);
        gPad->SetLeftMargin(0.10);
        gPad->SetTickx();
        gPad->SetTicky();
		gPad->Update();
	}

    subPad->cd(8);
    TLegend* leg  = new TLegend(0.1, 0.1, 0.95, 0.9, "", "NDC");
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.09);
    leg->SetLineWidth(0);
    leg->AddEntry(CosHists[0], "DATA", "le");
    leg->AddEntry(CosHistsNEUT[0], "NEUT v5.3.3", "l");
    leg->Draw();
    gPad->Update();
    
    subPad->Update();
    c->Update();
    c->SaveAs("~/Desktop/tn328_updated.png");
}
