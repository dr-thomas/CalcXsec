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

#include "TMatrixF.h"


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

	//calc xsec relative covariance matrix
	Float_t** xsecCovMat = new Float_t*[19];
	for(int ii=0; ii<19; ii++){
		xsecCovMat[ii] = new Float_t[19];
	}
	for(int ii=0; ii<19; ii++){
		for(int jj=0; jj<19; jj++){
			xsecCovMat[ii][jj] = 0.;
		}
	}

	for(int iToy=0; iToy<nToys; iToy++){
		for(int ii=0; ii<19; ii++){
			Float_t xsec_i = nData[iToy][ii]/(nSel[iToy][ii]/nGen[iToy][ii]);
			xsec_i*=1.0/binWidth[ii];
			xsec_i*=1.0/intFlux[iToy];
			xsec_i*=1.0/nTargets[iToy];
			for(int jj=0; jj<19; jj++){
				Float_t xsec_j = nData[iToy][jj]/(nSel[iToy][jj]/nGen[iToy][jj]);
				xsec_j*=1.0/binWidth[jj];
				xsec_j*=1.0/intFlux[iToy];
				xsec_j*=1.0/nTargets[iToy];

				xsecCovMat[ii][jj] += (xsec_i-xsecStat[ii]->GetMean())/xsecStat[ii]->GetMean()*(xsec_j-xsecStat[jj]->GetMean())/xsecStat[jj]->GetMean();
			}
		}
	}
	for(int ii=0; ii<19; ii++){
		for(int jj=0; jj<19; jj++){
			xsecCovMat[ii][jj] *= (Float_t)1./nToys;
		}
	}

	TFile* outFDR = new TFile("DataRelease.root", "RECREATE");
	//TODO: print central values and cov matrix to text file, root file
	TH1F* xsecDataRelease = new TH1F("xsecDataRelease", "cross_section_central_values", 19,0,19);
	TH2F* covDataRelease = new TH2F("covDataRelease", "cross_section_covariance_matrix", 19,0,19,19,0,19);
	
	cout << "bin i: xsec[i]" << endl;
	for(int ii=0; ii<19; ii++){
		xsecDataRelease->SetBinContent(ii+1, xsecStat[ii]->GetMean());
		cout << ii+1 << ": " << xsecStat[ii]->GetMean() << endl;
	}
	cout << "bin i, bin j: cov[i][j]" << endl;
	for(int ii=0; ii<19; ii++){
		for(int jj=0; jj<19; jj++){
			covDataRelease->SetBinContent(ii+1,jj+1, xsecCovMat[ii][jj]);
			cout << ii+1 << ", " << jj+1 << ": " << xsecCovMat[ii][jj] << endl;
		}
	}

	TMatrixF xsecCov(19,19);
	for(int ii=0; ii<19; ii++){
		for(int jj=0; jj<19; jj++){
			xsecCov(ii,jj) = xsecCovMat[ii][jj];
		}
	}
	xsecDataRelease->Write();
	covDataRelease->Write();
	outFDR->Write();
	outFDR->Close();


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
	TFile* inFnuisanceNeut5_4_0 = new TFile("~/Downloads/neut_5.4.0_T2K-P0D_H2O.root", "OPEN");
	TH1F* nuisNeutResNeut5_4_0 = (TH1F*) inFnuisanceNeut5_4_0->Get("T2K_CC0pi_XSec_2DPcos_anu_P0D_MC");
	TFile* inFnuisanceNeut5_4_1 = new TFile("~/Downloads/neut_5.4.1_T2K-P0D_AntiNuMuCC0pi_H2O.root", "OPEN");
	TH1F* nuisNeutResNeut5_4_1 = (TH1F*) inFnuisanceNeut5_4_1->Get("T2K_CC0pi_XSec_2DPcos_anu_P0D_MC");
	TFile* inFnuisanceGENIE = new TFile("~/Downloads/genie_r-2-12-10_T2K-P0D_AntiNuMuCC0pi_H2O.root", "OPEN");
	TH1F* nuisNeutResGENIE = (TH1F*) inFnuisanceGENIE->Get("T2K_CC0pi_XSec_2DPcos_anu_P0D_MC");
	//TFile* inFnuisanceNuWro = new TFile("~/Downloads/nuwro_18.02.1_T2K-P0D_AntiNuMuCC0pi_H2O.root", "OPEN");
	TFile* inFnuisanceNuWro = new TFile("~/Downloads/nuwro_18.02.1_T2K-P0D_AntiNuMuCC0pi_H2O-1.root", "OPEN");
	TH1F* nuisNeutResNuWro = (TH1F*) inFnuisanceNuWro->Get("T2K_CC0pi_XSec_2DPcos_anu_P0D_MC");

	//TFile* inFnuisanceNuWro1 = new TFile("~/Downloads/nuwro_18.02.1_T2K-P0D_AntiNuMuCC0pi_H2O-1.root", "OPEN");
	//TFile* inFnuisanceNuWro1 = new TFile("~/Downloads/nuwro_18.02.1_T2K-P0D_AntiNuMuCC0pi_H2O-3.root", "OPEN");
	TFile* inFnuisanceNuWro1 = new TFile("~/Downloads/nuwro_18.02.1_T2K-P0D_AntiNuMuCC0pi_H2O-4.root", "OPEN");
	TH1F* nuisNeutResNuWro1 = (TH1F*) inFnuisanceNuWro1->Get("T2K_CC0pi_XSec_2DPcos_anu_P0D_MC");

	//chi2s
	xsecCov.Invert();

	TVector resNeut(19);
	for(int ii=0; ii<19; ii++){
		resNeut(ii) = (nuisNeutResNeut5_4_1->GetBinContent(ii+1) - xsecStat[ii]->GetMean())/xsecStat[ii]->GetMean();
	}
	Float_t chi2resNeut = 0.;
	TVector resNeut1 = xsecCov*resNeut;
	for(int ii=0; ii<19; ii++){
		chi2resNeut += resNeut(ii)*resNeut1(ii);
	}
	cout << "chi2 to NEUT: " << chi2resNeut << endl;

	TVector resGenie(19);
	for(int ii=0; ii<19; ii++){
		resGenie(ii) = (nuisNeutResGENIE->GetBinContent(ii+1) - xsecStat[ii]->GetMean())/xsecStat[ii]->GetMean();
	}
	Float_t chi2resGenie = 0.;
	TVector resGenie1 = xsecCov*resGenie;
	for(int ii=0; ii<19; ii++){
		chi2resGenie += resGenie(ii)*resGenie1(ii);
	}
	cout << "chi2 to GENIE: " << chi2resGenie << endl;

	TVector resNuWro(19);
	for(int ii=0; ii<19; ii++){
		resNuWro(ii) = (nuisNeutResNuWro1->GetBinContent(ii+1) - xsecStat[ii]->GetMean())/xsecStat[ii]->GetMean();
	}
	Float_t chi2resNuWro = 0.;
	TVector resNuWro1 = xsecCov*resNuWro;
	for(int ii=0; ii<19; ii++){
		chi2resNuWro += resNuWro(ii)*resNuWro1(ii);
	}
	cout << "chi2 to NuWro: " << chi2resNuWro << endl;



	TH1F** CosHists = new TH1F*[7];
	TH1F** CosHistsNEUT = new TH1F*[7];
	TH1F** CosHistsNEUTSF = new TH1F*[7];
	TH1F** CosHistsNuisNEUT = new TH1F*[7];
	TH1F** CosHistsNuisNEUT5_4_1 = new TH1F*[7];
	TH1F** CosHistsNuisNEUT5_4_0 = new TH1F*[7];
	TH1F** CosHistsNuisGENIE = new TH1F*[7];
	TH1F** CosHistsNuisNuWro = new TH1F*[7];
	TH1F** CosHistsNuisNuWro1 = new TH1F*[7];

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
		labelStr = "CosHistsNuisNEUT5_4_1";
		labelStr += ii;
		CosHistsNuisNEUT5_4_1[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);
		labelStr = "CosHistsNuisNEUT5_4_0";
		labelStr += ii;
		CosHistsNuisNEUT5_4_0[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);
		labelStr = "CosHistsNuisGENIE";
		labelStr += ii;
		CosHistsNuisGENIE[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);
		labelStr = "CosHistsNuisNuWro";
		labelStr += ii;
		CosHistsNuisNuWro[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);
		labelStr = "CosHistsNuisNuWro1";
		labelStr += ii;
		CosHistsNuisNuWro1[ii] = new TH1F(labelStr,"",tempNbins,CosBinDraw[ii+1]);

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
			CosHistsNuisNEUT5_4_1[ii]->SetBinContent(jj+1, nuisNeutResNeut5_4_1->GetBinContent(nFilled+1));
			CosHistsNuisNEUT5_4_0[ii]->SetBinContent(jj+1, nuisNeutResNeut5_4_0->GetBinContent(nFilled+1));
			CosHistsNuisGENIE[ii]->SetBinContent(jj+1, nuisNeutResGENIE->GetBinContent(nFilled+1));
			CosHistsNuisNuWro[ii]->SetBinContent(jj+1, nuisNeutResNuWro->GetBinContent(nFilled+1));
			CosHistsNuisNuWro1[ii]->SetBinContent(jj+1, nuisNeutResNuWro1->GetBinContent(nFilled+1));
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
	gStyle->SetTitleY(0.23);
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

		CosHistsNuisNEUT[ii]->Scale(1e41);
		CosHistsNuisNEUT[ii]->SetLineColor(kRed);
		CosHistsNuisNEUT[ii]->SetLineWidth(3);
		CosHistsNuisNEUT[ii]->SetLineStyle(7);

		CosHistsNuisNEUT5_4_1[ii]->Scale(1e41);
		CosHistsNuisNEUT5_4_1[ii]->SetLineColor(kBlue);
		CosHistsNuisNEUT5_4_1[ii]->SetLineWidth(3);
		CosHistsNuisNEUT5_4_1[ii]->SetLineStyle(8);

		CosHistsNuisNEUT5_4_0[ii]->Scale(1e41);
		CosHistsNuisNEUT5_4_0[ii]->SetLineColor(kRed);
		CosHistsNuisNEUT5_4_0[ii]->SetLineWidth(3);
		CosHistsNuisNEUT5_4_0[ii]->SetLineStyle(8);

		CosHistsNuisGENIE[ii]->Scale(1e41);
		CosHistsNuisGENIE[ii]->SetLineColor(kGreen);
		CosHistsNuisGENIE[ii]->SetLineWidth(3);
		CosHistsNuisGENIE[ii]->SetLineStyle(9);

		CosHistsNuisNuWro[ii]->Scale(1e41);
		CosHistsNuisNuWro[ii]->SetLineColor(kRed);
		CosHistsNuisNuWro[ii]->SetLineWidth(3);
		CosHistsNuisNuWro[ii]->SetLineStyle(7);

		CosHistsNuisNuWro1[ii]->Scale(1e41);
		//CosHistsNuisNuWro1[ii]->SetLineColor(kBlue);
		CosHistsNuisNuWro1[ii]->SetLineColor(kRed);
		CosHistsNuisNuWro1[ii]->SetLineWidth(3);
		CosHistsNuisNuWro1[ii]->SetLineStyle(7);

		float maxval = CosHistsNuisNuWro[ii]->GetMaximum();
		maxval *= 1.2;
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

		//CosHists[ii]->Draw("PE0");
		//CosHistsNuisNEUT[ii]->Draw("same hist");
		CosHistsNuisNEUT5_4_0[ii]->Draw("same hist");
		CosHistsNuisNEUT5_4_1[ii]->Draw("same hist");
		//CosHistsNuisGENIE[ii]->Draw("same hist");
		//CosHistsNuisNuWro[ii]->Draw("same hist");
		//CosHistsNuisNuWro1[ii]->Draw("same hist");

        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.10);  
        gPad->SetRightMargin(0.05);
        gPad->SetLeftMargin(0.10);
        gPad->SetTickx();
        gPad->SetTicky();
		gPad->Update();
	}

    subPad->cd(8);
    TLegend* leg  = new TLegend(0., 0., 0.7, 0.9, "", "NDC");
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.09);
    leg->SetLineWidth(0);
    leg->AddEntry(CosHists[0], "DATA", "le");

    //Regularized
    //leg->AddEntry(CosHistsNuisNEUT5_4_1[0], "#splitline{NEUT LFG+2p2h}{#chi^{2}=22.2}", "l");//v5.4.0
    //leg->AddEntry(CosHistsNuisGENIE[0], "#splitline{GENIE BRRFG}{#chi^{2}=26.0}", "l");//v2.12.10
    //leg->AddEntry(CosHistsNuisNuWro1[0], "#splitline{NuWro LFG+2p2h}{#chi^{2}=16.8}", "l");//v18.02.1
	
    //Unregularized
    leg->AddEntry(CosHistsNuisNEUT5_4_1[0], "#splitline{NEUT LFG+2p2h}{#chi^{2}=25.1}", "l");//v5.4.0
    leg->AddEntry(CosHistsNuisNEUT5_4_0[0], "#splitline{old NEUT LFG+2p2h}{#chi^{2}=25.1}", "l");//v5.4.0
    //leg->AddEntry(CosHistsNuisGENIE[0], "#splitline{GENIE BRRFG}{#chi^{2}=28.4}", "l");//v2.12.10
    //leg->AddEntry(CosHistsNuisNuWro1[0], "#splitline{NuWro LFG+2p2h}{#chi^{2}=18.4}", "l");//v18.02.1

    leg->Draw();
    gPad->Update();
    
    subPad->Update();
    c->Update();
    c->SaveAs("~/Desktop/tn328_updated.eps");
}
