#include <iostream>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void DrawXsec(Float_t** nData, Float_t** nSel, Float_t** nGen, Float_t* binWidth, Float_t* intFlux, Float_t* nTargets, Float_t nTargetsNomMC, int nToys){
	//calc xsec
	suffStat** xsecStat = new suffStat*[19];
	suffStat** xsecStatNEUT = new suffStat*[19];
	for(int ii=0; ii<19; ii++){
		xsecStat[ii] = new suffStat(1e-39);
		xsecStatNEUT[ii] = new suffStat(1e-39);
	}

	for(int iToy=0; iToy<nToys; iToy++){
		for(int ii=0; ii<19; ii++){
			Float_t xsec = nData[iToy][ii]/(nSel[iToy][ii]/nGen[iToy][ii]);
			xsec*=1.0/binWidth[ii];
			xsec*=1.0/intFlux[iToy];
			xsec*=1.0/nTargets[iToy];
			xsecStat[ii]->Fill(xsec);
			xsec = nGen[iToy][ii]/(binWidth[ii]*intFlux[iToy]);
			xsec*=1.0/nTargetsNomMC;
			xsecStatNEUT[ii]->Fill(xsec);
		}
	}

	//calc final cov matrix 

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

	leg = new TLegend(0.1,0.1,0.9,0.9);
	leg->AddEntry(CosHist1,"Post Fit MC","pel");
	leg->AddEntry(CosHistNEUT1,"Pre Fit MC (NEUT)","l");

	CosHist1->Draw("PE0");
	CosHistNEUT1->Draw("same");
	c->Print("./plots/CosByMom1.pdf");

	c = new TCanvas;
	leg->Draw();
	c->Print("./plots/CosByMomLegend.pdf");

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

	CosHist2->Draw("PE0");
	CosHistNEUT2->Draw("same");
	c->Print("./plots/CosByMom2.pdf");

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

	CosHist3->Draw("PE0");
	CosHistNEUT3->Draw("same");
	c->Print("./plots/CosByMom3.pdf");

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

	CosHist4->Draw("PE0");
	CosHistNEUT4->Draw("same");
	c->Print("./plots/CosByMom4.pdf");

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

	CosHist5->Draw("PE0");
	CosHistNEUT5->Draw("same");
	c->Print("./plots/CosByMom5.pdf");

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

	CosHist6->Draw("PE0");
	CosHistNEUT6->Draw("same");
	c->Print("./plots/CosByMom6.pdf");

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

	CosHist7->Draw("PE0");
	CosHistNEUT7->Draw("same");
	c->Print("./plots/CosByMom7.pdf");

}
