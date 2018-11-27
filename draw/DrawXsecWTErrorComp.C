#include <iostream>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void DrawXsecWTErrorComp(Float_t** nData, Float_t** nDataB, Float_t** nSel, Float_t** nSelB, Float_t** nGen, Float_t** nGenB, Float_t* binWidth, Float_t* intFlux, Float_t* nTargets, Float_t nTargetsNomMC, int nToys){

	//calc xsec
	suffStat** xsecStat = new suffStat*[19];
	suffStat** xsecStatB = new suffStat*[19];
	suffStat** xsecStatNEUT = new suffStat*[19];
	for(int ii=0; ii<19; ii++){
		xsecStat[ii] = new suffStat(1e-39);
		xsecStatB[ii] = new suffStat(1e-39);
		xsecStatNEUT[ii] = new suffStat(1e-39);
	}
	suffStat** xsecMomStat = new suffStat*[7];
	suffStat** xsecMomStatB = new suffStat*[7];
	suffStat** xsecMomStatNEUT = new suffStat*[7];
	for(int ii=0; ii<7; ii++){
		xsecMomStat[ii] = new suffStat(1e-39);
		xsecMomStatB[ii] = new suffStat(1e-39);
		xsecMomStatNEUT[ii] = new suffStat(1e-39);
	}
	Float_t* singleBin = new Float_t[400];
	Float_t* singleBinB = new Float_t[400];
	Float_t* singleBinN = new Float_t[400];
	for(int ii=0; ii<nToys; ii++){
		singleBin[ii]=0.;
		singleBinB[ii]=0.;
		singleBinN[ii]=0.;
	}

	//TODO: single diff shit is sloppy af, should re-write

	for(int iToy=0; iToy<nToys; iToy++){
		Float_t xsecMom=0.;
		Float_t xsecMomB=0.;
		Float_t xsecMomN=0.;
		int nBinsMom=2;
		int nDrawnP=0;
		int pIndex=0;
		for(int ii=0; ii<19; ii++){
			Float_t xsec = nData[iToy][ii]/(nSel[iToy][ii]/nGen[iToy][ii]);
			Float_t xsecB = nDataB[iToy][ii]/(nSelB[iToy][ii]/nGenB[iToy][ii]);
			xsec*=1.0/binWidth[ii];
			xsec*=1.0/intFlux[iToy];
			xsec*=1.0/nTargets[iToy];
			//xsec*=(1902./2655.)/nTargets[iToy];
			xsecB*=1.0/binWidth[ii];
			xsecB*=1.0/intFlux[iToy];
			xsecB*=(1902./2655.)/nTargets[iToy];
			xsecStat[ii]->Fill(xsec);
			xsecStatB[ii]->Fill(xsecB);
			singleBin[iToy]+=xsec;
			singleBinB[iToy]+=xsecB;
			Float_t xsecN = nGen[iToy][ii]/(binWidth[ii]*intFlux[iToy]);
			xsecN*=1.0/nTargetsNomMC;
			xsecStatNEUT[ii]->Fill(xsecN);
			singleBinN[iToy]+=xsecN;
			//Momentum single (Note: not normalized by bin width)
			xsec*=binWidth[ii];
			xsecB*=binWidth[ii];
			xsecN*=binWidth[ii];
			if(nDrawnP<nBinsMom){
				xsecMom+=xsec;
				xsecMomB+=xsecB;
				xsecMomN+=xsecN;
				nDrawnP++;
				if(ii==18) {
					xsecMomStat[pIndex]->Fill(xsecMom);
					xsecMomStatB[pIndex]->Fill(xsecMomB);
					xsecMomStatNEUT[pIndex]->Fill(xsecMomN);
				}
			} else {
				nDrawnP=0;
				xsecMomStat[pIndex]->Fill(xsecMom);
				xsecMomStatB[pIndex]->Fill(xsecMomB);
				xsecMomStatNEUT[pIndex]->Fill(xsecMomN);
				xsecMom=0.;
				xsecMomB=0.;
				xsecMomN=0.;
				pIndex++;
				if(ii==16) nBinsMom=2;
				else nBinsMom=3;
				//add again
				xsecMom+=xsec;
				xsecMomB+=xsecB;
				xsecMomN+=xsecN;
				nDrawnP++;
			}
		}
	}

	//calc final cov matrix 

	cout << "xsec results in data (x+-error | NEUT): " << endl;
	for(int ii=0; ii<19; ii++){
		cout << xsecStat[ii]->GetMean() << " +- " << xsecStat[ii]->GetRMS();
		cout << " | " << xsecStatNEUT[ii]->GetMean() << endl;
	}

	suffStat* sbStat = new suffStat(1e-39);
	suffStat* sbStatB = new suffStat(1e-39);
	suffStat* sbStatN = new suffStat(1e-39);
	for(int ii=0; ii<nToys; ii++){
		sbStat->Fill(singleBin[ii]);
		sbStatB->Fill(singleBinB[ii]);
		sbStatN->Fill(singleBinN[ii]);
	}
	cout << "Single Bin Results:" << endl;
	cout << "Data: " << sbStat->GetMean() << " +- " << sbStat->GetRMS() << endl;
	cout << "Data B Fit: " << sbStatB->GetMean() << " +- " << sbStatB->GetRMS() << endl;
	cout << "NEUT: " << sbStatN->GetMean() << endl;

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

	TH1F* CosHistB1 = new TH1F("CosHistB1","",2,CosBinDraw[1]);
	TH1F* CosHistB2 = new TH1F("CosHistB2","",3,CosBinDraw[2]);
	TH1F* CosHistB3 = new TH1F("CosHistB3","",3,CosBinDraw[3]);
	TH1F* CosHistB4 = new TH1F("CosHistB4","",3,CosBinDraw[4]);
	TH1F* CosHistB5 = new TH1F("CosHistB5","",3,CosBinDraw[5]);
	TH1F* CosHistB6 = new TH1F("CosHistB6","",3,CosBinDraw[6]);
	TH1F* CosHistB7 = new TH1F("CosHistB7","",2,CosBinDraw[7]);

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

	//B fit
	CosHistB1->SetBinContent(1,xsecStatB[0]->GetMean());
	CosHistB1->SetBinError(1,xsecStatB[0]->GetRMS());
	CosHistB1->SetBinContent(2,xsecStatB[1]->GetMean());
	CosHistB1->SetBinError(2,xsecStatB[1]->GetRMS());

	CosHistB2->SetBinContent(1,xsecStatB[2]->GetMean());
	CosHistB2->SetBinError(1,xsecStatB[2]->GetRMS());
	CosHistB2->SetBinContent(2,xsecStatB[3]->GetMean());
	CosHistB2->SetBinError(2,xsecStatB[3]->GetRMS());
	CosHistB2->SetBinContent(3,xsecStatB[4]->GetMean());
	CosHistB2->SetBinError(3,xsecStatB[4]->GetRMS());

	CosHistB3->SetBinContent(1,xsecStatB[5]->GetMean());
	CosHistB3->SetBinError(1,xsecStatB[5]->GetRMS());
	CosHistB3->SetBinContent(2,xsecStatB[6]->GetMean());
	CosHistB3->SetBinError(2,xsecStatB[6]->GetRMS());
	CosHistB3->SetBinContent(3,xsecStatB[7]->GetMean());
	CosHistB3->SetBinError(3,xsecStatB[7]->GetRMS());

	CosHistB4->SetBinContent(1,xsecStatB[8]->GetMean());
	CosHistB4->SetBinError(1,xsecStatB[8]->GetRMS());
	CosHistB4->SetBinContent(2,xsecStatB[9]->GetMean());
	CosHistB4->SetBinError(2,xsecStatB[9]->GetRMS());
	CosHistB4->SetBinContent(3,xsecStatB[10]->GetMean());
	CosHistB4->SetBinError(3,xsecStatB[10]->GetRMS());

	CosHistB5->SetBinContent(1,xsecStatB[11]->GetMean());
	CosHistB5->SetBinError(1,xsecStatB[11]->GetRMS());
	CosHistB5->SetBinContent(2,xsecStatB[12]->GetMean());
	CosHistB5->SetBinError(2,xsecStatB[12]->GetRMS());
	CosHistB5->SetBinContent(3,xsecStatB[13]->GetMean());
	CosHistB5->SetBinError(3,xsecStatB[13]->GetRMS());

	CosHistB6->SetBinContent(1,xsecStatB[14]->GetMean());
	CosHistB6->SetBinError(1,xsecStatB[14]->GetRMS());
	CosHistB6->SetBinContent(2,xsecStatB[15]->GetMean());
	CosHistB6->SetBinError(2,xsecStatB[15]->GetRMS());
	CosHistB6->SetBinContent(3,xsecStatB[16]->GetMean());
	CosHistB6->SetBinError(3,xsecStatB[16]->GetRMS());

	CosHistB7->SetBinContent(1,xsecStatB[17]->GetMean());
	CosHistB7->SetBinError(1,xsecStatB[17]->GetRMS());
	CosHistB7->SetBinContent(2,xsecStatB[18]->GetMean());
	CosHistB7->SetBinError(2,xsecStatB[18]->GetRMS());

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

	CosHistB1->SetMarkerStyle(8);
	CosHistB1->SetMarkerSize(1);
	CosHistB1->SetMarkerColor(kGreen+2);

	CosHistNEUT1->SetLineColor(kRed+2);
	CosHistNEUT1->SetLineStyle(7);


	leg = new TLegend(0.1,0.1,0.9,0.9);
	leg->AddEntry(CosHist1,"New Result","pel");
	leg->AddEntry(CosHistB1,"Old Result","pel");

	CosHist1->Draw("PE0");
	CosHistB1->Draw("same PE0");
	c->Print("~/Desktop/plots/CosByMom1.pdf");

	c = new TCanvas;
	leg->Draw();
	c->Print("~/Desktop/plots/CosByMomLegend.pdf");

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

	CosHistB2->SetMarkerStyle(8);
	CosHistB2->SetMarkerSize(1);
	CosHistB2->SetMarkerColor(kGreen+2);

	CosHistNEUT2->SetLineColor(kRed+2);
	CosHistNEUT2->SetLineStyle(7);

	CosHist2->Draw("PE0");
	CosHistB2->Draw("same PE0");
	c->Print("~/Desktop/plots/CosByMom2.pdf");

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

	CosHistB3->SetMarkerStyle(8);
	CosHistB3->SetMarkerSize(1);
	CosHistB3->SetMarkerColor(kGreen+2);

	CosHistNEUT3->SetLineColor(kRed+2);
	CosHistNEUT3->SetLineStyle(7);


	CosHist3->Draw("PE0");
	CosHistB3->Draw("same PE0");
	c->Print("~/Desktop/plots/CosByMom3.pdf");

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

	CosHistB4->SetMarkerStyle(8);
	CosHistB4->SetMarkerSize(1);
	CosHistB4->SetMarkerColor(kGreen+2);

	CosHistNEUT4->SetLineColor(kRed+2);
	CosHistNEUT4->SetLineStyle(7);


	CosHist4->Draw("PE0");
	CosHistB4->Draw("same PE0");
	c->Print("~/Desktop/plots/CosByMom4.pdf");

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

	CosHistB5->SetMarkerStyle(8);
	CosHistB5->SetMarkerSize(1);
	CosHistB5->SetMarkerColor(kGreen+2);

	CosHistNEUT5->SetLineColor(kRed+2);
	CosHistNEUT5->SetLineStyle(7);


	CosHist5->Draw("PE0");
	CosHistB5->Draw("same PE0");
	c->Print("~/Desktop/plots/CosByMom5.pdf");

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

	CosHistB6->SetMarkerStyle(8);
	CosHistB6->SetMarkerSize(1);
	CosHistB6->SetMarkerColor(kGreen+2);

	CosHistNEUT6->SetLineColor(kRed+2);
	CosHistNEUT6->SetLineStyle(7);

	CosHist6->Draw("PE0");
	CosHistB6->Draw("same PE0");
	c->Print("~/Desktop/plots/CosByMom6.pdf");

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

	CosHistB7->SetMarkerStyle(8);
	CosHistB7->SetMarkerSize(1);
	CosHistB7->SetMarkerColor(kGreen+2);

	CosHistNEUT7->SetLineColor(kRed+2);
	CosHistNEUT7->SetLineStyle(7);

	CosHist7->Draw("PE0");
	CosHistB7->Draw("same PE0");
	c->Print("~/Desktop/plots/CosByMom7.pdf");

	// single diff
	Float_t PBinsDraw[8]={400,530,670,800,1000,1380,2010,3410};
	Float_t* PBinWidths = new Float_t[7];
	for(int ii=0; ii<7; ii++){
		PBinWidths[ii]=(PBinsDraw[ii+1]-PBinsDraw[ii]);
	}
	TH1F* MomFitResult = new TH1F("MomFitResult","",7,PBinsDraw);
	TH1F* MomFitResultNEUT = new TH1F("MomFitResultNEUT","",7,PBinsDraw);
	for(int ii=0; ii<7; ii++){
		MomFitResult->SetBinContent(ii+1, xsecMomStat[ii]->GetMean()/PBinWidths[ii]);
		MomFitResult->SetBinError(ii+1, xsecMomStat[ii]->GetRMS()/PBinWidths[ii]);
		MomFitResultNEUT->SetBinContent(ii+1, xsecMomStatNEUT[ii]->GetMean()/PBinWidths[ii]);
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

	leg = new TLegend(0.7,0.7,0.9,0.9);
	leg->AddEntry(MomFitResult,"Post Fit MC","pel");
	leg->AddEntry(MomFitResultNEUT,"NEUT RFG+RPA","l");

	MomFitResult->Draw("PE0");
	leg->Draw("same");
	c->Print("~/Desktop/plots/MomOverlay.pdf");
}
