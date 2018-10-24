#include <iostream>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void DrawXsecSFFitComp(Float_t** nData, Float_t** nDataSF, Float_t** nSel, Float_t** nGen, Float_t* binWidth, Float_t* intFlux, Float_t* nTargets, Float_t nTargetsNomMC, int nToys){

	//calc xsec
	suffStat** xsecStat = new suffStat*[19];
	suffStat** xsecStatSF = new suffStat*[19];
	suffStat** xsecStatNEUT = new suffStat*[19];
	for(int ii=0; ii<19; ii++){
		xsecStat[ii] = new suffStat(1e-39);
		xsecStatSF[ii] = new suffStat(1e-39);
		xsecStatNEUT[ii] = new suffStat(1e-39);
	}
	suffStat** xsecMomStat = new suffStat*[7];
	suffStat** xsecMomStatSF = new suffStat*[7];
	suffStat** xsecMomStatNEUT = new suffStat*[7];
	for(int ii=0; ii<7; ii++){
		xsecMomStat[ii] = new suffStat(1e-39);
		xsecMomStatSF[ii] = new suffStat(1e-39);
		xsecMomStatNEUT[ii] = new suffStat(1e-39);
	}
	Float_t* singleBin = new Float_t[400];
	Float_t* singleBinSF = new Float_t[400];
	Float_t* singleBinN = new Float_t[400];
	for(int ii=0; ii<nToys; ii++){
		singleBin[ii]=0.;
		singleBinSF[ii]=0.;
		singleBinN[ii]=0.;
	}

	//TODO: single diff shit is sloppy af, should re-write

	for(int iToy=0; iToy<nToys; iToy++){
		Float_t xsecMom=0.;
		Float_t xsecMomSF=0.;
		Float_t xsecMomN=0.;
		int nBinsMom=2;
		int nDrawnP=0;
		int pIndex=0;
		for(int ii=0; ii<19; ii++){
			Float_t xsec = nData[iToy][ii]/(nSel[iToy][ii]/nGen[iToy][ii]);
			Float_t xsecSF = nDataSF[iToy][ii]/(nSel[iToy][ii]/nGen[iToy][ii]);
			xsec*=1.0/binWidth[ii];
			xsec*=1.0/intFlux[iToy];
			xsec*=1.0/nTargets[iToy];
			xsecSF*=1.0/binWidth[ii];
			xsecSF*=1.0/intFlux[iToy];
			xsecSF*=1.0/nTargets[iToy];
			xsecStat[ii]->Fill(xsec);
			xsecStatSF[ii]->Fill(xsecSF);
			singleBin[iToy]+=xsec;
			singleBinSF[iToy]+=xsecSF;
			Float_t xsecN = nGen[iToy][ii]/(binWidth[ii]*intFlux[iToy]);
			xsecN*=1.0/nTargetsNomMC;
			xsecStatNEUT[ii]->Fill(xsecN);
			singleBinN[iToy]+=xsecN;
			//Momentum single (Note: not normalized by bin width)
			xsec*=binWidth[ii];
			xsecSF*=binWidth[ii];
			xsecN*=binWidth[ii];
			if(nDrawnP<nBinsMom){
				xsecMom+=xsec;
				xsecMomSF+=xsecSF;
				xsecMomN+=xsecN;
				nDrawnP++;
				if(ii==18) {
					xsecMomStat[pIndex]->Fill(xsecMom);
					xsecMomStatSF[pIndex]->Fill(xsecMomSF);
					xsecMomStatNEUT[pIndex]->Fill(xsecMomN);
				}
			} else {
				nDrawnP=0;
				xsecMomStat[pIndex]->Fill(xsecMom);
				xsecMomStatSF[pIndex]->Fill(xsecMomSF);
				xsecMomStatNEUT[pIndex]->Fill(xsecMomN);
				xsecMom=0.;
				xsecMomSF=0.;
				xsecMomN=0.;
				pIndex++;
				if(ii==16) nBinsMom=2;
				else nBinsMom=3;
				//add again
				xsecMom+=xsec;
				xsecMomSF+=xsecSF;
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
	suffStat* sbStatSF = new suffStat(1e-39);
	suffStat* sbStatN = new suffStat(1e-39);
	for(int ii=0; ii<nToys; ii++){
		sbStat->Fill(singleBin[ii]);
		sbStatSF->Fill(singleBinSF[ii]);
		sbStatN->Fill(singleBinN[ii]);
	}
	cout << "Single Bin Results:" << endl;
	cout << "Data: " << sbStat->GetMean() << " +- " << sbStat->GetRMS() << endl;
	cout << "Data SF Fit: " << sbStatSF->GetMean() << " +- " << sbStatSF->GetRMS() << endl;
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

	TH1F* CosHistSF1 = new TH1F("CosHistSF1","",2,CosBinDraw[1]);
	TH1F* CosHistSF2 = new TH1F("CosHistSF2","",3,CosBinDraw[2]);
	TH1F* CosHistSF3 = new TH1F("CosHistSF3","",3,CosBinDraw[3]);
	TH1F* CosHistSF4 = new TH1F("CosHistSF4","",3,CosBinDraw[4]);
	TH1F* CosHistSF5 = new TH1F("CosHistSF5","",3,CosBinDraw[5]);
	TH1F* CosHistSF6 = new TH1F("CosHistSF6","",3,CosBinDraw[6]);
	TH1F* CosHistSF7 = new TH1F("CosHistSF7","",2,CosBinDraw[7]);

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

	//SF fit
	CosHistSF1->SetBinContent(1,xsecStatSF[0]->GetMean());
	CosHistSF1->SetBinError(1,xsecStatSF[0]->GetRMS());
	CosHistSF1->SetBinContent(2,xsecStatSF[1]->GetMean());
	CosHistSF1->SetBinError(2,xsecStatSF[1]->GetRMS());

	CosHistSF2->SetBinContent(1,xsecStatSF[2]->GetMean());
	CosHistSF2->SetBinError(1,xsecStatSF[2]->GetRMS());
	CosHistSF2->SetBinContent(2,xsecStatSF[3]->GetMean());
	CosHistSF2->SetBinError(2,xsecStatSF[3]->GetRMS());
	CosHistSF2->SetBinContent(3,xsecStatSF[4]->GetMean());
	CosHistSF2->SetBinError(3,xsecStatSF[4]->GetRMS());

	CosHistSF3->SetBinContent(1,xsecStatSF[5]->GetMean());
	CosHistSF3->SetBinError(1,xsecStatSF[5]->GetRMS());
	CosHistSF3->SetBinContent(2,xsecStatSF[6]->GetMean());
	CosHistSF3->SetBinError(2,xsecStatSF[6]->GetRMS());
	CosHistSF3->SetBinContent(3,xsecStatSF[7]->GetMean());
	CosHistSF3->SetBinError(3,xsecStatSF[7]->GetRMS());

	CosHistSF4->SetBinContent(1,xsecStatSF[8]->GetMean());
	CosHistSF4->SetBinError(1,xsecStatSF[8]->GetRMS());
	CosHistSF4->SetBinContent(2,xsecStatSF[9]->GetMean());
	CosHistSF4->SetBinError(2,xsecStatSF[9]->GetRMS());
	CosHistSF4->SetBinContent(3,xsecStatSF[10]->GetMean());
	CosHistSF4->SetBinError(3,xsecStatSF[10]->GetRMS());

	CosHistSF5->SetBinContent(1,xsecStatSF[11]->GetMean());
	CosHistSF5->SetBinError(1,xsecStatSF[11]->GetRMS());
	CosHistSF5->SetBinContent(2,xsecStatSF[12]->GetMean());
	CosHistSF5->SetBinError(2,xsecStatSF[12]->GetRMS());
	CosHistSF5->SetBinContent(3,xsecStatSF[13]->GetMean());
	CosHistSF5->SetBinError(3,xsecStatSF[13]->GetRMS());

	CosHistSF6->SetBinContent(1,xsecStatSF[14]->GetMean());
	CosHistSF6->SetBinError(1,xsecStatSF[14]->GetRMS());
	CosHistSF6->SetBinContent(2,xsecStatSF[15]->GetMean());
	CosHistSF6->SetBinError(2,xsecStatSF[15]->GetRMS());
	CosHistSF6->SetBinContent(3,xsecStatSF[16]->GetMean());
	CosHistSF6->SetBinError(3,xsecStatSF[16]->GetRMS());

	CosHistSF7->SetBinContent(1,xsecStatSF[17]->GetMean());
	CosHistSF7->SetBinError(1,xsecStatSF[17]->GetRMS());
	CosHistSF7->SetBinContent(2,xsecStatSF[18]->GetMean());
	CosHistSF7->SetBinError(2,xsecStatSF[18]->GetRMS());

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

	CosHistSF1->SetMarkerStyle(8);
	CosHistSF1->SetMarkerSize(1);
	CosHistSF1->SetMarkerColor(kGreen+2);

	CosHistNEUT1->SetLineColor(kRed+2);
	CosHistNEUT1->SetLineStyle(7);


	leg = new TLegend(0.1,0.1,0.9,0.9);
	leg->AddEntry(CosHist1,"Nominal Fit","pel");
	leg->AddEntry(CosHistSF1,"SF Tuned Fit","pel");
	leg->AddEntry(CosHistNEUT1,"Pre Fit MC (NEUT)","l");

	CosHist1->Draw("PE0");
	CosHistSF1->Draw("same PE0");
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

	CosHistSF2->SetMarkerStyle(8);
	CosHistSF2->SetMarkerSize(1);
	CosHistSF2->SetMarkerColor(kGreen+2);

	CosHistNEUT2->SetLineColor(kRed+2);
	CosHistNEUT2->SetLineStyle(7);

	CosHist2->Draw("PE0");
	CosHistSF2->Draw("same PE0");
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

	CosHistSF3->SetMarkerStyle(8);
	CosHistSF3->SetMarkerSize(1);
	CosHistSF3->SetMarkerColor(kGreen+2);

	CosHistNEUT3->SetLineColor(kRed+2);
	CosHistNEUT3->SetLineStyle(7);


	CosHist3->Draw("PE0");
	CosHistSF3->Draw("same PE0");
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

	CosHistSF4->SetMarkerStyle(8);
	CosHistSF4->SetMarkerSize(1);
	CosHistSF4->SetMarkerColor(kGreen+2);

	CosHistNEUT4->SetLineColor(kRed+2);
	CosHistNEUT4->SetLineStyle(7);


	CosHist4->Draw("PE0");
	CosHistSF4->Draw("PE0");
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

	CosHistSF5->SetMarkerStyle(8);
	CosHistSF5->SetMarkerSize(1);
	CosHistSF5->SetMarkerColor(kGreen+2);

	CosHistNEUT5->SetLineColor(kRed+2);
	CosHistNEUT5->SetLineStyle(7);


	CosHist5->Draw("PE0");
	CosHistSF5->Draw("same PE0");
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

	CosHistSF6->SetMarkerStyle(8);
	CosHistSF6->SetMarkerSize(1);
	CosHistSF6->SetMarkerColor(kGreen+2);

	CosHistNEUT6->SetLineColor(kRed+2);
	CosHistNEUT6->SetLineStyle(7);

	CosHist6->Draw("PE0");
	CosHistSF6->Draw("same PE0");
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

	CosHistSF7->SetMarkerStyle(8);
	CosHistSF7->SetMarkerSize(1);
	CosHistSF7->SetMarkerColor(kGreen+2);

	CosHistNEUT7->SetLineColor(kRed+2);
	CosHistNEUT7->SetLineStyle(7);

	CosHist7->Draw("PE0");
	CosHistSF7->Draw("same PE0");
	CosHistNEUT7->Draw("same");
	c->Print("./plots/CosByMom7.pdf");

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
	leg->AddEntry(MomFitResultNEUT,"Pre Fit MC (NEUT)","l");

	MomFitResult->Draw("PE0");
	MomFitResultNEUT->Draw("same");
	leg->Draw("same");
	c->Print("./plots/MomOverlay.pdf");
}
