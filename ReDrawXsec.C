#include <TFile.h>
#include <TMatrixF.h>
#include <TMatrixD.h>
#include <TVectorF.h>

#include "./draw/DrawXsec.C"

void ReDrawXsec(){

	//TODO: something is broken between here and return

	TFile* inF = new TFile("./outCalcXsec.root", "OPEN");
	TMatrixF* nDataM = (TMatrixF*)inF->Get("nDataM");
	TMatrixF* nSelM = (TMatrixF*)inF->Get("nSelM");
	TMatrixF* nGenM = (TMatrixF*)inF->Get("nGenM");
	TMatrixF* nGenSFM = (TMatrixF*)inF->Get("nGenSFM");

	TVectorF* binWidthV = (TVectorF*)inF->Get("binWidthV");
	TVectorF* intFluxV = (TVectorF*)inF->Get("intFluxV");
	TVectorF* nTargetsV = (TVectorF*)inF->Get("nTargetsV");

	int nToys = nDataM->GetNrows();
	if(nToys==0){
		return;
	}

	cout << "here" << endl;
	return;


	Float_t** nData = new Float_t*[nToys];
	Float_t** nSel = new Float_t*[nToys];
	Float_t** nGen = new Float_t*[nToys];
	Float_t** nGenSF = new Float_t*[nToys];

	for(int ii=0; ii<nToys; ii++){
		nData[ii] = new Float_t[19];
		nSel[ii] = new Float_t[19];
		nGen[ii] = new Float_t[19];
		nGenSF[ii] = new Float_t[19];
	}

	Float_t* binWidth = new Float_t[19];
	for(int jj=0; jj<19; jj++){
		binWidth[jj] = (*binWidthV)(jj);
	}

	Float_t* intFlux = new Float_t[nToys];
	Float_t* nTargets = new Float_t[nToys];

	for(int ii=0; ii<nToys; ii++){
		intFlux[ii] = (*intFluxV)(ii);
		nTargets[ii] = (*nTargetsV)(ii);
		for(int jj=0; jj<19; jj++){
			nData[ii][jj] = (*nDataM)(ii,jj);
			nSel[ii][jj] = (*nSelM)(ii,jj);
			nGen[ii][jj] = (*nGenM)(ii,jj);
			nGenSF[ii][jj] = (*nGenSFM)(ii,jj);
		}
	}

	Float_t nTargetsNomMC=6.46e+28; 
	DrawXsec(nData,nSel,nGen,nGenSF,binWidth,intFlux,nTargets,nTargetsNomMC,nToys);
}
