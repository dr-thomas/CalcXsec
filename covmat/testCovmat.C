#include "./covmat.hxx"
#include "./covmat.cxx"

#include <TMatrixDSym.h>
#include <TMatrix.h>
#include <TFile.h>
#include <vector>
#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>

void testCovmat(){
	TFile* finFluxCov = new TFile("./flux_covariance_banff_13av1.1.root", "OPEN");
	TAxis* numubBins = (TAxis*)finFluxCov->Get("nd5_anumode_numub_bins");

	TMatrixDSym* covIn   = (TMatrixDSym*)finFluxCov->Get("total_flux_cov");
	int nBins = numubBins->GetNbins();
	TMatrixDSym cov(nBins);

	for(int i=0;i<nBins;i++){
		for(int j=0;j<nBins;j++){
			cov(i, j) = (*covIn)(i,j);
      }
    }
    finFluxCov->Close();

	covMatD* mat = new covMatD(nBins);
	for (int ii=0; ii<nBins; ii++) {
		for (int jj=0; jj<nBins; jj++) {
			mat->SetMat(ii, jj, cov(ii,jj));
		}
	}

	TCanvas* c = new TCanvas;

	TMatrixDSym* newmat = mat->GetTMatrix();
	mat->Decompose();
	TMatrixD* othermat = mat->GetDecompTMatrix();

	TMatrixD testTrsp = (*othermat);
	testTrsp.T();
	TMatrixD testMat = (*othermat)*testTrsp;

	bool wasSame=true;
	for (int ii=0; ii<nBins; ii++) {
		for (int jj=0; jj<nBins; jj++) {
			if(TMath::Abs((*newmat)(ii,jj)-(testMat)(ii,jj))>1e-6){
				wasSame=false;
			}
		}
	}


	if(!wasSame && nBins>0){
		cout << "matrix decmposion failed!  dumping plots and bailing out..." << endl;

		newmat->Draw("colz");
		c->Print("./plots/cov.pdf");
		c = new TCanvas;
		testMat.Draw("colz");
		c->Print("./plots/shouldBeCov.pdf");

		return;
	}

	int nThrows = 1000;
	TH1D* someResult = new TH1D("someResult","",1000,-5,5);
	for(int ii=0; ii<nThrows; ii++) {
		mat->Throw();
		someResult->Fill(mat->varVec[nBins-1]);
	}
	if(TMath::Abs((someResult->GetRMS())-TMath::Sqrt((*newmat)(nBins-1,nBins-1)))>0.03){
		cout << "variation is different than expected, bailing out..." << endl;
	}

	cout << "all tests passed, ahhhh yeah!" << endl;

}

