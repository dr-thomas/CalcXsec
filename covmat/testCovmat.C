#include "./covmat.hxx"
#include "./covmat.cxx"

#include <TMatrixDSym.h>
#include <TMatrix.h>
#include <TFile.h>
#include <vector>
#include <TAxis.h>

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
	TMatrixDSym* newmat = mat->GetTMatrix();
	newmat->Print();
	mat->Decompose();
	newmat = mat->GetTMatrix();
	newmat->Print();
}
