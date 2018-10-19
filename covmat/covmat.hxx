#ifndef COVMAT_H
#define COVMAT_H

#include <TMatrixDSym.h>

// covMatD has covariance matrix throwing utilities 
//
//Oct. 18, 2018
//Author: Thomas Campbell <thomascampbell1@gmail.com>
//

class covMatD {

	public:
		covMatD(Int_t dimension);
		~covMatD();
		void SetMat(Int_t ii, Int_t jj, Double_t val);
		void SetVec(Int_t ii, Double_t val);
		TMatrixDSym* GetTMatrix();
		void SetSeed(int in);
		void Decompose();

	private:
		Double_t** mat;
		Double_t** matDecomp;
		Double_t* vec;
		Double_t* resCen;
		Double_t* resErr;
		Int_t dim;
		Int_t seed;
};
#endif
