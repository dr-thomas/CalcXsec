#ifndef COVMAT_H
#define COVMAT_H

#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include <TRandom3.h>

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
		TMatrixDSym* GetTMatrix();
		TMatrixD* GetDecompTMatrix();
		void SetSeed(int in);
		void Decompose();
		void Throw();

		Double_t* varVec;

	private:
		Double_t** mat;
		Double_t** matDecomp;
		Int_t dim;
		TRandom3 randN;
};
#endif
