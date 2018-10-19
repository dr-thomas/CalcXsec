#define COVMAT_C
#include "covmat.hxx"

#include <iostream>
#include <TMatrix.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>

covMatD::covMatD(Int_t dimension){
	dim = dimension;
	vec = new Double_t[dim];
	mat = new Double_t*[dim];
	matDecomp = new Double_t*[dim];
	for (int ii=0; ii<dim; ii++) {
		mat[ii] = new Double_t[dim];
		matDecomp[ii] = new Double_t[dim];
	}
}
covMatD::~covMatD(){
	for (int ii=0; ii<dim; ii++) {
		delete mat[ii];
		delete matDecomp[ii];
	}
	delete mat;
	delete matDecomp;
	delete vec;
}

// Setting stuff
void covMatD::SetMat(Int_t ii, Int_t jj, Double_t val) {
	if (ii>=dim || jj>=dim) {
		return;
	}
	mat[ii][jj]=val;
}
void covMatD::SetVec(Int_t ii, Double_t val) {
	if (ii>=dim) { 
		return;
	}
	vec[ii]=val;
}

// Retrieving stuff
TMatrixDSym* covMatD::GetTMatrix() {
	TMatrixDSym* out = new TMatrixDSym(dim);
	for (int ii=0; ii<dim; ii++) {
		for (int jj=0; jj<dim; jj++) {
			(*out)(ii,jj) = mat[ii][jj];
		}
	}
	return out;
}

// Throwing stuff
void covMatD::SetSeed(int in){
	seed = in;
}
void covMatD::Decompose() {
	TMatrixDSym dmat(dim);
	for (int ii=0; ii<dim; ii++) {
		for (int jj=0; jj<dim; jj++) {
			dmat(ii,jj) = mat[ii][jj];
		}
	}
	TDecompChol chol(dmat);

	chol.Decompose();
	TMatrix cholU = chol.GetU();
	cholU.T();
	for (int ii=0; ii<dim; ii++) {
		for (int jj=0; jj<dim; jj++) {
			matDecomp[ii][jj]=cholU(ii,jj);
		}
	}
}

/* Throwing stuff
 * what to return: set vector as central values, return vecotor of 
 * averages and uncertainties, or maybe a result struct <- that one
 *
 * want ability to add multiple covariance matricies and specify 
 * similar paramters (eg correlated throwing of flux parameters from 
 * fit cov and BANFF matrix)?  -> maybe not nesseary?? 
 */

