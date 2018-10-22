#define COVMAT_C
#include "covmat.hxx"

#include <iostream>
#include <TMatrix.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TRandom3.h>

// covMatD has covariance matrix throwing utilities 
//
// Note: be very careful with varVec, directly accessing an array, no warning if
//       the user tries to access an elemet that does not exist! (TODO, should 
//       probably do something about this...)
//
//Oct. 18, 2018
//Author: Thomas Campbell <thomascampbell1@gmail.com>
//

covMatD::covMatD(Int_t dimension){
	dim = dimension;
	mat = new Double_t*[dim];
	varVec = new Double_t[dim];
	matDecomp = new Double_t*[dim];
	for (int ii=0; ii<dim; ii++) {
		mat[ii] = new Double_t[dim];
		matDecomp[ii] = new Double_t[dim];
	}
	randN.SetSeed(21344213);
	wasDecomposed=false;
}
covMatD::~covMatD(){
	for (int ii=0; ii<dim; ii++) {
		delete []mat[ii];
		delete []matDecomp[ii];
	}
	delete []mat;
	delete []matDecomp;
	delete []varVec;
}

// Setting stuff
void covMatD::SetMat(Int_t ii, Int_t jj, Double_t val) {
	if (ii>=dim || jj>=dim) {
		cerr << "he's a very bad man (add index larger than covMat dim)" << endl;
		return;
	}
	mat[ii][jj]=val;
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
TMatrixD* covMatD::GetDecompTMatrix() {
	TMatrixD* out = new TMatrixD(dim,dim);
	for (int ii=0; ii<dim; ii++) {
		for (int jj=0; jj<dim; jj++) {
			(*out)(ii,jj) = matDecomp[ii][jj];
		}
	}
	return out;
}

// Throwing stuff
void covMatD::SetSeed(int in){
	randN.SetSeed(in);
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
	wasDecomposed=true;
}

// Throw updates the variation vector which is to be directly accessed.
// Important Reminder: weight = 1. + variation 
void covMatD::Throw() {
	if(!wasDecomposed){
		cerr << "matrix was not decomposed!!!!!!!" << endl;
		cerr << "Throws will return 0 variations" << endl;
	}
	Double_t* randVec = new Double_t[dim];
	for(int ii=0; ii<dim; ii++){
		randVec[ii]=randN.Gaus();
		varVec[ii]=0.;
	}

	for(int ii=0; ii<dim; ii++){
		for(int jj=0; jj<dim; jj++){
			varVec[ii]+=matDecomp[ii][jj]*randVec[jj];
		}
	}
	delete []randVec;
}
