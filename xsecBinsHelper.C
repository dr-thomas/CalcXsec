#include <iostream>
#include "Rtypes.h"

// xsecBinsHelper has utilities for and the definition of the bins
//
//Oct. 18, 2018
//Author: Thomas Campbell <thomascampbell1@gmail.com>
//
//TODO:
// - bin width
// - make into proper header with prototype and function defs separate


class xsecBinsHelper {
	public:
		xsecBinsHelper(){
			Float_t initPBins[10]={0,400,530,670,800,1000,1380,2010,3410,50000};
			Float_t initCosBins[9][5]={
				{-1,1,1,1,1},//only one bin here
				{-1,0.84,0.94,1,1},//only three bins here
				{-1,0.85,0.92,0.96,1},//Col. Forbins...
				{-1,0.88,0.93,0.97,1},
				{-1,0.90,0.94,0.97,1},
				{-1,0.91,0.95,0.97,1},
				{-1,0.94,0.97,0.98,1},
				{-1,0.95,0.98,1,1},//only three bins here
				{-1,1,1,1,1}//only one bin here
			};
			PBins = new Float_t[10];
			for(int ii=0; ii<10; ii++){
				PBins[ii] = initPBins[ii];
			}
			CosBins = new Float_t*[9];
			for(int ii=0; ii<9; ii++){
				CosBins[ii] = new Float_t[5];
			}
			for(int ii=0; ii<9; ii++){
				for(int jj=0; jj<5; jj++){
					CosBins[ii][jj] = initCosBins[ii][jj];
				}
			}

			Int_t initBinMap[36] = {
				0,0,0,0,
				1,2,3,3,
				4,5,6,7,
				8,9,10,11,
				12,13,14,15,
				16,17,18,19,
				20,21,22,23,
				24,25,26,26,
				27,27,27,27
			};
			BinMap = new Int_t[36];
			for(int ii=0; ii<36; ii++){
				BinMap[ii] = initBinMap[ii];
			}
		}
		~xsecBinsHelper(){
			delete []PBins;
			for(int ii=0; ii<9; ii++){
				delete []CosBins[ii];
			}
			delete []BinMap;
		}

		void PrintBins() {
			cout << "Printing Bins: " << endl;
			for (int ii=0; ii<9; ii++) {
				cout << PBins[ii] << " < p < " << PBins[ii+1] << endl;
				for (int jj=0; jj<5; jj++) {
					cout << CosBins[ii][jj] << " ";
				}
				cout << endl;
			}
		}

		int GetPBin(Float_t p) {
			int ii=0;
			while(p>PBins[ii]){
				ii++;
			}
			ii--;
			if(ii<0 || ii>9){
				return -1;
			}
			return ii;
		}
		int GetCosBin(Float_t cos, int PBin) {
			int ii=0;
			while(cos>CosBins[PBin][ii]){
				ii++;
			}
			ii--;
			if(ii<0 || ii>5){
				return -1;
			}
			return ii;
		}

		int GetBin(Float_t p, Float_t cos) {
			int PBin=GetPBin(p);
			int CosBin=GetCosBin(cos, PBin);
			return BinMap[CosBin+PBin*4];
		}
		int GetBinIPS(Float_t p, Float_t cos) {
			int bin=GetBin(p,cos);
			int PBin = GetPBin(p);
			int skip[9] = {0,1,4,8,12,16,20,24,27};
			for(int ii=0; ii<9; ii++) {
				if(bin==skip[ii]){
					return -1;
				}
			}
			return (bin-PBin-1);
		}

	private:
		Float_t* PBins;
		Float_t** CosBins;
		Int_t* BinMap;
};
