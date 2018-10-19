#ifndef XSECBINSHELPER_H
#define XSECBINSHELPER_H

// xsecBinsHelper has utilities for and the definition of the bins
//
//Oct. 18, 2018
//Author: Thomas Campbell <thomascampbell1@gmail.com>
//

class xsecBinsHelper {
	public:
		xsecBinsHelper();
		~xsecBinsHelper();

		void PrintBins();
		int GetPBin(Float_t p);
		int GetCosBin(Float_t cos, int PBin);
		int GetBin(Float_t p, Float_t cos);
		int GetBinIPS(Float_t p, Float_t cos);

	private:
		Float_t* PBins;
		Float_t** CosBins;
		Int_t* BinMap;
};
#endif
