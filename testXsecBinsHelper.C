#include <iostream>
#include "./xsecBinsHelper.C"

//here are tests and development scratch work for the xsecBinsHelper class
//
//Oct. 18, 2018
//Author: Thomas Campbell <thomascampbell1@gmail.com>

bool wasBinGood(Float_t, Float_t, int, xsecBinsHelper*);
bool wasBinGoodIPS(Float_t, Float_t, int, xsecBinsHelper*);

void testXsecBinsHelper(){
	cout << "calculating xsec" << endl;
	xsecBinsHelper* helper = new xsecBinsHelper();
	helper->PrintBins();

	//some test data (validated by hand ... :< ) 
	if(!wasBinGood(555,0.97,7,helper)) return;
	if(!wasBinGood(1550,0.945,21,helper)) return;
	if(!wasBinGood(5000,0.945,27,helper)) return;
	if(!wasBinGood(50,0.945,0,helper)) return;
	
	if(!wasBinGoodIPS(50,0.95,-1,helper)) return;
	if(!wasBinGoodIPS(900,0.99,10,helper))	return;
	if(!wasBinGoodIPS(350,0.4,-1,helper))	return;
	if(!wasBinGoodIPS(690,0.85,-1,helper))	return;
	if(!wasBinGoodIPS(1400,0.5,-1,helper))	return;
	if(!wasBinGoodIPS(420,0.85,0,helper))	return;
	if(!wasBinGoodIPS(650,0.97,4,helper))	return;
	if(!wasBinGoodIPS(3300,0.99,18,helper))	return;

	cout << "all tests pass, nice work!" << endl;
	return;
}

bool wasBinGood(Float_t p, Float_t cos, int bin, xsecBinsHelper* helper) {
	if(helper->GetBin(p,cos)!=bin){
		cout << "wasBinGood Fail! bailing out..." << p << " " << cos << endl;
		return false;
	}
	return true;
}
bool wasBinGoodIPS(Float_t p, Float_t cos, int bin, xsecBinsHelper* helper) {
	if(helper->GetBinIPS(p,cos)!=bin){
		cout << "wasBinGoodIPS Fail! bailing out..." << p << " " << cos << endl;
		return false;
	}
	return true;
}
