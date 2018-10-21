#include <iostream>
#include <TRandom3.h>
#include <TMath.h>

#include "./suffstat.hxx"
#include "./suffstat.cxx"

void testSuffStat(){
	cout << "testing suffStat" << endl;

	suffStat* testStat = new suffStat();

	TRandom3 randN;
	randN.SetSeed(32478126);
	for(int ii=0; ii<100000; ii++){
		testStat->Fill(randN.Gaus());
	}
	if(TMath::Abs(testStat->GetMean()-0.)>1e-2){
		cout << "Bad Mean " << testStat->GetMean() << endl;
		return;
	}
	if(TMath::Abs(testStat->GetRMS()-1.)>3e-2){
		cout << "Bad RMS " << testStat->GetRMS() << endl;
		return;
	}

	cout << "all tests passed!" << endl;
}

