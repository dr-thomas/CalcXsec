#define SUFFSTAT_C
#include <TMath.h>
#include <iostream>

suffStat::suffStat(Float_t scaleIn){
	sx=0.;
	sxx=0.;
	sxxx=0.;
	n=0.;
	scale=scaleIn;
}
suffStat::~suffStat(){}

void suffStat::Fill(Float_t in){
	in*=1.0/scale;
	sx+=in;
	sxx+=(in*in);
	sxxx+=(in*in*in);
	n+=1.;
}

Float_t suffStat::GetMean(){
	return (sx/n)*scale;
}

Float_t suffStat::GetRMS(){
	return (TMath::Sqrt(n*sxx-sx*sx)/n)*scale;
}

void suffStat::Dump(){
	cout << "Dumping SuffStat Innards:" << endl;
	cout << "Sx:" << sx << endl;
	cout << "Sx2:" << sxx << endl;
	cout << "Sx3:" << sxxx << endl;
}
void suffStat::Reset(){
	sx=0.;
	sxx=0.;
	sxxx=0.;
	n=0.;
}
