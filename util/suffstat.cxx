#define SUFFSTAT_C
#include <TMath.h>

suffStat::suffStat(){
	sx=0.;
	sxx=0.;
	sxxx=0.;
	n=0.;
}
suffStat::~suffStat(){}

void suffStat::Fill(Float_t in){
	sx+=in;
	sxx+=(in*in);
	sxxx+=(in*in*in);
	n+=1.;
}

Float_t suffStat::GetMean(){
	return (sx/n);
}

Float_t suffStat::GetRMS(){
	return (TMath::Sqrt(sxx-sx*sx)/n);
}
