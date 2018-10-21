#ifndef SUFFSTAT_H
#define SUFFSTAT_H

class suffStat {
	public:
		suffStat(Float_t scale=1.);
		~suffStat();

		void Fill(Float_t);
		Float_t GetMean();
		Float_t GetRMS();
		void Dump();
	private:
		Float_t sx;
		Float_t sxx;
		Float_t sxxx;
		Float_t n;
		Float_t scale;
};
#endif
