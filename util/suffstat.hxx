#ifndef SUFFSTAT_H
#define SUFFSTAT_H

class suffStat {
	public:
		suffStat();
		~suffStat();

		void Fill(Float_t);
		Float_t GetMean();
		Float_t GetRMS();
	private:
		Float_t sx;
		Float_t sxx;
		Float_t sxxx;
		Float_t n;
};
#endif
