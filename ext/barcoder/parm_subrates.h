#ifndef PARM_SUBRATES_H
#define PARM_SUBRATES_H

#include "parm.h"
#include "MbVector.h"

using namespace std;

class SubRates : public Parm {

	public:
                	              SubRates(MbRandom *rn, Model *mp, int ns, double a0);  
								  SubRates(SubRates &b);
				       SubRates   &operator=(SubRates &b);
						 double   change(void);
					     double   getLnPriorProbability(void);
						 string   getParmName(void) { return parmName; }
						   void   print(void);
						   void   clone(SubRates &b);
			   MbVector<double>   getVal(void) { return rates; }
			             double   getVal(int i) { return rates[i]; }
                	              
	private:
                            int   numStates;
							int   numRates;
			   MbVector<double>   rates;
			   MbVector<double>   a;
			   MbVector<double>   dirParm;
						 double   alpha0;
                            
};

#endif
