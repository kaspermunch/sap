#ifndef PARM_BASEFREQS_H
#define PARM_BASEFREQS_H

#include "parm.h"
#include "MbVector.h"

using namespace std;



class BaseFreqs : public Parm {

	public:
                	              BaseFreqs(MbRandom *rn, Model *mp, int ns, double a0);  
								  BaseFreqs(BaseFreqs &b);
					  BaseFreqs   &operator=(BaseFreqs &b);
						 double   change(void);
						   void   print(void);
						   void   clone(BaseFreqs &b);
					     double   getLnPriorProbability(void);
						 string   getParmName(void) { return parmName; }
			   MbVector<double>   getVal(void) { return freqs; }
			             double   getVal(int i) { return freqs[i]; }
                	              
	private:
                            int   numStates;
			   MbVector<double>   freqs;
			   MbVector<double>   a;
			   MbVector<double>   dirParm;
						 double   alpha0;
                            
};

#endif
