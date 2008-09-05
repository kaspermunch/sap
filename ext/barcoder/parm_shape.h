#ifndef PARM_SHAPE_H
#define PARM_SHAPE_H

#include "parm.h"
#include "MbVector.h"

using namespace std;

class GammaShape : public Parm {

	public:
                	              GammaShape(MbRandom *rn, Model *mp, double tn, double lam, int nc, double x);  
								  GammaShape(GammaShape &b);
					 GammaShape   &operator=(GammaShape &b);
						 double   change(void);
						   void   print(void);
						   void   clone(GammaShape &p);
					     double   getLnPriorProbability(void);
						 string   getParmName(void) { return parmName; }
						    int   getNumCats(void) { return numCats; }
			             double   getVal(void) { return alpha; }
			   MbVector<double>   getRates(void) { return r; }
                	              
	private:
					     double   alpha;
						 double   lambda;
						 double   tuning;
						    int   numCats;
			   MbVector<double>   r;
                            
};

#endif
