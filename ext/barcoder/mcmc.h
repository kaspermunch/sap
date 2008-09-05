#ifndef MCMC_H
#define MCMC_H

#include <fstream>

class MbRandom;
class Model;
class Settings;
class Mcmc {

	public:
                	              Mcmc(Model *m, MbRandom *r, Settings *sp, IoManager *outputPtr);  

	private:
						    int   flip(int x);
	                       void   sampleChainState(int n, ofstream &parmOut, ofstream &treeOut);
						 double   lnL;
						  Model   *modelPtr;
					   MbRandom   *rnd;
						    int   chainLength;
						    int   printFreq;
						    int   sampleFreq;
                            
};

#endif
