#ifndef PARM_H
#define PARM_H

#include <string>

using namespace std;



class BaseFreqs;
class GammaShape;
class MbRandom;
class Model;
class SubRates;
class Tree;
class Parm {

	public:
                                  Parm(MbRandom *rn, Model *mp);
				           Parm   &operator=(Parm &p);
				 virtual double   change(void)=0;
		         virtual double   getLnPriorProbability(void)=0;
				 virtual string   getParmName(void)=0;
			       virtual void   print(void)=0;

	private:
				         double   lnPriorProbability;

	protected:
				       MbRandom   *ranPtr;
					      Model   *modelPtr;
						 string   parmName;

};

#endif
