#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include "parm.h"



class Alignment;
class CondLike;
class Constraints;
class MbRandom;
class MbTransitionMatrix;
class Parm;
class Settings;
class TiProbs;
class Tree;
class Model {

	public:
                	              Model(MbRandom *rn, Alignment *al, Constraints *cs, Settings *sp);  
								  ~Model(void);
						   void   copyParms(int from, int to);
						   void   flipActiveParm(void) { (activeParm == 0 ? activeParm = 1 : activeParm = 0); }
						   void   flipAllActiveCls(void);
						   void   flipAllActiveTis(void);
						    int   getActiveParm(void) { return activeParm; }
					  Alignment   *getAlignmentPtr(void) { return alignmentPtr; }
					       Parm   *getCurParameter(int i);
					  BaseFreqs   *getCurBaseFreqs(void);
					   SubRates   *getCurSubRates(void);
					 GammaShape   *getCurGammaShape(void);
						   Tree   *getCurTree(void);
						    int   getNumParms(int i) { return parms[i].size(); }
					     double   getLnL(void) { return lnL; }
						 double   lnLikelihood(void);
						    int   pickParm(void);
						   void   restoreQ(Parm *p);
						   void   setAllUpdateCls(bool tf);
						   void   setAllUpdateTis(bool tf);
						   void   setLnL(double x) { lnL = x; }
						   void   upDateTiProbs(void);
						   void   upDateQ(Parm *p);

	private:
	                        int   activeParm;
					  Alignment   *alignmentPtr;
					Constraints   *constraintsPtr;
					   CondLike   *condLikes;
						 double   lnL;
						 double   oldLnL;
		         vector<Parm *>   parms[2];
				 vector<double>   parmProb;
					   MbRandom   *ranPtr;
			 MbTransitionMatrix   *tiMatrix;
						TiProbs   *tiProbs;

};

#endif
