#ifndef TIPROBS_H
#define TIPROBS_H

#include <vector>
#include "MbMatrix.h"



class GammaShape;
class MbTransitionMatrix;
class Tree;
class TiProbs {

	public:
                	              TiProbs(MbTransitionMatrix *tm, int nt, int ngc);  
								  ~TiProbs(void);
	           MbMatrix<double>   &getTiMatrix(int whichSpace, int whichNode, int whichGammaCat);
			               void   print(void);
			               void   upDateTiProbs(Tree *t, GammaShape *r);
			

	private:
	         MbTransitionMatrix   *tmatrixPtr;
						    int   numTaxa;
							int   numNodes;
							int   numGammaCats;
	           MbMatrix<double>   **tis[2];
                            
};

#endif
