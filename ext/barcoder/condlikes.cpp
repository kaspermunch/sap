#include "condlikes.h"
#include "alignment.h"



CondLike::CondLike(Alignment *al, int ngc) {

	/* initialze some basic variables */
	numGammaCats = ngc;
	numChar      = al->getNumChar();
	numTaxa      = al->getNumTaxa();
	numStates    = 4;
	
	/* allocate the conditional likelihood vector */
	cls = new double[2 * 2*numTaxa * numChar * numGammaCats * numStates];
	for (int i=0; i<2*2*numTaxa*numChar*numGammaCats*numStates; i++)
		cls[i] = 0.0;
		
	/* and also allocate a vector containing the number of sites of each pattern */
	numSitesOfPat = new int[numChar];
	
	/* allocate and initialize a matrix containing pointers to the appropriate position of the conditional likelihood vector */
	clsPtr = new double***[2];
	clsPtr[0] = new double**[2 * 2*numTaxa];
	clsPtr[1] = clsPtr[0] + (2*numTaxa);
	for (int i=0; i<2; i++)
		for (int j=0; j<2*numTaxa; j++)
			clsPtr[i][j] = new double*[numChar];
	for (int i=0; i<2; i++)
		for (int j=0; j<2*numTaxa; j++)
			for (int k=0; k<numChar; k++)
				clsPtr[i][j][k] = &cls[ i*(2*numTaxa*numChar*numGammaCats*numStates) + j*(numChar*numGammaCats*numStates) + k*(numGammaCats*numStates) ];
				
	/* initialize the conditional likelihoods */
	for (int i=0; i<numTaxa; i++)
		{
		for (int c=0; c<numChar; c++)
			{
			numSitesOfPat[c] = al->getNumOfPattern(c);
			int theNuc = al->getNucleotide(i, c);
			int nucs[4];
			al->getPossibleNucs(theNuc, &nucs[0]);
			double *p0 = clsPtr[0][i][c];
			double *p1 = clsPtr[1][i][c];
			for (int k=0; k<numGammaCats; k++)
				{
				for (int s=0; s<numStates; s++)
					{
					if (nucs[s] == 1)
						{
						p0[s] = 1.0;
						p1[s] = 1.0;
						}
					}
				p0 += numStates;
				p1 += numStates;
				}
			}
		}
	
}



CondLike::~CondLike(void) {

	for (int i=0; i<2; i++)
		for (int j=0; j<2*numTaxa; j++)
			delete [] clsPtr[i][j];
	delete [] clsPtr[0];
	delete [] clsPtr;
	delete [] cls;
	delete [] numSitesOfPat;
	
}






