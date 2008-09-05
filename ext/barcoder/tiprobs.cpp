#include <iostream>
#include "parm_shape.h"
#include "parm_tree.h"
#include "tiprobs.h"
#include "MbMatrix.h"
#include "MbTransitionMatrix.h"

using namespace std;



TiProbs::TiProbs(MbTransitionMatrix *tm, int nt, int ngc) {

	tmatrixPtr   = tm;
	numTaxa      = nt;
	numNodes     = 2 * numTaxa;
	numGammaCats = ngc;
	
	tis[0] = new MbMatrix<double>*[numGammaCats];
	tis[1] = new MbMatrix<double>*[numGammaCats];
	tis[0][0] = new MbMatrix<double>[numGammaCats * numNodes];
	tis[1][0] = new MbMatrix<double>[numGammaCats * numNodes];
	for (int i=1; i<numGammaCats; i++)
		{
		tis[0][i] = tis[0][i-1] + numNodes;
		tis[1][i] = tis[1][i-1] + numNodes;
		}
	for (int k=0; k<numGammaCats; k++)
		{
		for (int i=0; i<numNodes; i++)
			{
			tis[0][k][i] = MbMatrix<double>(4, 4);
			tis[1][k][i] = MbMatrix<double>(4, 4);
			}
		}
		
}



TiProbs::~TiProbs(void) {

	delete [] tis[0][0];
	delete [] tis[1][0];
	delete [] tis[0];
	delete [] tis[1];
	
}



MbMatrix<double> &TiProbs::getTiMatrix(int whichSpace, int whichNode, int whichGammaCat) {

	return tis[ whichSpace ][ whichGammaCat ][ whichNode ];
	
}



void TiProbs::print(void) {

	for (int n=0; n<numNodes; n++)
		{
		cout << "Node " << n << ":" << endl;
		for (int k=0; k<numGammaCats; k++)
			{
			for (int i=0; i<4; i++)
				{
				cout << "   ";
				for (int j=0; j<4; j++)
					cout << fixed << setprecision(8) << tis[0][k][n][i][j] << " ";
				cout << "   ";
				for (int j=0; j<4; j++)
					cout << fixed << setprecision(8) << tis[1][k][n][i][j] << " ";
				cout << endl;
				}
			}
		}
		
}



void TiProbs::upDateTiProbs(Tree *t, GammaShape *r) {

	MbVector<double> rates = r->getRates();
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node *p = t->getDownPassNode(n);
		int idx = p->getIndex();

		if ( p->getUpdateTi() == true && p->getAnc() != NULL )
			{
			double v = p->getV();
			for (int k=0; k<numGammaCats; k++)
				{
				double vr = v * rates[k];
				MbMatrix<double> A = getTiMatrix( p->getActiveTi(), idx, k );
				A = tmatrixPtr->tiProbs(vr, A);
				}
			p->setUpdateTi( false );
			}
		}

}




