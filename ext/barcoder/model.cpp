#include "model.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "alignment.h"
#include "condlikes.h"
#include "MbRandom.h"
#include "MbTransitionMatrix.h"
#include "parm.h"
#include "parm_basefreqs.h"
#include "parm_shape.h"
#include "parm_subrates.h"
#include "parm_tree.h"
#include "settings.h"
#include "tiprobs.h"

using namespace std;



Model::Model(MbRandom *rn, Alignment *al, Constraints *cs, Settings *sp) {

	/* initialize the pointers */
	ranPtr         = rn;
	alignmentPtr   = al;
	constraintsPtr = cs;

	/* initialize the conditional likelihoods */
	condLikes = new CondLike( alignmentPtr, sp->getNumGammaCats() );
	
	/* set the current parameter space: all parameters share the same space */
	activeParm = 0;

	cout << "   Setting phylogenetic model" << endl;
	cout << "      Tree topology         = Discrete Uniform(1,# Topologies), subject to constraint(s)" << endl;
	cout << "      Branch lengths        = Independent Exponential(" << sp->getBrlenLambda() << ") distribution" << endl;
	cout << "      Substitution rates    = Dirichlet(1,1,1,1,1,1) distribution" << endl;
	cout << "      Base frequencies      = Dirichlet(1,1,1,1) distribution" << endl;
	cout << "      Gamma shape parameter = Expoenential(2) distribution" << endl << endl;
	
	/* add the tree parameter */ 
	parms[0].push_back( new Tree(ranPtr, this, alignmentPtr->getNumTaxa(), false, sp->getBrlenLambda(), alignmentPtr->getTaxonNames(), constraintsPtr) );
	parms[1].push_back( new Tree(ranPtr, this, alignmentPtr->getNumTaxa(), false, sp->getBrlenLambda(), alignmentPtr->getTaxonNames(), constraintsPtr) );
	parmProb.push_back( 0.7 );

	/* add the base frequency parameter */
	parms[0].push_back( new BaseFreqs(ranPtr, this, 4, 300.0) );
	parms[1].push_back( new BaseFreqs(ranPtr, this, 4, 300.0) );	
	parmProb.push_back( 0.1 );

	/* add the substitution rates parameter */
	parms[0].push_back( new SubRates(ranPtr, this, 4, 500.0) );
	parms[1].push_back( new SubRates(ranPtr, this, 4, 500.0) );
	parmProb.push_back( 0.1 );

	/* add the gamma shape parameter */
	parms[0].push_back( new GammaShape(ranPtr, this, log(2.0), 2.0, sp->getNumGammaCats(), 0.5) );
	parms[1].push_back( new GammaShape(ranPtr, this, log(2.0), 2.0, sp->getNumGammaCats(), 0.5) );
	parmProb.push_back( 0.1 );

	/* make certain the values for the two sets of parameters are equal */
	copyParms( 0, 1 );

	/* set up the transition probability matrix calculator */
	tiMatrix = new MbTransitionMatrix( getCurSubRates()->getVal(), getCurBaseFreqs()->getVal(), true );
	
	/* set up the transition probabilities */
	tiProbs = new TiProbs( tiMatrix, alignmentPtr->getNumTaxa(), sp->getNumGammaCats() );

	/* update all of the transition probabilities */
	setAllUpdateTis( true );
	tiProbs->upDateTiProbs( getCurTree(), getCurGammaShape() );

	/* calculate the log likelihood of the current state */
	setAllUpdateCls( true );	
	lnL = lnLikelihood();
	oldLnL = lnL;
	//cout << "lnL = " << lnL << endl;

	setAllUpdateTis( true );
	flipAllActiveTis();
	tiProbs->upDateTiProbs( getCurTree(), getCurGammaShape() );
	setAllUpdateCls( true );
	flipAllActiveCls();
	lnL = lnLikelihood();
	//cout << "lnL2 = " << lnL << endl;

	/* print values for the parameters */
#	if 0
	for (int i=0; i<2; i++)
		{
		for (vector<Parm *>::iterator p=parms[i].begin(); p != parms[i].end(); p++)
			(*p)->print();
		}
#	endif

}


Model::~Model(void) {

	for (int i=0; i<2; i++)
		for (vector<Parm *>::iterator p=parms[i].begin(); p != parms[i].end(); p++)
			delete *p;
	delete condLikes;
	delete tiMatrix;
	delete tiProbs;
	

}



void Model::copyParms(int from, int to) {

	for (int i=0; i<parms[0].size(); i++)
		{
		Parm *pFrom = parms[from][i];
		Parm *pTo   = parms[ to ][i];
		(*pTo) = (*pFrom);
		}

}



void Model::flipAllActiveCls(void) {

	Tree *t = getCurTree();
	t->flipAllActiveCls();

}



void Model::flipAllActiveTis(void) {

	Tree *t = getCurTree();
	t->flipAllActiveTis();
	
}



Parm* Model::getCurParameter(int i) {

	return parms[ activeParm ][ i ];
	
}



BaseFreqs* Model::getCurBaseFreqs(void) {

	int numParms = parms[0].size();
	for (int i=0; i<numParms; i++)
		{
		Parm *p = parms[ activeParm ][i];
		BaseFreqs *derivedPtr = dynamic_cast<BaseFreqs *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
		}
	cout << "Problem finding BaseFreqs" << endl;
	return NULL;
	
}



SubRates* Model::getCurSubRates(void) {

	int numParms = parms[0].size();
	for (int i=0; i<numParms; i++)
		{
		Parm *p = parms[ activeParm ][i];
		SubRates *derivedPtr = dynamic_cast<SubRates *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
		}
	cout << "Problem finding SubRates" << endl;
	return NULL;

}



GammaShape* Model::getCurGammaShape(void) {

	int numParms = parms[0].size();
	for (int i=0; i<numParms; i++)
		{
		Parm *p = parms[ activeParm ][i];
		GammaShape *derivedPtr = dynamic_cast<GammaShape *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
		}
	cout << "Problem finding GammaShape" << endl;
	return NULL;

}



Tree* Model::getCurTree(void) {

	int numParms = parms[0].size();
	for (int i=0; i<numParms; i++)
		{
		Parm *p = parms[ activeParm ][i];
		Tree *derivedPtr = dynamic_cast<Tree *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
		}
	cout << "Problem finding Tree" << endl;
	return NULL;

}



double Model::lnLikelihood(void) {

	Tree *t = getCurTree();
	int nK = getCurGammaShape()->getNumCats();
	int nC = alignmentPtr->getNumChar();
	int nS = condLikes->getNumStates();
	MbMatrix<double> *lftP = new MbMatrix<double>[3 * nK];
	MbMatrix<double> *rhtP = lftP + nK;
	MbMatrix<double> *ancP = lftP + (2 * nK);
	
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node *p = t->getDownPassNode(n);
		
		if ( p->getUpdateCl() == false )
			continue;
		
		if ( p->getIsLeaf() == false )
			{
			if ( t->getRootPtr()->getLft() == p )
				{
				/* basal three-way split */
				int lftIdx = p->getLft()->getIndex();
				int rhtIdx = p->getRht()->getIndex();
				int ancIdx = p->getAnc()->getIndex();
				int theIdx = p->getIndex();
				for (int k=0; k<nK; k++)
					{
					lftP[k] = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx, k );
					rhtP[k] = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx, k );
					ancP[k] = tiProbs->getTiMatrix( p->getActiveTi(),           theIdx, k );
					}
				double *lftCl = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, 0 );
				double *rhtCl = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, 0 );
				double *ancCl = condLikes->getClPtr( p->getAnc()->getActiveCl(), ancIdx, 0 );
				double *theCl = condLikes->getClPtr( p->getActiveCl(),           theIdx, 0 );
				for (int c=0; c<nC; c++)
					{
					for (int k=0; k<nK; k++)
						{
						for (int i=0; i<nS; i++)
							{
							double sumL=0.0, sumR=0.0, sumA=0.0;
							for (int j=0; j<nS; j++)
								{
								sumL += lftCl[j] * lftP[k][i][j];
								sumR += rhtCl[j] * rhtP[k][i][j];
								sumA += ancCl[j] * ancP[k][i][j];
								}
							theCl[i] = sumL * sumR * sumA;
							}
						lftCl += nS;
						rhtCl += nS;
						ancCl += nS;
						theCl += nS;
						}
					}
				
				}
			else
				{
				/* interior node, dealing with only left and right */
				int lftIdx = p->getLft()->getIndex();
				int rhtIdx = p->getRht()->getIndex();
				int theIdx = p->getIndex();
				for (int k=0; k<nK; k++)
					{
					lftP[k] = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx, k );
					rhtP[k] = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx, k );
					}
				double *lftCl = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, 0 );
				double *rhtCl = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, 0 );
				double *theCl = condLikes->getClPtr( p->getActiveCl(),           theIdx, 0 );
				for (int c=0; c<nC; c++)
					{
					for (int k=0; k<nK; k++)
						{
						for (int i=0; i<nS; i++)
							{
							double sumL=0.0, sumR=0.0;
							for (int j=0; j<nS; j++)
								{
								sumL += lftCl[j] * lftP[k][i][j];
								sumR += rhtCl[j] * rhtP[k][i][j];
								}
							theCl[i] = sumL * sumR;
							}
						lftCl += nS;
						rhtCl += nS;
						theCl += nS;
						}
					}
				/* end, interior node conditional likelihood calculation */
				}
			}
			
		p->setUpdateCl( false );
		}

	delete [] lftP;
	
	Node *b = t->getRootPtr()->getLft();
	int theIdx = b->getIndex();
	double *theCl = condLikes->getClPtr( b->getActiveCl(), theIdx, 0 );
	double rateCatProb = 1.0 / nK;
	MbVector<double> bf = getCurBaseFreqs()->getVal();
	double logLike = 0.0;
	for (int c=0; c<nC; c++)
		{
		double like=0.0;
		for (int k=0; k<nK; k++)
			{
			for (int i=0; i<nS; i++)
				{
				like += theCl[i] * bf[i] * rateCatProb;
				}
			theCl += nS;
			}
		logLike += log(like) * condLikes->getNumOfPattern(c);
		}
		
	return logLike;
	
}



int Model::pickParm(void) {

	double u = ranPtr->uniformRv();
	double sum = 0.0;
	for (int i=0; i<parms[0].size(); i++)
		{
		sum += parmProb[i];
		if (u < sum)
			return i;
		}
	cout << "Problem picking a parameter at random" << endl;
	return -1;
	
}



void Model::upDateTiProbs(void) {

	tiProbs->upDateTiProbs( getCurTree(), getCurGammaShape() );

}



void Model::upDateQ(Parm *p) {

	{
	BaseFreqs *derivedPtr = dynamic_cast<BaseFreqs *>(p);
	if ( derivedPtr != 0 )
		{
		tiMatrix->updateQ( getCurSubRates()->getVal(), derivedPtr->getVal() );
		goto exitFunction;
		}
	}
	{
	SubRates *derivedPtr = dynamic_cast<SubRates *>(p);
	if ( derivedPtr != 0 )
		{
		tiMatrix->updateQ( derivedPtr->getVal(), getCurBaseFreqs()->getVal() );
		goto exitFunction;
		}
	}

	exitFunction:
		;
		
}


void Model::restoreQ(Parm *p) {

	{
	BaseFreqs *derivedPtr = dynamic_cast<BaseFreqs *>(p);
	if ( derivedPtr != 0 )
		{
		tiMatrix->restoreQ();
		goto exitFunction;
		}
	}
	{
	SubRates *derivedPtr = dynamic_cast<SubRates *>(p);
	if ( derivedPtr != 0 )
		{
		tiMatrix->restoreQ();
		goto exitFunction;
		}
	}

	exitFunction:
		;
		
}



void Model::setAllUpdateCls(bool tf) {

	Tree *t = getCurTree();
	t->setAllUpdateCls( tf );

}



void Model::setAllUpdateTis(bool tf) {

	Tree *t = getCurTree();
	t->setAllUpdateTis( tf );

}




