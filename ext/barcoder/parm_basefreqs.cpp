#include <iostream>
#include "parm_basefreqs.h"
#include "parm.h"
#include "model.h"
#include "settings.h"
#include "MbRandom.h"



BaseFreqs::BaseFreqs(MbRandom *rn, Model *mp, int ns, double a0) : Parm(rn, mp) {

	parmName = "base frequencies";
	numStates = ns;
	alpha0    = a0;
	freqs     = MbVector<double>(numStates);
	a         = MbVector<double>(numStates);
	dirParm   = MbVector<double>(numStates);
	for (int i=0; i<numStates; i++)
		{
		freqs[i]   = 1.0/numStates;
		a[i]       = 1.0;
		dirParm[i] = 1.0;
		}
#	if 0
	rnd->dirichletRv(a, freqs);
#	endif

}



BaseFreqs::BaseFreqs(BaseFreqs &b) : Parm(b.ranPtr, b.modelPtr) {

	numStates = b.numStates;
	alpha0    = b.alpha0;
	freqs     = MbVector<double>(numStates);
	a         = MbVector<double>(numStates);
	dirParm   = MbVector<double>(numStates);
	for (int i=0; i<numStates; i++)
		{
		freqs[i]   = b.freqs[i];
		a[i]       = b.a[i];
		dirParm[i] = b.dirParm[i];
		}
	
}



double BaseFreqs::change(void) {
	
	/* set all of the flags */
	modelPtr->setAllUpdateCls( true );
	modelPtr->setAllUpdateTis( true );
	modelPtr->flipAllActiveCls();
	modelPtr->flipAllActiveTis();
	
	/* change the state */
	MbVector<double> bfNew(numStates);
	MbVector<double> bfOld(numStates);
	MbVector<double> aNew(numStates);
	MbVector<double> aOld(numStates);
	for (int i=0; i<numStates; i++)
		{
		aNew[i] = freqs[i] * alpha0;
		bfOld[i] = freqs[i];
		}
	ranPtr->dirichletRv(aNew, bfNew);
	double sum = 0.0;
	for (int i=0; i<numStates; i++)
		{
		if (bfNew[i] < 0.001)
			bfNew[i] = 0.001;
		sum += bfNew[i];
		}
	for (int i=0; i<numStates; i++)
		bfNew[i] /= sum;
	for (int i=0; i<numStates; i++)
		{
		aOld[i] = bfNew[i] * alpha0;
		freqs[i] = bfNew[i];
		}
	return ranPtr->lnDirichletPdf(aNew, bfNew) - ranPtr->lnDirichletPdf(aOld, bfOld);
	
}



double BaseFreqs::getLnPriorProbability(void) {

	return ranPtr->lnDirichletPdf(dirParm, freqs);

}



void BaseFreqs::print(void) {

	std::cout << "Base frequencies   = " << freqs << std::endl;

}



BaseFreqs &BaseFreqs::operator=(BaseFreqs &b) {

	if (this != &b)
		{
		cout << "assignment operator" << endl;
		Parm::operator=(b);
		}
	return *this;

}



void BaseFreqs::clone(BaseFreqs &b) {

	numStates = b.numStates;
	alpha0    = b.alpha0;
	freqs.inject(b.freqs);
	a.inject(b.a);
	dirParm.inject(b.dirParm);
	
}

