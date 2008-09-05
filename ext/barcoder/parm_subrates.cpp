#include <iostream>
#include "model.h"
#include "parm_subrates.h"
#include "parm.h"
#include "model.h"
#include "settings.h"
#include "MbRandom.h"
#include "MbVector.h"



SubRates::SubRates(MbRandom *rn, Model *mp, int ns, double a0) : Parm(rn, mp) {

	parmName = "substitution rate parameters";
	numStates = ns;
	numRates  = numStates * (numStates - 1) / 2;
	alpha0    = a0;
	rates     = MbVector<double>(numRates);
	a         = MbVector<double>(numRates);
	dirParm   = MbVector<double>(numRates);
	for (int i=0; i<numRates; i++)
		{
		rates[i]   = 1.0/numRates;
		a[i]       = 1.0;
		dirParm[i] = 1.0;
		}
#	if 0
	if ( isFixed == false )
		rnd->dirichletRv(a, rates);
#	endif

}



SubRates::SubRates(SubRates &b) : Parm(b.ranPtr, b.modelPtr) {

	numStates = b.numStates;
	numRates  = b.numRates;
	alpha0    = b.alpha0;
	rates     = MbVector<double>(numRates);
	a         = MbVector<double>(numRates);
	dirParm   = MbVector<double>(numRates);
	for (int i=0; i<numRates; i++)
		{
		rates[i]   = b.rates[i];
		a[i]       = b.a[i];
		dirParm[i] = b.dirParm[i];
		}
	
}



SubRates &SubRates::operator=(SubRates &b) {

	if (this != &b)
		{
		Parm::operator=(b);
		}
	return *this;

}



double SubRates::change(void) {

	modelPtr->setAllUpdateCls( true );
	modelPtr->setAllUpdateTis( true );
	modelPtr->flipAllActiveCls();
	modelPtr->flipAllActiveTis();

	MbVector<double> srNew(numRates);
	MbVector<double> srOld(numRates);
	MbVector<double> aNew(numRates);
	MbVector<double> aOld(numRates);
	for (int i=0; i<numRates; i++)
		{
		aNew[i] = rates[i] * alpha0;
		srOld[i] = rates[i];
		}
	ranPtr->dirichletRv(aNew, srNew);
	double sum = 0.0;
	for (int i=0; i<numRates; i++)
		{
		if (srNew[i] < 0.001)
			srNew[i] = 0.001;
		sum += srNew[i];
		}
	for (int i=0; i<numRates; i++)
		srNew[i] /= sum;
	for (int i=0; i<numRates; i++)
		{
		aOld[i] = srNew[i] * alpha0;
		rates[i] = srNew[i];
		}
	return ranPtr->lnDirichletPdf(aNew, srNew) - ranPtr->lnDirichletPdf(aOld, srOld);
	
}



double SubRates::getLnPriorProbability(void) {

	return ranPtr->lnDirichletPdf(dirParm, rates);

}



void SubRates::print(void) {

	std::cout << "Substitution rates = (";
	for (int i=0; i<numRates; i++)
		{
		std::cout << rates[i];
		if (i+1 != numRates)
			std::cout << ",";
		}
	std::cout << ")" << std::endl;

}



void SubRates::clone(SubRates &b) {

	numStates = b.numStates;
	numRates  = b.numRates;
	alpha0    = b.alpha0;
	rates.inject(b.rates);
	a.inject(b.a);
	dirParm.inject(b.dirParm);
	
}
