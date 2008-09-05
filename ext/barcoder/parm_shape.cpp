#include <iostream>
#include <cmath>
#include "parm_shape.h"
#include "parm.h"
#include "model.h"
#include "settings.h"
#include "MbRandom.h"


GammaShape::GammaShape(MbRandom *rn, Model *mp, double tn, double lam, int nc, double x) : Parm(rn, mp) {

	parmName = "gamma shape parameter";
	alpha     = x;
	lambda    = lam;
	tuning    = tn;
	numCats   = nc;
	alpha = ranPtr->exponentialRv(lambda);
	r = MbVector<double>(numCats);
	if (numCats == 1)
		r[0] = 1.0;
	else
		ranPtr->discretizeGamma( r, alpha, alpha, numCats, false );

}



GammaShape::GammaShape(GammaShape &b) : Parm(b.ranPtr, b.modelPtr) {

	alpha   = b.alpha;
	lambda  = b.lambda;
	tuning  = b.tuning;
	numCats = b.numCats;
	r = MbVector<double>(numCats);
	for (int i=0; i<numCats; i++)
		r[i] = b.r[i];
	
}



double GammaShape::change(void) {

	modelPtr->setAllUpdateCls( true );
	modelPtr->setAllUpdateTis( true );
	modelPtr->flipAllActiveCls();
	modelPtr->flipAllActiveTis();

	double oldAlpha = alpha;
	alpha = oldAlpha * exp( tuning*(ranPtr->uniformRv()-0.5) );
	if (numCats == 1)
		r[0] = 1.0;
	else
		ranPtr->discretizeGamma( r, alpha, alpha, numCats, false );
	return log(alpha) - log(oldAlpha);
	
}



double GammaShape::getLnPriorProbability(void) {

	return ranPtr->lnExponentialPdf(lambda, alpha);

}



void GammaShape::print(void) {

	std::cout << "Gamma shape        = " << alpha << std::endl;;

}



GammaShape &GammaShape::operator=(GammaShape &b) {

	if (this != &b)
		{
		Parm::operator=(b);
		}
	return *this;

}



void GammaShape::clone(GammaShape &b) {

	alpha   = b.alpha;
	lambda  = b.lambda;
	tuning  = b.tuning;
	numCats = b.numCats;
	r.inject(b.r);
	
}
