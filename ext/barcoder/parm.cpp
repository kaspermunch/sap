#include "parm.h"
#include "parm_subrates.h"
#include "parm_shape.h"
#include "parm_basefreqs.h"
#include "parm_tree.h"



Parm::Parm(MbRandom *rn, Model *mp) {

	ranPtr   = rn;
	modelPtr = mp;
	
}


Parm &Parm::operator=(Parm &b) {

	if (this != &b)
		{
		/* copy base class data */
		ranPtr   = b.ranPtr;
		modelPtr = b.modelPtr;
		
		/* We need to downcast the object to the derived class pointers.
		   This allows us to call the correct assignment functions in the
		   derived class. */
		
		/* check to see if the object is a substitution rate parameter */
		{
		SubRates *thatDerivedPtr = dynamic_cast<SubRates *>(&b);
		SubRates *thisDerivedPtr = dynamic_cast<SubRates *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}

		/* check to see if the object is a gamma shape parameter */
		{
		GammaShape *thatDerivedPtr = dynamic_cast<GammaShape *>(&b);
		GammaShape *thisDerivedPtr = dynamic_cast<GammaShape *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}

		/* check to see if the object is a base frequency parameter */
		{
		BaseFreqs *thatDerivedPtr = dynamic_cast<BaseFreqs *>(&b);
		BaseFreqs *thisDerivedPtr = dynamic_cast<BaseFreqs *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}

		/* check to see if the object is a tree parameter */
		{
		Tree *thatDerivedPtr = dynamic_cast<Tree *>(&b);
		Tree *thisDerivedPtr = dynamic_cast<Tree *>(this);
		if ( thatDerivedPtr != 0 && thisDerivedPtr != 0 )
			{
			thisDerivedPtr->clone( *thatDerivedPtr );
			goto exitOperator;
			}
		}
		
		cout << "Problem in Parameter assignment operator" << endl;
		exit(1);
			
		exitOperator:
			;
		}
	return *this;

}
