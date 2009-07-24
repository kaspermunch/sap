#include <iostream>
#include "alignment.h"
#include "constraints.h"
#include "iomanager.h"
#include "settings.h"
#include "MbRandom.h"
#include "mcmc.h"
#include "model.h"
#include "parm_tree.h"

using namespace std;



int main (int argc, char *argv[]) {

	cout << endl;
	cout << "   Barcoder 1.01" << endl << endl;
	cout << "   John P. Huelsenbeck" << endl;
	cout << "   Department of Integrative Biology" << endl;
	cout << "   University of California, Berkeley" << endl << endl;
	
	/* set parameters */
	Settings mySettings( argc, argv );
	
	/* initialize file managers */
	IoManager dataFileInfo( mySettings.getDataFilePathName() );
	IoManager constraintFileInfo( mySettings.getConstraintFilePathName() );
	IoManager outputFileInfo( mySettings.getOutPutFileName() );

	/* read the sequence data */
	Alignment alignment( dataFileInfo.getFilePath(), dataFileInfo.getFileName() );
	alignment.compress();

	/* read the constraints */
	Constraints myConstraints( constraintFileInfo.getFilePath(), constraintFileInfo.getFileName(), &alignment );
	
	/* instantiate random number generator */
	MbRandom myRandom(3);
	
	/* instantiate the phylogenetic model */
	Model myModel( &myRandom, &alignment, &myConstraints, &mySettings );
	
	/* run the Markov chain */
	Mcmc markovChain( &myModel, &myRandom, &mySettings, &outputFileInfo );
		
    return 0;

}
