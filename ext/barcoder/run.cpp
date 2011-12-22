#include <iostream>
#include "alignment.h"
#include "constraints.h"
#include "iomanager.h"
#include "settings.h"
#include "MbRandom.h"
#include "mcmc.h"
#include "model.h"
#include "parm_tree.h"

#include <cstdio>
#include <cstdlib> 
#include <cstring>

using namespace std;

extern "C" int runprogram (int argc, char **argv, char *outputBaseName);

int runprogram (int argc, char **argv, char *outputBaseName) {


        /* Redirect cout and cerr to files: ************/
        char *stdoutFile = NULL;
        char *stderrFile = NULL;

        stdoutFile = (char *) calloc(strlen(outputBaseName) + 5, sizeof(char));
        sprintf(stdoutFile,"%s.out", outputBaseName);  
        stderrFile = (char *) calloc(strlen(outputBaseName) + 5, sizeof(char));
        sprintf(stderrFile,"%s.err", outputBaseName);  

	ofstream outFile;
	outFile.open (stdoutFile);
	streambuf* savedStdout = cout.rdbuf();
	std::cout.rdbuf(outFile.rdbuf());

	ofstream errFile;
	errFile.open (stderrFile);
	streambuf* savedStderr = cerr.rdbuf();
	std::cerr.rdbuf(errFile.rdbuf());

	/* This does not redirect printf(). That can be done with
	   freopen but I don't know how the recover the old stdout and
	   stderr later */
	/* 
	std::freopen(stdoutFile, "w", stdout);
	std::freopen(stderrFile, "w", stderr);
	*/
	/*************************************************/

	cout << endl;
	cout << "   Barcoder 1.0" << endl << endl;
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
	MbRandom myRandom;
	
	/* instantiate the phylogenetic model */
	Model myModel( &myRandom, &alignment, &myConstraints, &mySettings );
	
	/* run the Markov chain */
	Mcmc markovChain( &myModel, &myRandom, &mySettings, &outputFileInfo );


	/* Reset cout and cerr *****************/
	cout.rdbuf(savedStdout);
	cerr.rdbuf(savedStderr);
	/***************************************/
		
    return 0;

}
