#include "settings.h"
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>

using namespace std;


Settings::Settings(int argc, char *argv[]) {

#	if 0
	/* set up fake command-line argument string */
	char *cmdString[9];
	cmdString[0] = "bc";
	cmdString[1] = "-i";
	cmdString[2] = "/Users/johnh/Desktop/barcoder2/sample.phylip";
	cmdString[3] = "-o";
	cmdString[4] = "/Users/johnh/Desktop/barcoder2/sample.out";
	cmdString[5] = "-c";
	cmdString[6] = "/Users/johnh/Desktop/barcoder2/sample.constr";
	cmdString[7] = "-l";
	cmdString[8] = "1000";
	argc = 9;
	argv = cmdString;
#	endif

	enum Mode { CONSTRAINTS_FILE, DATA_FILE, OUTPUT_FILE, CHAIN_LENGTH, PRINT_FREQ, SAMPLE_FREQ, BRLEN_PARM, NUM_GAMMA_CATS, NONE };

	/* set default values for parameters */
	constraintFilePathName = "";
	dataFilePathName       = "";
	outPutFileName         = "";
	chainLength            = 1000000;
	printFrequency         = 100;
	sampleFrequency        = 100;
	brlenLambda            = 10.0;
	numGammaCats           = 4;
	
	if (argc > 1)
		{
		if (argc % 2 == 0)
			{
			cout << "Usage:" << endl;
			cout << "   -i : Input file name" << endl;
			cout << "   -c : Constraints file name" << endl;
			cout << "   -o : Output file name" << endl;
			cout << "   -l : Number of MCMC cycles" << endl;
			cout << "   -p : Print frequency" << endl;
			cout << "   -s : Sample frequency" << endl;
			cout << "   -b : Exponential parameter for branch lengths" << endl;
			cout << "   -g : Number of gamma rate categories" << endl << endl;
			cout << "Example:" << endl;
			cout << "   bc1 -i <input file> -c <constraints file> -o <output file>" << endl;
			exit(0);
			}
			
		/* read the command-line arguments */
		int status = NONE;
		for (int i=1; i<argc; i++)
			{
			string cmd = argv[i];
			//cout << cmd << endl;
			if (status == NONE)
				{
				/* read the parameter specifier */
				if ( cmd == "-i" )
					status = DATA_FILE;
				else if ( cmd == "-c" )
					status = CONSTRAINTS_FILE;
				else if ( cmd == "-o" )
					status = OUTPUT_FILE;
				else if ( cmd == "-l" )
					status = CHAIN_LENGTH;
				else if ( cmd == "-p" )
					status = PRINT_FREQ;
				else if ( cmd == "-s" )
					status = SAMPLE_FREQ;
				else if ( cmd == "-b" )
					status = BRLEN_PARM;
				else if ( cmd == "-g" )
					status = NUM_GAMMA_CATS;
				else
					{
					cerr << "Could not interpret option \"" << cmd << "\"." << endl;
					exit(1);
					}
				}
			else
				{
				/* read the parameter */
				if ( status == DATA_FILE )
					dataFilePathName = argv[i];
				else if ( status == CONSTRAINTS_FILE )
					constraintFilePathName = argv[i];
				else if ( status == OUTPUT_FILE )
					outPutFileName = argv[i];
				else if ( status == CHAIN_LENGTH )
					chainLength = atoi(argv[i]);
				else if ( status == PRINT_FREQ )
					printFrequency = atoi(argv[i]);
				else if ( status == SAMPLE_FREQ )
					sampleFrequency = atoi(argv[i]);
				else if ( status == BRLEN_PARM )
					brlenLambda = atof(argv[i]);
				else if ( status == NUM_GAMMA_CATS )
					numGammaCats = atoi(argv[i]);
				else
					{
					cerr << "Unknown status reading command line information" << endl;
					exit(1);
					}
				status = NONE;
				}
			}
		}
	else
		{
		/* set the parameters from the standard in */
		cout << "   Program settings (Note: This program supports command line specification of settings)" << endl;
		cout << "      Data file name        = ";
		cin >> dataFilePathName;
		cout << "      Constraint file name  = ";
		cin >> constraintFilePathName;
		cout << "      Output file name      = ";
		cin >> outPutFileName;
		cout << "      Chain length          = ";
		cin >> chainLength;
		cout << "      Print frequency       = ";
		cin >> printFrequency;
		cout << "      Sample frequency      = ";
		cin >> sampleFrequency;
		cout << endl;
		}	

}
