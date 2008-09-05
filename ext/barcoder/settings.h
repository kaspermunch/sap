#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>

using namespace std;

class Settings {

	public:
                            Settings(int argc, char *argv[]);
				   double   getBrlenLambda(void) { return brlenLambda; }
				      int   getChainLength(void) { return chainLength; }
				   string   getConstraintFilePathName(void) { return constraintFilePathName; }
				   string   getDataFilePathName(void) { return dataFilePathName; }
				      int   getPrintFrequency(void) { return printFrequency; }
				      int   getNumGammaCats(void) { return numGammaCats; }
				   string   getOutPutFileName(void) { return outPutFileName; }
					  int   getSampleFrequency(void) { return sampleFrequency; }
					 void   setBrlenLambda(double x) { brlenLambda = x; }
					 void   setConstraintFilePathName(string s) { constraintFilePathName = s; }
					 void   setDataFilePathName(string s) { dataFilePathName = s; }
					 void   setNumGammaCats(int x) { numGammaCats = x; }
					 void   setOutPutFileName(string s) { outPutFileName = s; }

	private:
				   double   brlenLambda;
				      int   chainLength;
				   string   constraintFilePathName;
                   string   dataFilePathName;
				      int   numGammaCats;
				   string   outPutFileName;
					  int   printFrequency;
					  int   sampleFrequency;

};


#endif
