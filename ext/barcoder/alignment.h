#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>
#include <vector>

using namespace std;

class Alignment {

	public:
                	              Alignment(string filePath, string fileName);  
                	              ~Alignment(void);
                	       void   compress(void);
                	        int   getNumTaxa(void) { return numTaxa; }
                	        int   getNumChar(void) { return (compressedData == true ? numSitePatterns : numChar); }
                	       void   getPossibleNucs (int nucCode, int *nuc);
                	        int   getNucleotide(int i, int j);
							int   getTaxonIndex(string ns);
                	       void   listTaxa(void);
                	       void   print(void);
                	       void   uncompress(void);
                	     string   getTaxonName(int i);
				 vector<string>   getTaxonNames(void) { return taxonNames; }
						    int   getNumOfPattern(int i) { return patternCount[i]; }
                	              
	private:
                            int   numTaxa;
                            int   numChar;
                            int   numSitePatterns;
                 vector<string>   taxonNames;
                           bool   compressedData;
                            int   nucID(char nuc);
                            int   **matrix;
                            int   **compressedMatrix;
                            int   *patternCount;
                            
};

#endif
