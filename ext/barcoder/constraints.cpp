#include "constraints.h"
#include "alignment.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <fstream>
#include "MbBitfield.h"
#include "crossplatform.h"

using namespace std;



Bipartition::Bipartition(int nt) {
	
	numTaxa = nt;
	parts = new MbBitfield(numTaxa);
	mask  = new MbBitfield(numTaxa);
	
}



Bipartition::~Bipartition(void) {

	delete parts;
	delete mask;
	
}



void Bipartition::addToBipartition(int tn) {

	parts->setBit( tn );
	
}



void Bipartition::addToMask(int tn) {

	mask->setBit( tn );
	
}



int Bipartition::getBipartition(int tn) {

	if (parts->isBitSet( tn ) == true)
		return 1;
	else
		return 0;
	
}



int Bipartition::getMask(int tn) {

	if (mask->isBitSet( tn ) == true)
		return 1;
	else
		return 0;
	
}



void Bipartition::print(void) {

	cout << "Part -- " << (*parts) << endl;
	cout << "Mask -- " << (*mask) << endl;
	
}



void Bipartition::normalizeBipartition(void) {

#	if 0
	if (parts->isBitSet(0) == true)
		{
		parts->flipBits();
		(*parts) &= (*mask);
		}
#	else
	if (parts->isBitSet(0) == true)
		{
		for (int i=0; i<parts->dim(); i++)
			{
			if (mask->isBitSet(i) == true)
				parts->flipBit(i);
			}
		}
#	endif
}



Constraints::Constraints(string filePath, string fileName, Alignment *al) {

	alignmentPtr = al;
	numTaxa      = alignmentPtr->getNumTaxa();
	
	string filePathAndName = filePath + PATH_SEPERATOR + fileName;
	
	cout << "   Reading constraints" << endl;
	if (fileName == "")
		{
		cout << "      Constraints file name = NO FILE TO READ" << fileName << endl;
		cout << "      Number of constraints = 0" << endl << endl;
		return;
		}
	else
		cout << "      Constraints file name = " << fileName << endl;

	/* open the file */
	ifstream myFileStream(filePathAndName.c_str());
	if (!myFileStream) 
		{
		cerr << "      ERROR: Cannot open file \"" + fileName + "\"" << endl;
		exit(1);
		}

	// read the constraints, one line at a time
	string word = "";
	char ch;
	int line = 0;
	int whichSide = 0;
	int wordNum = 0;
	int numOnSide[2] = { 0, 0 };
	int numConstraints = 0;
	while ( (ch = myFileStream.get()) != EOF )
		{
		//cout << ch;
		if ( ch == ' ' || ch == '\t' || ch == '\n' )
			{
			if (wordNum == 0)
				{
				biparts.push_back( new Bipartition(numTaxa) );
				numConstraints++;
				}
			/* interpret the word */
			if ( word == "|" )
				{
				whichSide++;
				if (whichSide > 1)
					{
					cerr << "Problem reading constraint" << endl;
					exit (1);
					}
				}
			else
				{
				int whichTaxon = alignmentPtr->getTaxonIndex( word );
				if (whichSide == 0)
					{
					//cout << "adding " << whichTaxon << " to bipartition " << biparts[line] << endl;
					biparts[line]->addToBipartition(whichTaxon);
					}
				biparts[line]->addToMask(whichTaxon);
				numOnSide[whichSide]++;
				}
			//cout << line << " " << wordNum << " word = " << word << endl;
			word = "";
			wordNum++;
			}
		else
			{
			word += ch;
			}
		if (ch == '\n')
			{
			/* end of line */
			line++;
			whichSide = 0;
			wordNum = 0;
			if (numOnSide[0] < 2 || numOnSide[1] < 2)
				{
				numConstraints--;
				cout << "      WARNING: Ignoring constraint on line " << line << ", which is consistent with all trees" << endl;
				Bipartition *bp = biparts[numConstraints];
				delete bp;
				biparts.pop_back();
				}
			}
		}
		
	// close the file
	myFileStream.close();

	// normalize the taxon bipartitions
	for (vector<Bipartition *>::iterator p=biparts.begin(); p != biparts.end(); p++)
		(*p)->normalizeBipartition();
		
	cout << "      Number of constraints = " << biparts.size() << endl << endl;

#	if 0
	print();
#	endif

}



Constraints::~Constraints(void) {

	for (vector<Bipartition *>::iterator p=biparts.begin(); p != biparts.end(); p++)
		delete (*p);
		
}



void Constraints::print(void) {

	cout << "Constraints -- " << endl;
	int i = 0;
	for (vector<Bipartition *>::iterator p=biparts.begin(); p != biparts.end(); p++)
		{
		cout << "Bipartition " << ++i << endl;
		(*p)->print();
		}

}





