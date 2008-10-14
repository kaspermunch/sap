#include "alignment.h"
#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <fstream>
#include "crossplatform.h"



Alignment::Alignment(string filePath, string fileName) {

	string filePathAndName = filePath + PATH_SEPERATOR + fileName;

	/* open the file */
	ifstream seqStream(filePathAndName.c_str());
	if (!seqStream) 
		{
		cerr << "Cannot open file \"" + fileName + "\"" << endl;
		exit(1);
		}
		
	string linestring = "";
	int line = 0;
	string theSequence = "";
	int taxonNum = 0;
	matrix = NULL;
	numTaxa = numChar = 0;
	while( getline(seqStream, linestring).good() )
		{
		istringstream linestream(linestring);
		int ch;
		string word = "";
		int wordNum = 0;
		int siteNum = 0;
		do
			{
			word = "";
			linestream >> word;
			wordNum++;
			//cout << "word(" << wordNum << ") = " << word << endl;
			if (line == 0)
				{
				/* read the number of taxa/chars from the first line */
				int x;
				istringstream buf(word);
				buf >> x;
				if (wordNum == 1)
					numTaxa = x;
				else
					numChar = numSitePatterns = x;
				if (numTaxa > 0 && numChar > 0 && matrix == NULL)
					{	
					matrix = new int*[numTaxa];
					matrix[0] = new int[numTaxa * numChar];
					for (int i=1; i<numTaxa; i++)
						matrix[i] = matrix[i-1] + numChar;
					for (int i=0; i<numTaxa; i++)
						for (int j=0; j<numChar; j++)
							matrix[i][j] = 0;
					patternCount = new int[numChar];
					for (int i=0; i<numChar; i++)
						patternCount[i] = 1;
					compressedData = false;
					}
				}
			else
				{
				if (wordNum == 1)
					{
					taxonNames.push_back(word);
					taxonNum++;
					}
				else
					{
					for (int i=0; i<word.length(); i++)
						{
						char site = word.at(i);
						matrix[taxonNum-1][siteNum++] = nucID(site);
						}
					}
				}
			} while ( (ch=linestream.get()) != EOF );
			
		// NOTE: We probably do not need this bit of code.
		if (line == 0)
			{
			/* the first line should contain the number of taxa and the sequence length */
			istringstream buf(word);
			//buf >> genomeSize;
			}
		else
			{
			for (int i=0; i<word.length(); i++)
				{
				char site = word.at(i);
				if (tolower(site) == 'a' || tolower(site) == 'c' || tolower(site) == 'g' || tolower(site) == 't')
					theSequence += tolower(site);
				}
			}
		//cout << linestring << endl;
		line++;
		}	

	/* close the file */
	seqStream.close();

	cout << "   Reading alignment file" << endl;
	cout << "      File name             = " << fileName << endl;
	cout << "      Number of taxa        = " << numTaxa << endl;
	cout << "      Number of sites       = " << numChar << endl << endl;

}



Alignment::~Alignment(void) {

	delete [] matrix[0];
	delete [] matrix;
	delete [] patternCount;
	if (compressedData == true)
		{
		delete [] compressedMatrix[0];
		delete [] compressedMatrix;
		}

}



void Alignment::compress(void) {

	if (compressedData == false)
		{
		int *tempCount = new int[numChar];
		for (int j=0; j<numChar; j++)
			tempCount[j] = 1;
		
		for (int j=0; j<numChar; j++)
			{
			if (tempCount[j] > 0)
				{
				for (int k=j+1; k<numChar; k++)
					{
					if (tempCount[k] == 0)
						continue;
					bool isSame = true;
					for (int i=0; i<numTaxa; i++)
						{
						if (matrix[i][j] != matrix[i][k])
							{
							isSame = false;
							break;
							}
						}
					if (isSame == true)
						{
						tempCount[j]++;
						tempCount[k] = 0;
						}
					}
				}
			}
		numSitePatterns = 0;
		for (int j=0; j<numChar; j++)
			if (tempCount[j] > 0)
				numSitePatterns++;
				
		compressedMatrix = new int*[numTaxa];
		compressedMatrix[0] = new int[numTaxa * numSitePatterns];
		for (int i=1; i<numTaxa; i++)
			compressedMatrix[i] = compressedMatrix[i-1] + numSitePatterns;
		for (int i=0; i<numTaxa; i++)
			for (int j=0; j<numSitePatterns; j++)
				compressedMatrix[i][j] = 0;
		
		delete [] patternCount;
		patternCount = new int[numSitePatterns];
		for (int j=0, k=0; j<numChar; j++)
			{
			if (tempCount[j] > 0)
				{
				for (int i=0; i<numTaxa; i++)
					compressedMatrix[i][k] = matrix[i][j];
				patternCount[k] = tempCount[j];
				k++;
				}
			}
		compressedData = true;
		delete [] tempCount;
		}
		
}



int Alignment::getNucleotide(int i, int j) {

	if (compressedData == true)
		return compressedMatrix[i][j];
	else
		return matrix[i][j];
		
}



/*-------------------------------------------------------------------
|
|   GetPossibleNucs: 
|
|   This function initializes a vector, nuc[MAX_NUM_STATES]. The four elements
|   of nuc correspond to the four nucleotides in alphabetical order.
|   We are assuming that the nucCode is a binary representation of
|   the nucleotides that are consistent with the observation. For
|   example, if we observe an A, then the nucCode is 1 and the 
|   function initalizes nuc[0] = 1 and the other elements of nuc
|   to be 0.
|
|   Observation    nucCode        nuc
|        A            1           1000
|        C            2           0100
|        G            4           0010
|        T            8           0001
|        R            5           1010
|        Y           10           0101
|        M            3           1100
|        K           12           0011
|        S            6           0110
|        W            9           1001
|        H           11           1101
|        B           14           0111
|        V            7           1110
|        D           13           1011
|        N - ?       15           1111
|
-------------------------------------------------------------------*/
void Alignment::getPossibleNucs (int nucCode, int *nuc) {

	if (nucCode == 1)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 2)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 3)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 4)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 5)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 6)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 7)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 8)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 9)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 10)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 11)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 12)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 13)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 14)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 15)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 16)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}

}



void Alignment::listTaxa(void) {

	int i = 1;
	for (vector<string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		cout << setw(4) << i++ << " -- " << (*p) << endl;

}



string Alignment::getTaxonName(int i) {

	return taxonNames[i];

}



int Alignment::getTaxonIndex(string ns) {

	int taxonIndex = -1;
	int i = 0;
	for (vector<string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		{
		if ( (*p) == ns )
			{
			taxonIndex = i;
			break;
			}
		i++;
		}
	return taxonIndex;

}



/*-------------------------------------------------------------------
|
|   NucID: 
|
|   Take a character, nuc, and return an integer:
|
|       nuc        returns
|        A            1 
|        C            2     
|        G            4      
|        T U          8     
|        R            5      
|        Y           10       
|        M            3      
|        K           12   
|        S            6     
|        W            9      
|        H           11      
|        B           14     
|        V            7      
|        D           13  
|        N - ?       15       
|
-------------------------------------------------------------------*/
int Alignment::nucID(char nuc) {

	char		n;
	
	if (nuc == 'U' || nuc == 'u')
		n = 'T';
	else
		n = nuc;

	if (n == 'A' || n == 'a')
		{
		return 1;
		}
	else if (n == 'C' || n == 'c')
		{
		return 2;
		}
	else if (n == 'G' || n == 'g')
		{
		return 4;
		}
	else if (n == 'T' || n == 't')
		{
		return 8;
		}
	else if (n == 'R' || n == 'r')
		{
		return 5;
		}
	else if (n == 'Y' || n == 'y')
		{
		return 10;
		}
	else if (n == 'M' || n == 'm')
		{
		return 3;
		}
	else if (n == 'K' || n == 'k')
		{
		return 12;
		}
	else if (n == 'S' || n == 's')
		{
		return 6;
		}
	else if (n == 'W' || n == 'w')
		{
		return 9;
		}
	else if (n == 'H' || n == 'h')
		{
		return 11;
		}
	else if (n == 'B' || n == 'b')
		{
		return 14;
		}
	else if (n == 'V' || n == 'v')
		{
		return 7;
		}
	else if (n == 'D' || n == 'd')
		{
		return 13;
		}
	else if (n == 'N' || n == 'n')
		{
		return 15;
		}
	else if (n == '-')
		{
		return 15;
		}
	else if (n == '?')
		{
		return 15;
		}
	else
		return -1;
		
}



void Alignment::print(void) {

	int **x;
	if (compressedData == false)
		x = matrix;
	else
		x = compressedMatrix;
		
	cout << "                ";
	for (int i=0; i<numTaxa; i++)
		cout << setw(3) << i;
	cout << endl;
	cout << "----------------";
	for (int i=0; i<numTaxa; i++)
		cout << "---";
	cout << endl;	
	for (int j=0; j<numSitePatterns; j++)
		{
		cout << setw(4) << j << " -- ";
		cout << setw(4) << patternCount[j] << " -- ";
		for (int i=0; i<numTaxa; i++)
			{
			cout << setw(3) << x[i][j];
			}
		cout << endl;
		}
		
	int sum = 0;
	for (int j=0; j<numSitePatterns; j++)
		sum += patternCount[j];
	cout << "Number of sites = " << sum << endl;
	
}



void Alignment::uncompress(void) {

	if (compressedData == true)
		{
		delete [] compressedMatrix[0];
		delete [] compressedMatrix;
		delete [] patternCount;
		patternCount = new int[numChar];
		for (int j=0; j<numChar; j++)
			patternCount[j] = 1;
		numSitePatterns = numChar;
		compressedData = false;
		}

}

