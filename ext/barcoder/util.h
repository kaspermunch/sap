#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <string>
#include <fstream>

using namespace std;



string getLineFromFile(string fileName, int lineNum) {

	/* open the file */
	ifstream fileStream(fileName.c_str());
	if (!fileStream) 
		{
		cerr << "Cannot open file \"" + fileName + "\"" << endl;
		exit(1);
		}

	string linestring = "";
	int line = 0;
	while( getline(fileStream, linestring).good() )
		{
		line++;
		if (line == lineNum)
			break;
		}

	/* close the file */
	fileStream.close();
	
	if (line != lineNum)
		{
		cerr << "The file \"" + fileName + "\" has " << line << " lines. Could not find line " << lineNum << endl;
		exit(1);
		}
	
	return linestring;

}



int flip(int x) {
	
	if (x == 0)
		return 1;
	else
		return 0;
		
}
#endif