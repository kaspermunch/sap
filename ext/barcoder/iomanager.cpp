#include "iomanager.h"
#include <iostream>
#include <dirent.h>
#include "crossplatform.h"


IoManager::IoManager(void) {

	setFileName("");
	setFilePath("");
	setCurDirectory( findCurrentDirectory() );
	setFilePath( getCurDirectory() );

}



IoManager::IoManager(string s) {

	setFileName("");
	setFilePath("");
	parsePathFileNames(s);
	setCurDirectory( findCurrentDirectory() );
	if ( getFilePath() == "" )
		setFilePath( getCurDirectory() );
	
}



bool IoManager::parsePathFileNames(string s) {

	string delimiter = PATH_SEPERATOR;

	if ( s.length() == 0)
		{
		/* The string that is supposed to hold the
		   path/file information is empty. Simply
		   return: the constructor calling this
		   function will fill in default values
		   for the path and file name. */
		//cerr << "No path or file name provided" << endl;
		return false;
		}
		
	/* Find the location of the last "/". This is where
	   we will divide the path/file string into two. */
	int location = s.find_last_of( delimiter );
	
	if ( location == -1 )
		{
		/* There is no path in this string. We 
		   must have only the file name, and the
		   file should be in our current directory.*/
		fileName = s;
		filePath = "";
		}
	else if ( location == (int)s.length() - 1 )
		{
		/* It looks like the last character is "/", which
		   means that no file name has been provided. */
		s.erase( location );
		fileName = "";
		filePath = s;
		cerr << "Only a path name was provided" << endl;
		return false;
		}
	else
		{
		/* We can divide the path into the path and the
		   file. */
		fileName = s.substr( location+1, s.length()-location-1 );
		s.erase( location );
		filePath = s;
		}

	return true;

}



#define	MAX_DIR_PATH	2048
string IoManager::findCurrentDirectory(void) {

	string delimiter = PATH_SEPERATOR;

	char cwd[MAX_DIR_PATH+1];
	//if ( !getcwd(cwd, MAX_DIR_PATH+1) )
	if ( !GETCWD(cwd, MAX_DIR_PATH+1) )
		{
		cerr << "Problem finding the current director" << endl;
		return "";
		}
	string curdir = cwd;
	
	if ( curdir.at( curdir.length()-1 ) == delimiter[0] )
		curdir.erase( curdir.length()-1 );
	
	return curdir;

}



string IoManager::getFilePathName(void) {

	string delimiter = PATH_SEPERATOR;

	return filePath + delimiter + fileName;

}



bool IoManager::isDirectoryPresent(const string mp) {

	/* attempt to open the directory */
	DIR *dir = opendir( mp.c_str() );
	if ( !dir )
		return false;
		
	/* close the directory */
	if ( closedir(dir) == -1 )
		cerr << "Problem closing directory" << endl;
		
	return true;

}



bool IoManager::isFilePresent(const string mp, const string mf) {

	/* open the directory */
	DIR *dir = opendir( mp.c_str() );
	if ( !dir )
		{
		cerr << "Could not find path to directory" << endl;
		return false;
		}

	/* read the directory's contents */
	struct dirent *dirEntry;
	bool foundFile = false;
	while ( (dirEntry = readdir(dir)) != NULL ) 
		{
		string temp = dirEntry->d_name;
		if ( temp == mf )
			foundFile = true;
		}

	/* close the directory */
	if ( closedir(dir) == -1 )
		{
		cerr << "Problem closing directory" << endl;
		return false;
		}

	return foundFile;
	
}



bool IoManager::listDirectoryContents(void) {

	/* open the directory */
	DIR *dir = opendir( filePath.c_str() );
	if ( !dir )
		{
		cerr << "Could not find path to directory" << endl;
		return false;
		}

	/* read the directory's contents */
	struct dirent *dirEntry;
	while ( (dirEntry = readdir(dir)) != NULL ) 
		{
		cout << dirEntry->d_name << endl;
		}

	/* close the directory */
	if ( closedir(dir) == -1 )
		{
		cerr << "Problem closing directory" << endl;
		return false;
		}

	return true;
	
}



bool IoManager::openFile(ifstream &strm) {
	
	string delimiter = PATH_SEPERATOR;

	/* concactenate path and file name */
	string filePathName = filePath + delimiter + fileName;

	/* here we assume that the presence of the path/file has
	   been checked elsewhere */
	strm.open( filePathName.c_str(), ios::in );
	if ( !strm )
		return false;
	return true;
	
}



bool IoManager::openFile(ofstream &strm) {
	
	string delimiter = PATH_SEPERATOR;

	/* concactenate path and file name */
	string filePathName = filePath + delimiter + fileName;

	/* here we assume that the presence of the path/file has
	   been checked elsewhere */
	strm.open( filePathName.c_str(), ios::out );
	if ( !strm )
		return false;
	return true;
	
}



void IoManager::closeFile(ifstream &strm) {

	strm.close();
	
}



bool IoManager::testDirectory(void) {

	return isDirectoryPresent(filePath);

}



bool IoManager::testFile(void) {

	return isFilePresent(filePath, fileName);

}







