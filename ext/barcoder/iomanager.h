#ifndef IOMANAGER_H
#define IOMANAGER_H

#include <string>
#include <fstream>

using namespace std;

class IoManager {

	public:
                            IoManager(void);
                            IoManager(string s);
                   string   getCurDirectory(void) { return curDirectory; }
                   string   getFileName(void) { return fileName; }
                   string   getFilePath(void) { return filePath; }
                   string   getFilePathName(void);
                     void   setCurDirectory(string s) { curDirectory = s; }
                     void   setFileName(string s) { fileName = s; }
                     void   setFilePath(string s) { filePath = s; }
                     void   closeFile(ifstream &strm);
                     bool   openFile(ifstream &strm);
                     bool   openFile(ofstream &strm);
                     bool   testDirectory(void);
                     bool   testFile(void);
                     bool   listDirectoryContents(void);
                     bool   parsePathFileNames(string s);

	private:
                   string   curDirectory;	
                   string   fileName;
                   string   filePath;
                   string   findCurrentDirectory(void);
                     bool   isDirectoryPresent(const string mp);
                     bool   isFilePresent(const string mp, const string mf);

};


#endif
