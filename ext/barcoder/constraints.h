#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <string>
#include <vector>

using namespace std;



class MbBitfield;
class Bipartition {

	public:
                	              Bipartition(int nt);  
								  ~Bipartition(void);
						   void   addToBipartition(int tn);
						   void   addToMask(int tn);
						   void   normalizeBipartition(void);
						    int   getBipartition(int tn);
							int   getMask(int tn);
					 MbBitfield   *getBipartitionPtr(void) { return parts; }
					 MbBitfield   *getMaskPtr(void) { return mask; }
						   void   print(void);

	private:
						    int   numTaxa;
					 MbBitfield   *parts;
					 MbBitfield   *mask;

};



class Alignment;
class Constraints {

	public:
                	              Constraints(string filePath, string fileName, Alignment *al);  
								  ~Constraints(void);
						    int   getNumBipartitions(void) { return biparts.size(); }
					Bipartition   *getBipartitionPtr(int i) { return biparts[i]; }
						   void   print(void);
                	              
	private:
						    int   numTaxa;
					  Alignment   *alignmentPtr;
		  vector<Bipartition *>   biparts;
                            
};

#endif
