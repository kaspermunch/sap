#ifndef CONDLIKES_H
#define CONDLIKES_H


class Alignment;
class CondLike {

	public:
                	              CondLike(Alignment *al, int ngc);  
								  ~CondLike(void);
						   void   print(void);
						    int   getNumTaxa(void) { return numTaxa; }
							int   getNumChar(void) { return numChar; }
							int   getNumGammaCats(void) { return numGammaCats; }
							int   getNumStates(void) { return numStates; }
						    int   getNumOfPattern(int i) { return numSitesOfPat[i]; }
						 double   *getClPtr(int whichSpace, int whichNode, int whichSite) { return clsPtr[ whichSpace ][ whichNode ][ whichSite ]; }
                	              
	private:
						    int   numTaxa;
							int   numChar;
							int   numGammaCats;
							int   numStates;
						 double   *cls;
						 double   ****clsPtr;
						    int   *numSitesOfPat;
                            
};

#endif
