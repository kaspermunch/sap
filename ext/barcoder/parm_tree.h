#ifndef PARM_TREE_H
#define PARM_TREE_H

#include "parm.h"
#include <string>
#include <vector>
#include <sstream>

using namespace std;



class Node {

	public:
                	              Node(void);  
								  Node(Node &n);
						   Node   &operator=(Node &n);
						   void   flipActiveCl(void) { (activeCl == 0 ? activeCl = 1 : activeCl = 0); }
						   void   flipActiveTi(void) { (activeTi == 0 ? activeTi = 1 : activeTi = 0); }
						    int   getActiveCl(void) { return activeCl; }
							int   getActiveTi(void) { return activeTi; }
						   bool   getUpdateCl(void) { return updateCl; }
						   bool   getUpdateTi(void) { return updateTi; }
						 double   getV(void) { return v; }
							int   getIndex(void) { return index; }
						 string   getName(void) { return name; }
						   bool   getIsLeaf(void) { return isLeaf; }
						   Node   *getLft(void) { return lft; }
						   Node   *getRht(void) { return rht; }
						   Node   *getSib(void) { return sib; }
						   Node   *getAnc(void) { return anc; }
						   bool   getFlag(void) { return flag; }
						   void   setV(double x) { v = x; }
						   void   setIndex(int x) { index = x; }
						   void   setName(string s) { name = s; }
						   void   setIsLeaf(bool tf) { isLeaf = tf; }
	                       void   setLft(Node *p) { lft = p; }
						   void   setRht(Node *p) { rht = p; }
	                       void   setSib(Node *p) { sib = p; }
	                       void   setAnc(Node *p) { anc = p; }
						   void   setFlag(bool tf) { flag = tf; }
						   void   setActiveCl(int x) { activeCl = x; }
						   void   setActiveTi(int x) { activeTi = x; }
						   void   setUpdateCl(bool tf) { updateCl = tf; }
						   void   setUpdateTi(bool tf) { updateTi = tf; }

	private:
						 double   v;
						    int   index;
						 string   name;
						   bool   isLeaf;
						   Node   *lft;
						   Node   *rht;
						   Node   *sib;
						   Node   *anc;
						   bool   flag;
						    int   activeCl;
							int   activeTi;
						   bool   updateCl;
						   bool   updateTi;
};



class Constraints;
class MbBitfield;
class Tree : public Parm {

	public:
                	              Tree(MbRandom *rn, Model *mp, int nt, bool ir, double lm, vector<string> names);  
                	              Tree(MbRandom *rn, Model *mp, int nt, bool ir, double lm, vector<string> names, Constraints *cs);  
								  Tree(Tree &t);
								  ~Tree(void);
					       Tree   &operator=(Tree &t);
						 double   change(void);
						   void   clone(Tree &t);
						   void   flipAllActiveCls(void);
						   void   flipAllActiveTis(void);
						   Node   *getDownPassNode(int i) { return downPassSequence[i]; }
					     double   getLnPriorProbability(void);
						    int   getNumNodes(void) { return numNodes; }
						 string   getNewick(void);
						 string   getParmName(void) { return parmName; }
						   Node   *getRootPtr(void) { return root; }
						   void   print(void);
						   void   setAllUpdateCls(bool tf);
						   void   setAllUpdateTis(bool tf);
                	              
	private:
	                       void   buildRandomTree(vector<string> names);
	                       void   buildRandomTree(vector<string> names, Constraints *cs);
						 double   changeBrlen(void);
						 double   changeNni(void);
						 double   changeSpr(void);
						   void   convertTree(void);
						    int   dex(Node *p);
						   Node   *findCommonAncestor(Node **nde, int n);
						    int   getDownPassSeq(void);
						    int   getDownPassSeq2(void);
						    int   getNumAbove(Node *p);
						   void   initializePartitions(void);
						   bool   isTreeConsistentWithConstraints(void);
						   void   listNodes(void);
						   Node   *nodeLeftOf(Node *p);
						   void   passDn(Node *p, int *x);
						   void   passDn2(Node *p, int *x);
						   void   printParts(void);
						   void   printParts(MbBitfield *msk);
						   Node   *rightMostDescendant(Node *p);
						   void   rootAtNode(Node *p);
						   void   rotateTree(int *nextInteriorNode);
						   void   showNodes(Node *p, int indent);
						   void   showNodes2(Node *p, int indent);
						   void   writeTree(Node *p, stringstream &ss);
						   bool   isRooted;
						 double   lambda;
						    int   numTaxa;
							int   numNodes;
						   Node   *nodes;
						   Node   **downPassSequence;
						   Node   *root;
					 MbBitfield   **parts;
					Constraints   *constraintsPtr;
					     double   tuningBrlen;
						 double   tuningNni;
						 double   tuningSpr;
};

#endif
