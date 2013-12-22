#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include <iostream>
#include <algorithm>
#include <cstring>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <limits>
#include <exception>

#include "cConstrainedNJlib.h"
#include "utils.h"
#include "interface.h"


extern "C" void initcConstrainedNJlib (void) { /* should be init + name of module as given in setup.py */
  char *dummy = NULL; // dummy body
}


/** Global variables: ******************************************************************************/

// Maximal number of allowed sequences:
#define MAX_OTUS 1000

// Maximal allowed sequence length:
#define MAX_ALIGN_LENGTH 2000

// Array of array alignment[<sequence_nr>][base_nr]
char **alignment;

// Array of arrays each specifying a constraint.  (name integers)
std::vector<std::set<int> > backboneSetsList;

// Array of tuples specifying allowed pairs.     (index integers)
std::list<int *> allowedPairsList;

// Array of otus not part of back bone constraint.
std::vector<int> unconstrainedOTUs;

int N, K, L;

// Sequence alignment length:
int alignmentLength;

// Number of sequences (we need this for the other dimension of the alignmnet)
int nrOTUs;

// list of length alignmentLength for indexing into the alignment
int *alignmentIndexList;

// Minimum number of aligned bases allowed. 
int minPairAlnBases = 20;

// List of nodes specifying the growing tree:
Node *nodeList;

// index list to mapping 0 -> L integers to indexes in the matrices.
std::vector<int> idxL;

// index list to mapping 0 -> L integers to indexes in the node list.
std::vector<int> idxN;

// matrix of d values
// double Qmatrix[MAX_OTUS][MAX_OTUS];
double **Qmatrix;

// matrix of distances  - replaces getDist
// double distMatrix[MAX_OTUS][MAX_OTUS];
double **distMatrix;

// Vector for holding a list of r values:
double *rSums;

// Boolean indicating if resamplingMatrix has been built
bool resamplingMatrixComputed = false;

// A matrix of [transition, transversion, length] lists we then can sample distances from
int ***resamplingMatrix;

// For cached distances:
//double k2pCacheMatrix[MAX_OTUS][MAX_OTUS];
double **k2pCacheMatrix;

/** A node class for the tree: ** (see .h file for decl of class) **********************************/

Node::Node() {
  left = NULL;
  right = NULL;
}

// Node::~Node() {
//   free(left);
//   free(right);     /* right children, or None */
//   free(name);
//   free(distLeft);  /* Length of edges to the children, if any */
//   free(distRight);
//   std::set<int> leafSet;     /* list of subnodes as set in each node */
// }

// The assign and get methods here we just keep in case we need to use the Node Class for
// something else later:
void Node::assignLeft (Node node) {
  left = &node;
}

void Node::assignRight (Node node) {
  right = &node;
}

Node Node::getLeft (void) {
  return *left;
}

Node Node::getRight (void) {
  return *right;
}

char *Node::getSubTreeString(char *buffer) {

  if (left != NULL && right != NULL) {
    //char leftBuffer[1000000];
    char *leftBuffer = new char[10000];
    //char rightBuffer[1000000];
    char *rightBuffer = new char[10000];
    if ((*left).leafSet.size() > (*right).leafSet.size()) {
      if (branchLengths) {
	sprintf(buffer, "(%s:%.5f,%s:%.5f)", (*left).getSubTreeString(leftBuffer), distLeft, (*right).getSubTreeString(rightBuffer), distRight);
      } else {
	sprintf(buffer, "(%s,%s)", (*left).getSubTreeString(rightBuffer), (*right).getSubTreeString(leftBuffer));
      }
    } else {
      if (branchLengths) {
	sprintf(buffer, "(%s:%.5f,%s:%.5f)", (*right).getSubTreeString(rightBuffer), distRight, (*left).getSubTreeString(leftBuffer), distLeft);
      } else {
	sprintf(buffer, "(%s,%s)", (*right).getSubTreeString(rightBuffer), (*left).getSubTreeString(leftBuffer));
      }
    }
    delete [] leftBuffer;
    delete [] rightBuffer;
  } else {
    sprintf(buffer,"%d",name);
  }
  return buffer;
}

void freeExternMemmory(void) {
  
  extern int N;
  extern Node *nodeList;
  extern int *alignmentIndexList;
  extern std::list<int *> allowedPairsList;
  extern double **distMatrix;
  extern double **Qmatrix;
    
  std::list<int *>::iterator it;
  for ( it = allowedPairsList.begin() ; it != allowedPairsList.end() ; ++it )
    delete [] *it;
  allowedPairsList.clear();

  for(int i=0; i<N;i++) {
    delete [] distMatrix[i]; 
    delete [] Qmatrix[i]; 
  }
  delete [] distMatrix;
  delete [] Qmatrix;
  delete [] rSums;
  delete [] alignmentIndexList;
  delete [] nodeList;
}

/** Calculating K2P distances or getting them from a cache: ****************************************/

K2Pstats k2pStats(int i, int j) {
  // Calculates the statistics required to calculate the k2p distance

  extern int alignmentLength;
  extern char **alignment;
  extern int *alignmentIndexList;

  K2Pstats k2p;
  char x, y;
  int extraGaps = 0;
  int indels = 0;
  int transitions = 0;
  int transversions = 0;
  int uninformative = 0;

  char ambigousWildCards[] = "MKWSBDHVN";
  int nrAmbigousWildCards = 9;

//   for (int b=0; b<alignmentLength; b++)
//     std::cout << alignmentIndexList[b] << ' ';
//   std::cout << std::endl;
    
//     A = adenine (purine), C = cytosine (pyrimidine), G = guanine (purine), T = thymine (pyrimidine), U = uracil
//     R = G A (purine), Y = T C (pyrimidine), K = G T (keto), M = A C (amino)
//     S = G C (strong bonds), W = A T (weak bonds)
//     B = G T C (all but A), D = G A T (all but C), H = A C T (all but G), V = G C A (all but T),
//     N = A G C T (any)

  int idx;
  for (idx=0; idx < alignmentLength; idx++) {

    x = alignment[i][alignmentIndexList[idx]];
    y = alignment[j][alignmentIndexList[idx]];

    // Skip if a char is a gap:
    if (x == '-' && y == '-') {
  	 extraGaps += 1;
  	 continue;
    } else if (x == '-' || y == '-') {
  	 indels += 1;
  	 continue;
    } else if (x == y) {
  	 continue;
    }
  
    // Skip if a char is an ambigous wild card:
    int i;
    for (i=0; i<nrAmbigousWildCards; i++) {
  	 if (x == ambigousWildCards[i]) {
  	   uninformative += 1;
  	   continue;
  	 }
  	 if (y == ambigousWildCards[i]) {
  	   uninformative += 1;
  	   continue;
  	 }
    }

    // Count transitions and transversions:
    switch(x) {
    case 'A':
      switch(y) {
      case 'G': case 'R':
	transitions += 1;
	break;
      case 'T': case 'C': case 'Y':
	transversions += 1;
	break;
      }
      break;
    case 'C': case  'Y':
      switch(y) {
      case 'T':
	transitions += 1;
	break;
      case 'A': case 'G': case 'R':
	transversions += 1;
	break;
      }
      break;
    case 'T':
      switch(y) {
      case 'C': case 'Y':
	transitions += 1;
	break;
      case 'A': case 'G': case 'R':
	transversions += 1;
	break;
      }
      break;
    case 'G':
      switch(y) {
      case 'A': case 'R':
	transitions += 1;
	break;
      case 'T': case 'C': case 'Y':
	transversions += 1;
	break;
      }
      break;
    case 'R':
      switch(y) {
      case 'Y':
	transversions += 1;
	break;
      }
      break;
    }  
  }  

  int alignedBases = alignmentLength - extraGaps - indels - uninformative;
  
  if (alignedBases < minPairAlnBases) {
//     freeExternMemmory();
    throw TooFewAlnBasesError;
  }

  k2p.transitions = transitions;
  k2p.transversions = transversions;
  k2p.alignedBases = alignedBases;

  return k2p;
}

double k2pDist(K2Pstats k2p) {
  // Calculates the k2p distance from the summary statistics

  extern int alignmentLength;

  double Q, P, distance;

//   if (k2p.alignedBases == alignmentLength && k2pCacheMatrix[k2p.transitions][k2p.transversions] != -1) {
//     distance = k2pCacheMatrix[k2p.transitions][k2p.transversions];
//   } else {
    P = (double) k2p.transitions / k2p.alignedBases;
    Q = (double) k2p.transversions / k2p.alignedBases;    
    if (2 * P - Q >= 1 ||  2 * Q >= 1) {
      distance = 10000;
    } else {
      distance = 0.5 * log(1 / (1 - 2 * P - Q)) + 0.25 * log(1 / (1 - 2 * Q));
    }
//   }


 if (std::numeric_limits<double>::has_quiet_NaN && distance == std::numeric_limits<double>::quiet_NaN())
     distance = 10000;	

 if (std::numeric_limits<double>::has_signaling_NaN && distance == std::numeric_limits<double>::signaling_NaN())
     distance = 10000;	

 if ( !(distance == distance) ) {
     distance = 10000;
 }

  // Hack to avoid craches because the distance becomes inf:
  if ( distance == std::numeric_limits<double>::infinity() )
    distance = 10000;

//   k2pCacheMatrix[k2p.transitions][k2p.transversions] = distance;

  return distance;
}

double computeDistance(int i, int j) {
  // Calculates the k2p distance between i and j

  K2Pstats k2p = k2pStats(i, j);

  return k2pDist(k2p);
}

double resampleDistance(int i, int j) {

  std::cout << "RESAMPLING IS NOT IMPLEMENTED" << std::endl;
  exit(1);

//   extern double **resamplingMatrix;
// 
//   K2Pstats k2p = k2pStats(i, j);
// 
// //     transitions, transversions, alignedBases = resamplingMatrix[i][j]
// // 
// //     transitions = int(round(random.randint(1, alignedBases) * (float(transitions) / alignedBases)))
// //     transversions = int(round(random.randint(1, alignedBases) * (float(transversions) / alignedBases)))
// // 
// //     return k2pDist(transitions, transversions, alignedBases)
// 
//   return k2pDist(k2p);
}

/** Maintaining the Q matrix: ********************************************************************/

double calcRsum(int i) {
  extern double **distMatrix;
  extern std::vector<int> idxL;
  extern int L;

  double r = 0;
  for (int j=0; j<L; j++) {
    r += distMatrix[i][idxL[j]];
  }
  return r / (L - 2);
}

void updateMatrices(int old_i, int old_j) {

  extern double **distMatrix;
  extern double **Qmatrix;
  extern std::vector<int> idxL;
  extern int L, N;
  extern double *rSums;

  // Remove the node with the largest integer name from the list of nodes to be joined:
  for (int k=0; k<L; k++) {
    if (idxL[k] == old_j) {
      idxL.erase(idxL.begin() + k);
      break;
    }
  }
  // Decrement the number of nodes to be joined correspondingly:
  L -= 1;

  // Update the column old_i in the distance matrix that now repesents the joined pair:
  for (int k=0; k<L; k++)
    distMatrix[idxL[k]][old_i] = ( -distMatrix[old_i][old_j] + distMatrix[idxL[k]][old_j] + distMatrix[idxL[k]][old_i] ) / 2.0;
  // Updata the corresponding row:
  for (int k=0; k<L; k++)
    distMatrix[old_i][idxL[k]] = distMatrix[idxL[k]][old_i];
  // If there are only two nodes left we don't need to update the Qmatrix:
  if (L == 2)
    return;

  // Calculate r sums:
  for (int k=0; k<N; k++)
    rSums[k] = 0;
  // Recalculate relevant entries in Qmatrix:
  std::list<int *>::iterator it;
  for ( it=allowedPairsList.begin() ; it != allowedPairsList.end(); it++ ) {
    int k = **it;
    int l = *(*it+1);
    if (rSums[k] == 0)
      rSums[k] = calcRsum(k);
    if (rSums[l] == 0)
      rSums[l] = calcRsum(l);
    Qmatrix[k][l] = distMatrix[k][l] - (rSums[k] + rSums[l]);
    Qmatrix[l][k] = Qmatrix[k][l];
  }
}

/** Calculating and maintaining the list of pairs allowed to join: *********************************/

bool updateConstraints(int i, int j) {

  extern std::vector<std::set<int> > backboneSetsList;

  bool removedConstraint = false;

  int nrBackboneSets = backboneSetsList.size();
  for (int k=nrBackboneSets-1; k>=0; k--) {
    if (set_member(j+1, backboneSetsList[k])) {
      backboneSetsList[k].erase(j+1);
      if (backboneSetsList[k].size() == 1) {
	removedConstraint = true;
	backboneSetsList.erase(backboneSetsList.begin() + k);
      }
    }
  }

  return removedConstraint;
} 

void computeAllowedPairs(void) {

  extern int N;
  extern std::vector<std::set<int> > backboneSetsList;
  extern std::list<int *> allowedPairsList;

  std::list<int *>::iterator it;
  for ( it = allowedPairsList.begin() ; it != allowedPairsList.end() ; ++it )
    delete [] *it;
  allowedPairsList.clear();

  for (int i=0; i<N; i++) {
    int inConstraint = 0;
    std::set<int> constrainedSet;
    
    int nrBackboneSets = backboneSetsList.size();
    for (int j=0; j<nrBackboneSets; j++ ) {
      std::set<int> s = backboneSetsList[j];
      if (set_member(i+1, s)) {
	inConstraint = 1;
	std::vector<int> nodeIndeces;
	std::set<int>::iterator sit;
	for (sit=s.begin() ; sit != s.end(); sit++ ) {
	  if (!set_member(*sit, constrainedSet))
	    nodeIndeces.push_back(*sit - 1);
	}
	std::vector<int *> pairs = get_pairwise_combinations(nodeIndeces);
	int nrPairs = pairs.size();
	for (int k=0; k<nrPairs; k++) {
	  allowedPairsList.push_back(pairs[k]);
	}
	break;
      } else {
	constrainedSet = set_union(constrainedSet, s);	
      }
    }
    if (!inConstraint) {
      for (int j=0; j<N; j++) {
	if (i != j) {
	  //int *pair = (int *) calloc(2, sizeof(int));
	  int *pair = new int[2];
	  pair[0] = intMin(i, j);
	  pair[1] = intMax(i, j);
	  allowedPairsList.push_back(pair);
	}
      }      
    }
    allowedPairsList.sort(pair_sort);
    allowedPairsList.unique(pair_identity);
  }
}

void updateAllowedPairs(int i, int j) {

  extern std::list<int *> allowedPairsList;
  extern std::vector<std::set<int> > backboneSetsList;
  extern std::vector<int> unconstrainedOTUs;

  std::list<int *> newList;

//   if (findInIntVector(i, unconstrainedOTUs)) {
//     // If i is not in the backbone sets we don't want i to represent the child
//     // nodes. If they are both unconstrained swopping is not needed but does not matter.      
//     swapInts(&i, &j);
//   }

  bool removedConstraint = updateConstraints(i, j);

  std::list<int *>::iterator it;
  for ( it=allowedPairsList.begin() ; it != allowedPairsList.end(); it++ ) {
        
    if ((*it)[0] == j || (*it)[1] == j) {

      if (findInIntVector(j, unconstrainedOTUs))
	continue;
      
      if ((*it)[0] == i && (*it)[1] == j || 
	  (*it)[0] == j && (*it)[1] == i) { // we may have swopped them.
	continue;
      } else if ((*it)[0] == j) {
	//int *p = (int *) calloc(2, sizeof(int));
	int *p = new int[2];
	p[0] = intMin(i, (*it)[1]);
	p[1] = intMax(i, (*it)[1]);
	newList.push_back(p);
      } else if ((*it)[1] == j) {
	//int *p = (int *) calloc(2, sizeof(int));
	int *p = new int[2];
	p[0] = intMin(i, (*it)[0]);
	p[1] = intMax(i, (*it)[0]);
	newList.push_back(p);
      } else {
	std::cout << "WARNING: Not supposed to happen" << std::endl;
	exit(1);
      }
    } else {
      int *p = new int[2];
      p[0] = (*it)[0];
      p[1] = (*it)[1];
      newList.push_back(p);
      //newList.push_back((*it));
    }
  }

  if (removedConstraint) {
    // if we removed a set it means that i now represents all otus
    // from that set. So we find the first set in the updated
    // constraints list that includes i and then add all
    // combinations of i to items in this set:
    std::set<int> constrainedSet;

    int nrBackboneSets = backboneSetsList.size();
    for (int k=0; k<nrBackboneSets; k++ ) {
	std::set<int> s = backboneSetsList[k];
	if (set_member(i+1, s)) {

	  std::vector<int> nodeIndeces;
	  std::set<int>::iterator sit;
	  for (sit=s.begin() ; sit != s.end(); sit++ ) {
	    if (!set_member(*sit, constrainedSet))
	      nodeIndeces.push_back(*sit - 1);
	  }

	  std::vector<int *> pairs = get_pairwise_combinations(nodeIndeces);
	  std::vector<int *>::iterator lit;
	  for ( lit=pairs.begin() ; lit != pairs.end(); lit++ ) {
	    newList.push_back(*lit);
	  }
	  break;
	} else {
	  constrainedSet = set_union(constrainedSet, s);
	}
    }  
  }
  newList.sort(pair_sort);

  newList.unique(pair_identity);  
  for ( std::list<int *>::iterator it=allowedPairsList.begin() ; it != allowedPairsList.end(); it++ )
    delete *it;
  allowedPairsList.clear();
  allowedPairsList = newList;

//   newList.unique(pair_identity);  
//   allowedPairsList = newList;
}

/** Finding a pair to join and create a new node: **************************************************/

int *findPair(void) {

  extern double **Qmatrix;
  extern std::list<int *> allowedPairsList;

//   int mini = -1;
//   int minj = -1;
  double mind = 1000000000;
  std::vector<int *> pairList;

  std::list<int *>::iterator it;
  for ( it=allowedPairsList.begin() ; it != allowedPairsList.end(); it++ ) {
    int i = **it;
    int j = *(*it+1);

    if (Qmatrix[i][j] < mind) {
      mind = Qmatrix[i][j];
      std::vector<int *>::iterator pit;
      for ( pit = pairList.begin() ; pit != pairList.end() ; ++pit )
	delete [] *pit;
      pairList.clear();      
      //int *p = (int *) calloc(2, sizeof(int));
      int *p = new int[2];
      p[0] = intMin(i, j);
      p[1] = intMax(i, j);
      pairList.push_back(p);
    } else if (Qmatrix[i][j] == mind) {
      //int *p = (int *) calloc(2, sizeof(int));
      int *p = new int[2];
      p[0] = intMin(i, j);
      p[1] = intMax(i, j);
      pairList.push_back(p);
    }
  }
  if (pairList.size() == 0) {
    std::cout << "WARNING: no entries in pairlist" << std::endl;
    exit(1);
  }

  // Free the memory for the pairs we don't return:
  int pidx = rand() % pairList.size();
  int pairListSize = pairList.size();
  for ( int i=0; i<pairListSize; i++ ) {
    if (i != pidx)
      delete [] pairList[i];
  }
  return pairList[pidx];
}

void createParentNode(int i, int j) {

  extern double **distMatrix;
  extern Node *nodeList;
  extern int K;
  extern std::vector<int> idxN;

  double dik = (distMatrix[i][j] + calcRsum(i) - calcRsum(j)) / 2;
  double djk = distMatrix[i][j] - dik;

  Node node;  
  node.name = K+1;
  node.left = &nodeList[idxN[i]];
  node.distLeft = dik;
  node.right = &nodeList[idxN[j]];
  node.distRight = djk;
  node.leafSet = set_union(nodeList[idxN[i]].leafSet, nodeList[idxN[j]].leafSet);
  node.branchLengths = nodeList[idxN[i]].branchLengths;

  nodeList[K] = node;

  idxN[i] = K;
}

/** The main functions: ****************************************************************************/

//void initCache(int dim, double **cache) {
void initCache(int dim) {

  extern double **k2pCacheMatrix;

  //  k2pCacheMatrix = cache;
  k2pCacheMatrix = new double *[dim];
  for (int i=0; i<dim; i++) {
    k2pCacheMatrix[i] = new double[dim];
    for (int j=0; j<dim; j++) {
      k2pCacheMatrix[i][j] = -1;
    }
  }
}

void deleteCache(int dim) {
  extern double **k2pCacheMatrix;

  for(int i=0; i<dim;i++) {
    delete [] k2pCacheMatrix[i]; 
  }
  delete [] k2pCacheMatrix;
}

char *compute(int a_nrOTUs, char **a_alignment, int a_nrBackboneSets, char **a_backboneSetsList, int a_resample, int a_branchlengths) {

  // Called like this from Python: computeTree(alignment, constraintList, queryName)

//   // REMEMBER TO REMOVE THIS AGAIN.....
//   srand(7);
//   // REMEMBER TO REMOVE THIS AGAIN.....

  extern int nrOTUs, N, L, K, alignmentLength;
  extern char **alignment;
  extern Node *nodeList;
  extern int *alignmentIndexList;
  extern std::vector<std::set<int> > backboneSetsList;
  extern std::list<int *> allowedPairsList;
  extern std::vector<int> unconstrainedOTUs;
  extern double **distMatrix;
  extern double **Qmatrix;

  int i, j;

  int resample = a_resample;
  int branchlengths = a_branchlengths;

  // Populate teh global variables:
  nrOTUs = a_nrOTUs;
  alignment = a_alignment;

  // Sequence length of the alignment:
  alignmentLength = 0;
  for (i=0; alignment[0][i] != '\0'; i++)
    alignmentLength += 1;

  // The initial number of leaves
  N = nrOTUs;
  // The number of clusters created so far (including the first N "leaf culsters" we start with)
  K = N;
  // The number of nodes left to be joined:
  L = N;

  // Allocate for the list of nodes:
  nodeList = new Node[2*N-1];

  // Turn lists specifying unconstrained sets in to set objects:
  std::string s;
  std::string delim = " ";
  std::vector<std::string> splitList;
  int otu;
  
  backboneSetsList.clear();

  for (i=0; i<a_nrBackboneSets; i++) {
    s.assign(a_backboneSetsList[i]);
    splitList = stringSplit(s, delim);
    std::set<int> constSet;
    for (unsigned int j=0; j<splitList.size(); j++) {
      otu = atoi(splitList[j].c_str());
      constSet.insert(otu);
    }
    backboneSetsList.push_back(constSet);
  }

  std::vector<int>::iterator it;
  idxL.clear();
  idxN.clear();
  for (i=0; i<N; i++) {
    // Initialize the vector that keeps track of indeces in the distMatrix:
    idxL.push_back(i);
    // Initialize the vector that keeps track of indeces in the nodeList:
    idxN.push_back(i);
  }

  // Find the otus that are not part of the backbone constraint:
  std::set<int> unionSet;
  for (unsigned int i=0; i<backboneSetsList.size(); i++) {
    unionSet = set_union(unionSet, backboneSetsList[i]);
  }

  unconstrainedOTUs.clear();
  for (i=0; i<N; i++) {
    if (!set_member(i+1, unionSet))
      unconstrainedOTUs.push_back(i);
  }

  // The average distance to other leaves:
  rSums = new double[N];

  // Initilize the index list:
  alignmentIndexList = new int[alignmentLength];
  for (int i=0; i<alignmentLength; i++) {
    alignmentIndexList[i] = i;
  }

  // Create the distMatrix:
  distMatrix = new double *[N];
  for (int i=0; i<N; i++)
    distMatrix[i] = new double[N];

  int *tmpList = new int[alignmentLength];
  if (resample == 1) {

    // Resample the index list:
    for (int i=0; i<alignmentLength; i++) {
      j = rand() % alignmentLength;
      tmpList[i] = alignmentIndexList[j];
    }
    delete [] alignmentIndexList;
    alignmentIndexList = tmpList;

    // Fill in distMatrix
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
	if (i < j) {
	  distMatrix[i][j] = computeDistance(i, j);
	} else if (i > j) {
	  distMatrix[i][j] = distMatrix[j][i];
	} else {
	  distMatrix[i][j] = 0.0;
	}
      }
    }
//
//     if not resamplingMatrixComputed:
//         # Create resamplingMatrix:
//         for i in range(N):
//             resamplingMatrix.append([])
//             for j in range(N):
//                 resamplingMatrix[i].append(None)
//         # Fill in the resampling matrix if this has not been done already:
//         for i in range(N):
//             for j in range(N):
//                 if i < j:
//                     resamplingMatrix[i][j] = k2pStats(i, j)
//                 elif i > j:
//                     resamplingMatrix[i][j] = resamplingMatrix[j][i]
//                 else:
//                     resamplingMatrix[i][j] = [0, 0, 1]
//         resamplingMatrixComputed = true;
//                 
//     # Fill in distMatrix
//     for i in range(N):
//         for j in range(N):
//             if i < j:
//                 distMatrix[i][j] = resampleDistance(i, j)
//             elif i > j:
//                 distMatrix[i][j] = distMatrix[j][i]
//             else:
//                 distMatrix[i][j] = 0.0
  } else {
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
	if (i < j) {
	  distMatrix[i][j] = computeDistance(i, j);
	} else if (i > j) {
	  distMatrix[i][j] = distMatrix[j][i];
	} else {
	  distMatrix[i][j] = 0.0;
	}
      }
    }
  }

  // Fill in Create Qmatrix:
  Qmatrix =  new double *[N];
  for (int i=0; i<N; i++) {
    Qmatrix[i] = new double[N];
    for (int j=0; j<N; j++) {
      if (i <= j) {
	Qmatrix[i][j] = distMatrix[i][j] - (calcRsum(i) + calcRsum(j));
      } else {
	Qmatrix[i][j] = Qmatrix[j][i];
      }
    }
  }


  //printDistMatrix();
  //printQmatrix();

  // Create initial N leaf nodes:
  for (int i=0; i<N; i++) {
   Node node;  
   node.name = i+1;
   node.branchLengths = branchlengths;
   std::set<int> s;
   s.insert(i+1);
   node.leafSet = s;
   nodeList[i] = node;
  }

  // Compute the pairs that are allowed to join at first:
  computeAllowedPairs();

  // Main loop:
  while (K < 2 * N - 2) {
    int *p = findPair();

    if (findInIntVector(p[0], unconstrainedOTUs)) {
      // If i is not in the backbone sets we don't want i to represent the child
      // nodes. If they are both unconstrained swopping is not needed but does not matter.      
      swapInts(&p[0], &p[1]);
    }

    createParentNode(p[0], p[1]);
    updateAllowedPairs(p[0], p[1]);
    updateMatrices(p[0], p[1]);
    K += 1;
  }
  // only one join remains and this should be reflected in allowedPairsList..
  if (allowedPairsList.size() != 1) {
    std::cout << "WARNING: allowedPairsList is not one - it is " << allowedPairsList.size() << std::endl;
    exit(1);
  }
  int *p = allowedPairsList.front();
  double dij = distMatrix[p[0]][p[1]] / 2;
  // (We devide by two not to count this distance twice because we arbitrarily insert a root node.)

  Node node;  
  node.name = K+1;
  node.left = &nodeList[idxN[p[1]]];
  node.distLeft = dij;
  node.right = &nodeList[idxN[p[0]]];
  node.distRight = dij;
  node.leafSet = set_union(nodeList[idxN[p[1]]].leafSet, nodeList[idxN[p[0]]].leafSet);
  nodeList[K] = node;

  // Make a string representation:
  char *treeString = new char[10000];
  //char treeString[1000000];
  treeString = nodeList[K].getSubTreeString(treeString);
  //char *treeString = nodeList[K].getSubTreeString();

  // Free up memory:
  freeExternMemmory();
  
  return treeString;
  //return nodeList[K].getSubTreeString();
}

/***************************************************************************************************/

int main(void) {

//    int a_nrOTUs = 102;
// 
// 	char *a_alignment[102] = {		
// 		"TGAGCAGGTATAGTGGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCAGGATCTTTAATTGGAGATGATCAAATTTATAACACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCCGGGTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTGACAGCTCACGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCATTAAGACTTTTAATTCGTACAGAATTAGGAAACCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGACCAAATTTATAATACTATTGTCACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGTACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTGATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCACGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGGTCTTTAATTGGAGATGACCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATGGTAGGAACTTCATTAAGACTTTTAATTCGAACTGAATTAGGAAACCCAGGATCCTTAATTGGAGATGATCAAATTTATAATACTATTGTGACAGCTCATGCATTCATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTCTAATTCGTACTGAATTAGGAAATCCAGGATCTTTAATTGGAGACGATCAAATTTATAATACTATTGTTACAGCTCATGCCTTTATTATAATTTTTTT",
// 		"TGATCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCAGGATCATTAATTGGAGATGATCAAATTTATAATATTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGTACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGCATAGTAGGAACTTCTCTCAGACTTTTAATCCGAACTGAATTAGGTAACCCCGGATCTTTAATTGGTGATGATCAAATTTATAACACTATTGTCACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGTCTTTTAATTCGAACAGAATTAGGAAACCCTGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACATCATTAAGTTTATTAATTCGAACTGAATTAGGAAACCCAGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTCTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGTTTATTAATTCGTGCTGAATTAGGTAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATATTAGGAACTTCTTTAAGTCTTTTAATTCGAACTGAATTAGGAAATCCAGGATCTTTAATTGGTGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACATCATTAAGACTTTTAATTCGAGCTGAATTAGGAAATCCTGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATCGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCATTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTCTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATCGGGGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCCCTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCTGGTATAGTAGGAACTTCCTTAAGACTTTTAATTCGAACTGAATTAGGTAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCATTAAGATTATTAATTCGAGCTGAATTAGGTAACCCCGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACATCTCTTAGTTTATTAATTCGAACAGAATTAGGAAATCCTGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCATTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTTGGAACATCACTAAGACTTCTAATTCGAACTGAATTAGGAAACCCAGGATCTTTAATTGGTGATGATCAAATTTACAATACCATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGTCTTTTAATTCGAACAGAATTAGGAAATCCAGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTTGGAACATCTTTAAGACTTTTAATTCGAACTGAATTAGGAAACCCAGGATCATTAATTGGTGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATCTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGTACTTCTTTAAGTCTTCTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCCCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACCTCATTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCCTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGTGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTTGGAACCTCTTTAAGTCTTTTAATTCGAACTGAATTAGGTAATCCAGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACTGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCATTCATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCACTTAGTTTATTAATTCGAACTGAATTAGGTAATCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACAATTGTCACAGCACATGCTTTTATTATAATTTTTTT",
// 		"TGATCAGGAATAGTAGGAACATCTTTAAGATTATTAATTCGAACAGAATTAGGAAATCCTGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTCACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGATCAGGAATAGTAGGAACATCTTTAAGATTATTAATTCGAACAGAATTAGGAAATCCTGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTCACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACATCTTTAAGTCTTTTAATTCGAACTGAATTAGGAACCCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCCCTTAGTCTTATTATCCGAACAGAATTAGGAAATCCAGGTTCATTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCCCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATATTAGGAACTTCTTTAAGAATTTTAATTCGAATAGAATTAGGAACCCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACAATTGTAACAGCTCATGCATTTATTATAATTTTCTT",
// 		"TGAGCAGGAATACTAGGAACTTCTTTAAGAATTCTTATTCGAATAGAATTAGGAACTCCGGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTGGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCAGGATCTTTAATTGGAGATGATCAAATTTATAACACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCCGGGTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTGACAGCTCACGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCATTAAGACTTTTAATTCGTACAGAATTAGGAAACCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGACCAAATTTATAATACTATTGTCACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGTACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTGATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCACGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGGTCTTTAATTGGAGATGACCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATGGTAGGAACTTCATTAAGACTTTTAATTCGAACTGAATTAGGAAACCCAGGATCCTTAATTGGAGATGATCAAATTTATAATACTATTGTGACAGCTCATGCATTCATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTCTAATTCGTACTGAATTAGGAAATCCAGGATCTTTAATTGGAGACGATCAAATTTATAATACTATTGTTACAGCTCATGCCTTTATTATAATTTTTTT",
// 		"TGATCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCAGGATCATTAATTGGAGATGATCAAATTTATAATATTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGTACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGCATAGTAGGAACTTCTCTCAGACTTTTAATCCGAACTGAATTAGGTAACCCCGGATCTTTAATTGGTGATGATCAAATTTATAACACTATTGTCACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGTCTTTTAATTCGAACAGAATTAGGAAACCCTGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACATCATTAAGTTTATTAATTCGAACTGAATTAGGAAACCCAGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTCTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGTTTATTAATTCGTGCTGAATTAGGTAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATATTAGGAACTTCTTTAAGTCTTTTAATTCGAACTGAATTAGGAAATCCAGGATCTTTAATTGGTGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACATCATTAAGACTTTTAATTCGAGCTGAATTAGGAAATCCTGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATCGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCATTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTCTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATCGGGGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCCCTAAGACTTTTAATTCGAACTGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCTGGTATAGTAGGAACTTCCTTAAGACTTTTAATTCGAACTGAATTAGGTAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCATTAAGATTATTAATTCGAGCTGAATTAGGTAACCCCGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACATCTCTTAGTTTATTAATTCGAACAGAATTAGGAAATCCTGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCATTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTTGGAACATCACTAAGACTTCTAATTCGAACTGAATTAGGAAACCCAGGATCTTTAATTGGTGATGATCAAATTTACAATACCATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGTCTTTTAATTCGAACAGAATTAGGAAATCCAGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTTGGAACATCTTTAAGACTTTTAATTCGAACTGAATTAGGAAACCCAGGATCATTAATTGGTGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATCTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGTACTTCTTTAAGTCTTCTAATTCGAACAGAATTAGGAAACCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCCCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACCTCATTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCCTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCTGGATCTTTAATTGGTGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTTGGAACCTCTTTAAGTCTTTTAATTCGAACTGAATTAGGTAATCCAGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACTGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACAGAATTAGGAAATCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCATTCATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCACTTAGTTTATTAATTCGAACTGAATTAGGTAATCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACAATTGTCACAGCACATGCTTTTATTATAATTTTTTT",
// 		"TGATCAGGAATAGTAGGAACATCTTTAAGATTATTAATTCGAACAGAATTAGGAAATCCTGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTCACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGATCAGGAATAGTAGGAACATCTTTAAGATTATTAATTCGAACAGAATTAGGAAATCCTGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTCACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACATCTTTAAGTCTTTTAATTCGAACTGAATTAGGAACCCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGTATAGTAGGAACTTCCCTTAGTCTTATTATCCGAACAGAATTAGGAAATCCAGGTTCATTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCCCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATATTAGGAACTTCTTTAAGAATTTTAATTCGAATAGAATTAGGAACCCCAGGATCTTTAATTGGAGATGATCAAATTTATAATACAATTGTAACAGCTCATGCATTTATTATAATTTTCTT",
// 		"TGAGCAGGAATACTAGGAACTTCTTTAAGAATTCTTATTCGAATAGAATTAGGAACTCCGGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT",
// 		"TGAGCAGGAATAGTAGGAACTTCTTTAAGACTTTTAATTCGAACTGAATTAGGAAATCCCGGATCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTT"
// 	};
// 
// 
// 	 int a_nrBackboneSets = 1;
// 	 
// 	 char *a_backboneSetsList[1] = {
// 		 "40 8"
// 	};
// 
// 
// 
// //   int a_nrOTUs = 5;
// // 
// //   char *a_alignment[5] = {
// //     "CTATTGTACGCACAGTAGTACGACAGTATTGTACGCACAGTAGTACGACAGT", 
// //     "GTCTTGTACGCACAGTAGTACGACAGTATTGTACGCACAGTAGTACGACAGT", 
// //     "GCTTTGTACGCACAGTAGTACGACAGTATTGTACGCACAGTAGTACGACAGT", 
// //     "CGGTTGTACGCACAGTAGTACGACAGTATTGTACGCACAGTAGTACGACAGT", 
// //     "TAATTGTACGCACAGTAGTACGACAGTATTGTACGCACAGTAGTACGACAGT"
// //   };
// // 
// //   int a_nrBackboneSets = 2;
// // 
// //   char *a_backboneSetsList[2] = {
// //     "1 2 3", 
// //     "1 2 3 4"
// //   };
// 
//   int resample = 0;
// 
// //   for (int i=0; i<a_nrOTUs; i++) {
// //     std::cout << a_alignment[i] << std::endl;
// //   }
// 
//   initCache(200);
// 
//   char *treeStr = NULL;
//   for (int i=0; i<1000; i++) {
//     try
//       {
// 	treeStr = computeTree(a_nrOTUs, a_alignment, a_nrBackboneSets, a_backboneSetsList, resample);
//       }
//     catch ( std::exception& e)
//       {
// 	std::cout << "Exception caught: " << e.what() << std::endl;
//       }
//     
//     std::cout << treeStr << std::endl;
//   }
// 
// //   std::vector<int> elements;
// //   for (int i=0; i<4; i++) {
// //     elements.push_back(i);
// //   }
// //   std::vector<int*> pairList;
// //   pairList = get_pairwise_combinations(elements);
// 
// 
// /* 
//  *   int x, y;
//  *   x = 7;
//  *   y = 2;
//  *   printf("%d %d\n", x, y);
//  *   x = y;
//  *   printf("%d %d\n", x, y);
//  *   x = 4;
//  *   printf("%d %d\n", x, y);
//  */

  return 1;
}


// Cavendish:~/projects/sap/devel/trunk/ext/constrnj$ g++ -O4 -fPIC -c cConstrainedNJlib.cpp -arch i386Cavendish:~/projects/sap/devel/trunk/ext/constrnj$ g++ -O4 -fPIC -c utils.cpp -arch i386
// Cavendish:~/projects/sap/devel/trunk/ext/constrnj$ Cavendish:~/projects/sap/devel/trunk/ext/constrnj$ 
// Cavendish:~/projects/sap/devel/trunk/ext/constrnj$ g++ -dynamiclib -undefined suppress -flat_namespace *.o -o test.dylib -arch i386



