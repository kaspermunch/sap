#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "utils.h"

/** Stuff for working with sets and pairs: *********************************************************/

std::set<int> set_intersection(std::set<int> s1, std::set<int> s2) {
  // Calculates intersection of two sets of ints.
  std::set<int> s_inter;
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(s_inter, s_inter.begin()));
  return s_inter;

}

std::set<int> set_union(std::set<int> s1, std::set<int> s2) {
  // Calculate union of two sets of ints.
  std::set<int> s_union;
  set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(s_union, s_union.begin()));
  return s_union;

}

bool set_includes(std::set<int> s1, std::set<int> s2) {
  // Test if a set of ints is included in an other set of ints.
  bool answer;
  answer = includes(s1.begin(), s1.end(), s2.begin(), s2.end());
  return answer;
}

struct pair_compare {
  bool operator() (const int *first, const int *second) const
  { 
    if (first[0] < second[0]) {
      return true;
    } else if (first[0] == second[0] && first[1] < second[1]) {
      return true;
    } else {
      return false;
    }
  }
};

bool pair_identity (int *first, int *second) {
  // testing identiry of two pairs for use in list unique.
  if (first[0] == second[0] && first[1] == second[1]) {
    return true;
  } else {
    return false;
  }
}

bool pair_sort (int *first, int *second) {
  // Sorting of pairs.
  if (first[0] < second[0]) {
    return true;
  } else if (first[0] == second[0] && first[1] < second[1]) {
    return true;
  } else {
    return false;
  }
}

void printPairsSet(std::set<int *> s) {
  // Prints a set of pair tuples to the stdout: 
  std::set<int *>::iterator it;
  std::cout << '(';
  for (it=s.begin() ; it != s.end(); it++ ) {
    if (it != s.begin())
      std::cout << ", ";
    std::cout << '[' << **it << ',' << *(*it+1) << ']';
  }
  std::cout << ')' << std::endl;
}

/** Utility functions: *****************************************************************************/

int intMax(int a, int b) {

  return a > b ? a : b;
}

int intMin(int a, int b) {

  return a < b ? a : b;
}

bool findInIntVector(int n, std::vector<int> v) {

  int vSize = v.size();
  for (int i=0; i<vSize; i++) {
    if (n == v[i]) {
      return true;
    }
  }
  return false;
}

void swapInts(int *p1, int *p2) {
  int tmp;
  tmp = *p1;
  *p1 = *p2;
  *p2 = tmp;
}

std::vector<int*> get_pairwise_combinations(std::vector<int> elements) {
  std::vector<int*> pairList;
  for (int i=0; i<elements.size(); i++) {
    for (int j=i+1; j<elements.size(); j++) {
      int* pair = new int[2];
      pair[0] = intMin(elements[i], elements[j]);
      pair[1] = intMax(elements[i], elements[j]);
      pairList.push_back(pair);
    }
  }
  return pairList;
}

std::vector<std::string> stringSplit(std::string str, std::string delim) {
  int cutAt;
  std::vector<std::string> results;
  while( (cutAt = str.find_first_of(delim)) != str.npos ) {
    if(cutAt > 0) {
      results.push_back(str.substr(0,cutAt));
    }
    str = str.substr(cutAt+1);
  }
  if(str.length() > 0) {
    results.push_back(str);
  }
  return results;
}

void printQmatrix(void) {

  extern int N;
  extern double **Qmatrix;

  char buffer[1000];

  std::cout << "Q matrix:" << std::endl;
  for (int i=0; i<N; i++) {
    sprintf(buffer, "%-3d", i+1);
    std::cout << buffer;
    for (int j=0; j<N; j++) {
      sprintf(buffer, "%.3f " , Qmatrix[i][j]);
      std::cout << buffer;
    }
    sprintf(buffer, " %d ", i+1);
    std::cout << buffer << std::endl;
  }
  std::cout << "      ";
  for (int i=0; i<N; i++) {
    sprintf(buffer, "%-6d ", i+1);
    std::cout << buffer;
  }
  std::cout << std::endl;
}

void printDistMatrix(void) {

  extern int N;
  extern double **distMatrix;

  char buffer[1000];

  std::cout << "Distance matrix:" << std::endl;
  for (int i=0; i<N; i++) {
    sprintf(buffer, "%-3d", i+1);
    std::cout << buffer;
    for (int j=0; j<N; j++) {
      sprintf(buffer, "%.3f " , distMatrix[i][j]);
      std::cout << buffer;
    }
    sprintf(buffer, " %d ", i+1);
    std::cout << buffer << std::endl;
  }
  std::cout << "      ";
  for (int i=0; i<N; i++) {
    sprintf(buffer, "%-6d ", i+1);
    std::cout << buffer;
  }
  std::cout << std::endl;
}
