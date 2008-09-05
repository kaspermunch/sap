#ifndef UTILS_H
#define UTILS_H

// Here, we overload the << operator so it works on any set (python style)
template <typename TYPE>
inline std::ostream &operator<<(std::ostream &o, std::set<TYPE> s) {
     o << "set([";
     for (typename std::set<TYPE>::iterator it = s.begin(); it != s.end(); ++it) {
       if (it != s.begin() && it != s.end())
	 o << ", ";	  
       o<<*it;
     }
     o << "])";     
     return o;
}

template<class T, class Comp>
bool set_member(T x, std::set<T,Comp>& s)
// Test if int in in set of ints.
{
 return (s.count(x)==1 ? true : false);
}

// std::set<int> set_intersection(std::set<int>, std::set<int>);
std::set<int> set_union(std::set<int>, std::set<int>);
bool set_includes(std::set<int>, std::set<int>);
bool pair_identity (int *, int *);
bool pair_sort (int *, int *);
void printPairsSet(std::set<int *>);
int intMax(int, int);
int intMin(int, int );
void swapInts(int *, int *);
void printDistMatrix(void);
void printQmatrix(void);
bool findInIntVector(int, std::vector<int>);
std::vector<int*> get_pairwise_combinations(std::vector<int>);
std::vector<std::string> stringSplit(std::string str, std::string);


#endif
