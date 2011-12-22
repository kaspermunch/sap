#ifndef INTERFACE_H
#define INTERFACE_H

#include "cConstrainedNJlib.h"

extern "C" 
{

  //void initialize(int dim, double **cache) {
void initialize(int dim) {
  initCache(dim);
}

void cleanup(int dim) {
  deleteCache(dim);
}

char *computeTree(int a_nrOTUs, char **a_alignment, int a_nrBackboneSets, char **a_backboneSetsList, int a_resample, int a_branchlength) {
  return compute(a_nrOTUs, a_alignment, a_nrBackboneSets, a_backboneSetsList, a_resample, a_branchlength);
}

}
#endif
