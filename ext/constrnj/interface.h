#ifndef INTERFACE_H
#define INTERFACE_H

#include "cConstrainedNJlib.h"

void initialize(int dim) {
  initCache(dim);
}

char *computeTree(int a_nrOTUs, char **a_alignment, int a_nrBackboneSets, char **a_backboneSetsList, int a_resample) {
  return compute(a_nrOTUs, a_alignment, a_nrBackboneSets, a_backboneSetsList, a_resample);
}

#endif
