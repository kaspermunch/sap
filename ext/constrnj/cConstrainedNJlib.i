
%typemap(in) (int a_nrOTUs, char **a_alignment) {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int i;
    $1 = PyList_Size($input);
    $2 = (char **) malloc(($1+1)*sizeof(char *));
    for (i = 0; i < $1; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
	$2[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError,"list must contain strings");
	free($2);
	return NULL;
      }
    }
    $2[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

%typemap(in) (int a_nrBackboneSets, char **a_backboneSetsList) {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int i;
    $1 = PyList_Size($input);
    $2 = (char **) malloc(($1+1)*sizeof(char *));
    for (i = 0; i < $1; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
	$2[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError,"list must contain strings");
	free($2);
	return NULL;
      }
    }
    $2[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

%typemap(freearg) (int a_nrOTUs, char **a_alignment) {
  free((char *) $2);
}

%typemap(freearg) (int a_nrBackboneSets, char **a_backboneSetsList) {
  free((char *) $2);
}

%exception {
   Py_BEGIN_ALLOW_THREADS	
   try {
      $action
   } catch (std::exception &e) {
      Py_BLOCK_THREADS
      PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(e.what()));
      return NULL;
   }
   Py_END_ALLOW_THREADS
}

%module cConstrainedNJlib
%{
extern void initialize(int dim);
extern char *computeTree(int a_nrOTUs, char **a_alignment, int a_nrBackboneSets, char **a_backboneSetsList, int a_resample);
%}
%include interface.h
