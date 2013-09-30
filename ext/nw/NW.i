
// // This tells SWIG to treat an double * argument with name 'OutValue' as
// // an output value.  We'll append the value to the current result which 
// // is guaranteed to be a List object by SWIG.
// 
// %typemap(argout) double *OutValue {
//     PyObject *o, *o2, *o3;
//     o = PyFloat_FromDouble(*$1);
//     if ((!$result) || ($result == Py_None)) {
//         $result = o;
//     } else {
//         if (!PyTuple_Check($result)) {
//             PyObject *o2 = $result;
//             $result = PyTuple_New(1);
//             PyTuple_SetItem($result,0,o2);
//         }
//         o3 = PyTuple_New(1);
//         PyTuple_SetItem(o3,0,o);
//         o2 = $result;
//         $result = PySequence_Concat(o2,o3);
//         Py_DECREF(o2);
//         Py_DECREF(o3);
//     }
// }
// 
// // The typemap works as follows. First, a check is made to see if any previous result
// // exists. If so, it is turned into a tuple and the new output value is concatenated to
// // it. Otherwise, the result is returned normally. For the sample function spam(), there are
// // three output values--meaning that the function will return a 3-tuple of the results.
// // 
// // As written, the function must accept 4 arguments as input values, last two being pointers
// // to doubles. If these arguments are only used to hold output values (and have no meaningful
// // input value), an additional typemap can be written. For example:
//  
// %typemap(in,numinputs=0) double *OutValue(double temp) {
//     $1 = &temp;
// }
// 
// // By specifying numinputs=0, the input value is ignored. However, since the argument still
// // has to be set to some meaningful value before calling C, it is set to point to a local
// // variable temp. When the function stores its output value, it will simply be placed in this
// // local variable. As a result, the function can now be used as follows:
// 
// 
// // int spam(double a, double b, double *OutValue, double *OutValue);
// 
// // can then be used like this:
// // 
// // status, retval1, retval2 = spam(a, b)
// 
// %module NW
// 
// int align(double a, double b, double *OutValue, double *OutValue);
// 

// This tells SWIG to treat an char * argument with name 'OutValue' as
// an output value.  We'll append the value to the current result which 
// is guaranteed to be a List object by SWIG.

%module NW

%typemap(argout) char **OutValue {
    PyObject *o, *o2, *o3;
    o = PyString_FromString(*$1);
    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        if (!PyTuple_Check($result)) {
            PyObject *o2 = $result;
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
    free(*$1);
}

%typemap(in,numinputs=0) char **OutValue(char *temp) {
    $1 = &temp;
}

// This is apparently the C way:
int align(const char *a, const char *b, double gap_init, double gap_ext, double gap_flank, char **OutValue, char **OutValue);

// // and this is the C++ way:
// %{
// extern int align(const char *a, const char *b, char **OutValue, char **OutValue);
// %}




