/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"
#include "utilities.h"

/*********** LOCAL STUFF **********/
double loghalffact[50 * ABSMIGMAX];
int BITNUMBERTRUE[256];
unsigned long *iseed;

/* function prototype */
static int *ivector (long nl, long nh);
static void free_ivector (int *v, long nl, long nh);
static unsigned long *lvector (long nl, long nh);
static void free_lvector (unsigned long *v, long nl, long nh);
static double bessi0 (double x);
static double bessi1 (double x);

// stuff called by indexx() - sorting an index,  from NR 
#define NRANSI
#define NR_END 1l
#define FREE_ARG char*
int *
ivector (long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  v = (int *) malloc ((size_t) ((nh - nl + 1 + NR_END) * sizeof (int)));
  if (!v)
  {
    nrerror ("allocation failure in ivector()");
  }
  return v - nl + NR_END;
}

void
free_ivector (int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free ((FREE_ARG) (v + nl - NR_END));
} 

unsigned long *
lvector (long nl, long nh)
//long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;
  v =
    (unsigned long *)
    malloc ((unsigned int) ((nh - nl + 1 + NR_END) * sizeof (long)));
  if (!v)
    nrerror ("allocation failure in lvector()");
  return v - nl + NR_END;
}

void
free_lvector (unsigned long *v, long nl, long nh)
/* XFREE an unsigned long vector allocated with lvector() */
{
  free ((FREE_ARG) (v + nl - NR_END));
}

double
bessi0 (double x)
{
  double ax, ans;
  double y;
  if ((ax = fabs (x)) < 3.75)
  {
    y = x / 3.75;
    y *= y;
    ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                                                       + y * (0.2659732 +
                                                              y *
                                                              (0.360768e-1 +
                                                               y *
                                                               0.45813e-2)))));
  }
  else
  {
    y = 3.75 / ax;
    ans = (exp (ax) / sqrt (ax)) * (0.39894228 + y * (0.1328592e-1
                                                      + y * (0.225319e-2 +
                                                             y *
                                                             (-0.157565e-2 +
                                                              y *
                                                              (0.916281e-2
                                                               + y *
                                                               (-0.2057706e-1
                                                                +
                                                                y *
                                                                (0.2635537e-1
                                                                 +
                                                                 y *
                                                                 (-0.1647633e-1
                                                                  + y *
                                                                  0.392377e-2))))))));
  }
  return ans;
}

double
bessi1 (double x)
{
  double ax, ans;
  double y;
  if ((ax = fabs (x)) < 3.75)
  {
    y = x / 3.75;
    y *= y;
    ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                                                               + y *
                                                               (0.2658733e-1
                                                                +
                                                                y *
                                                                (0.301532e-2
                                                                 +
                                                                 y *
                                                                 0.32411e-3))))));
  }
  else
  {
    y = 3.75 / ax;
    ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
                                                   - y * 0.420059e-2));
    ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
                                                 + y * (0.163801e-2 +
                                                        y * (-0.1031555e-1 +
                                                             y * ans))));
    ans *= (exp (ax) / sqrt (ax));
  }
  return x < 0.0 ? -ans : ans;
}


/********** GLOBAL FUNCTIONS ***********/


#ifdef USE_MYASSERT
//#  ifndef ASSIGNMENT
void
myassert (int isit)
{
  if (!isit)
  {
    /*if (outfile != NULL) can't do this unless outfile is global
       f_close(outfile); */
    exit (-9);
  }
}

//#  endif /* ASSIGNMENT */
#endif /* USE_MYASSERT */


/* Here is how to add an error message.
 * 1. Add an element to enum with prefix of "IMERR_XXX"
 * 2. Add the corresponding error message string to simerrmsg.
 * 3. Call function IM_err (IMERR_XXX, ...) where ... is like the way that
 *    printf function arguments are used.
 */
static const char *simerrmsg[] = {
  /* 0 */ "",
  /* 1 */ "cannot open a file for reading",
  /* 2 */ " memory error ",
  /* 3 */ "problem finding .ti file(s)",
  /* 4 */ "cannot create file",
  /* 5 */ "cannot open file for appending",
  /* 6 */ " problem with opening or closing output file",
  /* 7 */ "",
  /* 8 */ "incompatibility on command line",
  /* 9 */ "not enough information provided on command line",
  /* 10 */ "problem with command line formatting ",
  /* 11 */ "",
  /* 12 */ "missing population string in input file",
  /* 13 */ "some problem in the population tree string ",
  /* 14 */ "",
  /* 15 */ "problem w/ number of loci indicated in data file",
  /* 16 */ "",
  /* 17 */ "",
  /* 18 */ "mutation range priors to constraining - not able to set starting values",
  /* 19 */ "product of mutation scalars not equal to 1",
  /* 20 */ "",
  /* 21 */ "Too many migration events found when storing edges",
  /* 22 */ "to much migration called for in checkmig - reduce maximum value of migration parameter",
  /* 23 */ "problem in calculating HPD interval",
  /* 24 */ "",
  /* 25 */ "error reading data, too many lines or line too long",
  /* 26 */ "error in data",
  /* 27 */ "",
  /* 28 */ "problem reading mcf file",
  /* 29 */ "cannot load mcf file, split times in file not compatiable with t prior",
  /* 30 */ "",
  /* 31 */ "",
  /* 32 */ "",
  /* 33 */ "problem with root ",
  /* 34 */ "",
  /* 35 */ "input file invalid ",
  /* 36 */ "data not compatible with infinite sites model",
  /* 37 */ "",
  /* 38 */ "likelihoods do not add up for stepwise model",
  /* 39 */ "",
  /* 40 */ "",
  /* 41 */ "error in lowergamma function",
  /* 42 */ "error in uppergamma function",
  /* 43 */ "",
  /* 44 */ "",
  /* 45 */ "error calculating prior of t  in multi_t_prior_func",
  /* 45 */ "",
  /* 46 */ "",
  /* 47 */ "",
  /* 48 */ "",
  /* 49 */ "",
  /* 50 */ "",
  /* 51 */ "error in NR functions",
  /* 52 */ "GSL error",
  /* 53 */ "",
  /* 54 */ "",
  /* 55 */ "Potential bug",
  /* 56 */ "Invalid input file: gene name",
  /* 57 */ "Assignment Error"
};

void
IM_err (int i, const char *fmt, ...)
{
  va_list args;
  va_start (args, fmt);
  fprintf (stderr, "IMamp: %s - ", simerrmsg[i]);
  vfprintf (stderr, fmt, args);
  fprintf (stderr, "\n");
  va_end (args);
  exit (i);
}


void
nrerror (const char error_text[])
/* Numerical Recipes standard error handler */
{
  IM_err (IMERR_NUMERICALRECIPES, " error text: %s", error_text);
}

/* for a large value x,  cosh[x] = sinh[x] 
also (cosh[x] + sinh[x] = exp[x]  so for a large value of x  cosh[x] = sinh[x] = Exp[x]/2
so rather than return floating error,  for the log of a cosh of a large number  just return x - log[2]  */
double
mylogcosh (double x)
{
  if (x < 100)
    return log (cosh (x));
  else
    return x - LOG2;
}

double
mylogsinh (double x)
{
  if (x < 100)
    return log (sinh (x));
  else
    return x - LOG2;
}


/***********************************/
/* RANDOM NUMBER RELATED FUNCTIONS */
/***********************************/
static int z_rndu = 137;
static int x_rndu = 11, y_rndu = 23;    // these had been in uniform() but I moved them out to make it easier to reset the random number generator 
void
resetseeds (int seed)           // use for debugging
{
  z_rndu = 137;
  x_rndu = 11;
  y_rndu = 23;
  z_rndu = 170 * (seed % 178) + 137;
  iseed = malloc (sizeof (unsigned long));
  *iseed = ULONG_MAX / 2 + (unsigned long) seed;        // just set it so that iseed points to a large number - probably not necessary
}

void
setseeds (int seed)
{
  z_rndu = 170 * (seed % 178) + 137;
  iseed = malloc (sizeof (unsigned long));
  *iseed = ULONG_MAX / 2 + (unsigned long) seed;        // just set it so that iseed points to a large number - probably not necessary
}

void
unsetseeds ()
{
  XFREE (iseed);
  iseed = NULL;
  return;
}

double
uniform ()
{

/* Rasmus Note:
    U(0,1): AS 183: Appl. Stat. 31:188-190
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.

   Suggested to me by Z. Yang who also provided me with
   the source code used here.
*/
/* different compilers were giving slightly diffferent numbers unless there are float casts */
  // turn this back on  static int x_rndu=11, y_rndu=23;
  double r;
  x_rndu = 171 * (x_rndu % 177) - 2 * (x_rndu / 177);
  y_rndu = 172 * (y_rndu % 176) - 35 * (y_rndu / 176);
  z_rndu = 170 * (z_rndu % 178) - 63 * (z_rndu / 178);
  if (x_rndu < 0)
    x_rndu += 30269;
  if (y_rndu < 0)
    y_rndu += 30307;
  if (z_rndu < 0)
    z_rndu += 30323;
  r = x_rndu / 30269.0 + y_rndu / 30307.0 + z_rndu / 30323.0;
  return (r - (long) r);
}

/* $\lambda e^{-\lambda t}$ where $\lambda$ equal to c */
double
expo (double c)
{
  /*float temp;d
     temp = log(uniform());
     temp = (1.0/(float) c) * temp;
     return (double) - temp; */
  return -log (uniform ()) / c;
}
__forceinline int
randposint (int lessthanval)    // actually should be rand_nonneg_int 
{
  return (int) floor (uniform () * lessthanval);
}

/* for binary random numbers  -	quick - based on Press et al irbit2() */
#define IB1 1
#define IB18 131072
#define MASK 19
__forceinline int
bitran (void)
{
  if (*iseed & IB18)
  {
    *iseed = ((*iseed ^ MASK) << 1) | IB1;
    return 1;
  }
  else
  {
    *iseed <<= 1;
    return 0;
  }
}


#undef MASK
#undef IB18
#undef IB1

#define ONEOVERSQRT2PI  0.3989422803
double
normprob (double mean, double stdev, double val)
{
  return ONEOVERSQRT2PI * exp (-SQR ((val - mean) / stdev) / 2) / stdev;
}


#undef ONEOVERSQRT2PI
double
normdev (double mean, double stdev)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  double rescale;
  if (iset == 0)
  {
    do
    {
      v1 = 2.0 * uniform () - 1.0;
      v2 = 2.0 * uniform () - 1.0;
      rsq = v1 * v1 + v2 * v2;
    }
    while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt (-2.0 * log (rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    rescale = ((v2 * fac * stdev) + mean);
    return rescale;
  }
  else
  {
    iset = 0;
    rescale = ((gset * stdev) + mean);
    return rescale;
  }
}                               /* normdev */


/* gets a poisson random variable.  checked it out with simulations 
if condition is 0, pick an even number
if condition is 1, pick an odd number
if condition is 2,  pick any number other than 0 (because they must be different)
if condition is 3, pick any number other than 1
if condition is -1 pick any number 
If param > minpp the normal distribution is used as an approximation

Also stuck in some stuff to deal with the case when condition is 1 and param << 1
*/

#define USENORMAL  25.0
#define MINPP 0.25
int
poisson (double param, int condition)
{
  long double randnum, raised, rcheck;
  int i;
  int stop;
  if (param < MINPP && condition == 1)  /* i.e. need an odd number but the parameter is a small value */
  {
    randnum = uniform ();
    rcheck = param / sinh (param);
    if (randnum < rcheck)       /* param / sinh(param) is prob of # being 1, given it must be odd */
      i = 1;
    else

    {
      rcheck = param * (6 + param * param) / (6 * sinh (param));
      if (randnum < rcheck)
        i = 3;

      else
        i = 5;
    }
    return (i);
  }

  /*if (param < 0.25 && condition == 2)  a function for cases when  need any number not zero,  but the parameter is a small value  */
  /* seems to work - check  file Poisson_low_value_parameter_simulation_check.nb */
  if (param < MINPP && condition == 2)
  {
    raised = exp (-param);
    rcheck = raised = param * raised / (1 - raised);
    randnum = uniform ();
    i = 1;
    while (randnum > rcheck)
    {
      raised *= param / (i + 1);
      rcheck += raised;
      i++;
    }
    return (i);
  }

  do
  {
    if (param >= USENORMAL)     // used normal approximation
    {
      i = IMAX (0, POSROUND (normdev (param, param)));

      // maybe should not cast as (int) before POSROUND ??
      //i = IMAX(0,POSROUND((int) normdev(param,param)));
    }
    else
    {
      raised = exp (-param);
      randnum = uniform ();
      for (i = 0; randnum > raised; i++)
        randnum *= uniform ();
    }

    switch (condition)
    {
    case -1:
      stop = 1;
      break;
    case 0:
      stop = !ODD (i);
      break;
    case 1:
      stop = ODD (i);
      break;
    case 2:
      stop = (i != 0);
      break;                    //*
    case 3:
      stop = (i != 1);
      break;
    }
  } while (!stop); 

  if (i > ADDMIGMAX)            /* trap numbers that get too large to handle (and that don't make sense) */
  {
    if (condition == 1)
      i = ADDMIGMAX - 1;
    else
      i = ADDMIGMAX;
  }
  return (i);
}                               /* pickpoisson */


#undef USENORMAL
#undef MINPP
int
geometric (double p)
/* returns a geometrically distributed random integer variable >= 0 */
/* this distribution is given by prob(k) = p*(1-p)^(k-1)  where k = 1,2,...  which has a mode of 1 and an expectation of 1/p */
/* it is a bit different from prob(k) = p*(1-p)^k   where k  = 0,1, 2... , which has a mode of 0 and an expectation of (1-p)/p */
/* the variance of these two different versions is the same, i.e. (1-p)/p^2 */
/* checked this in various ways - seems ok */
{
  return (int) ceil (log (uniform ()) / log (1.0 - p));
}

/***********************************/
/* SORTING FUNCTIONS */
/***********************************/

/* heap sort is slightly faster than shell sort */
void
hpsort (struct gtreeevent *lptr, int n)
{
  unsigned long i, ir, j, l;
  struct gtreeevent t;
  if (n < 2)
    return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;)
  {
    if (l > 1)
    {
      t = *(lptr + --l);
    }
    else
    {
      t = *(lptr + ir);
      *(lptr + ir) = *(lptr + 1);
      if (--ir == 1)
      {
        *(lptr + 1) = t;
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir)
    {
      if (j < ir && (lptr + j)->time < (lptr + (j + 1))->time)
        j++;
      if (t.time < (lptr + j)->time)
      {
        *(lptr + i) = *(lptr + j);
        i = j;
        j <<= 1;
      }
      else
        j = ir + 1;
    }
    *(lptr + i) = t;
  }
}                               /*hpsort */


/* heap sort is slightly faster than shell sort */
void
hpsortreg (double list[], int n)
{
  unsigned long i, ir, j, l;
  double t;
  if (n < 2)
    return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;)
  {
    if (l > 1)
    {
      t = list[--l];
    }
    else
    {
      t = list[ir];
      list[ir] = list[1];
      if (--ir == 1)

      {
        list[1] = t;
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir)
    {
      if (j < ir && list[j] < list[j + 1])
        j++;
      if (t < list[j])

      {
        list[i] = list[j];
        i = j;
        j <<= 1;
      }
      else
        j = ir + 1;
    }
    list[i] = t;
  }
}                               /*hpsort */


/* heap sort is slightly faster than shell sort */
void
hpsortmig (struct migstruct *lptr, int n)
{
  unsigned long i, ir, j, l;
  double t;
  if (n < 2)
    return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;)
  {
    if (l > 1)
    {
      t = (lptr + --l)->mt;
    }
    else
    {
      t = (lptr + ir)->mt;
      (lptr + ir)->mt = (lptr + 1)->mt;
      if (--ir == 1)

      {
        (lptr + 1)->mt = t;
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir)
    {
      if (j < ir && (lptr + j)->mt < (lptr + (j + 1))->mt)
        j++;
      if (t < (lptr + j)->mt)
      {
        (lptr + i)->mt = (lptr + j)->mt;
        i = j;
        j <<= 1;
      }
      else
      {
        j = ir + 1;
      }
    }
    (lptr + i)->mt = t;
  }
}                               /*hpsortmig */
void
shellhist (struct hlists *hptr, int length)
{
  double aln = 1.442695022, tiny = 1.0e-5;
  struct hlists h;
  static int nn, m, lognb2, i, j, k, l;
  lognb2 = (int) floor (log ((double) length) * aln + tiny);
  m = length;
  for (nn = 1; nn <= lognb2; nn++)
  {
    m = m / 2;
    k = length - m;
    for (j = 0; j <= k - 1; j++)
    {
      i = j;
    reloop:l = i + m;
      if (((hptr + l)->p < (hptr + i)->p)
          || (((hptr + l)->p == (hptr + i)->p)
              && (((hptr + l)->v < (hptr + i)->v))))
      {
        h = *(hptr + i);
        *(hptr + i) = *(hptr + l);
        *(hptr + l) = h;
        i = i - m;
        if (i >= 0)
          goto reloop;
      }
    }
  }
}                               /* shellhist */


/* quicksort of an index of locations */
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp
#define M 7
#define NSTACK 50
void
indexx (unsigned long n, struct gtreeevent *arr, unsigned long *indx)
{
  unsigned long i, indxt, ir = n, itemp, j, k, l = 1;
  int jstack = 0, *istack;
  double a;
  istack = ivector (1, NSTACK);
  for (j = 1; j <= n; j++)
    indx[j] = j;
  for (;;)
  {
    if (ir - l < M)
    {
      for (j = l + 1; j <= ir; j++)
      {
        indxt = indx[j];
        a = arr[indxt].time;
        for (i = j - 1; i >= l; i--)
        {
          if (arr[indx[i]].time <= a)
            break;
          indx[i + 1] = indx[i];
        }
        indx[i + 1] = indxt;
      }
      if (jstack == 0)
        break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else
    {
      k = (l + ir) >> 1;
      SWAP (indx[k], indx[l + 1]);
      if (arr[indx[l]].time > arr[indx[ir]].time)
      {
        SWAP (indx[l], indx[ir]);
      }
      if (arr[indx[l + 1]].time > arr[indx[ir]].time)
      {
        SWAP (indx[l + 1], indx[ir]);
      }
      if (arr[indx[l]].time > arr[indx[l + 1]].time)
      {
        SWAP (indx[l], indx[l + 1]);
      }
      i = l + 1;
      j = ir;
      indxt = indx[l + 1];
      a = arr[indxt].time;
      for (;;)
      {
        do
          i++;
        while (arr[indx[i]].time < a);

        do
          j--;
        while (arr[indx[j]].time > a);

        if (j < i)
          break;
        SWAP (indx[i], indx[j]);
      }
      indx[l + 1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack > NSTACK)
        nrerror ("NSTACK too small in indexx.");
      if (ir - i + 1 >= j - l)
      {
        istack[jstack] = ir;
        istack[jstack - 1] = i;
        ir = j - 1;
      }
      else
      {
        istack[jstack] = j - 1;
        istack[jstack - 1] = l;
        l = i;
      }
    }
  }
  free_ivector (istack, 1, NSTACK);
}


#undef SWAP

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp
//regular quicksort
//void sort(unsigned long n, float arr[]) line from nr  
// slightly modified,  works on array of doubles 
void
sort (double arr[], unsigned long n)
{
  unsigned long i, ir = n, j, k, l = 1, *istack;
  int jstack = 0;
  double a, temp;
  istack = lvector (1, NSTACK);
  for (;;)
  {
    if (ir - l < M)
    {
      for (j = l + 1; j <= ir; j++)
      {
        a = arr[j];
        for (i = j - 1; i >= l; i--)
        {
          if (arr[i] <= a)
            break;
          arr[i + 1] = arr[i];
        }
        arr[i + 1] = a;
      }
      if (jstack == 0)
        break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else
    {
      k = (l + ir) >> 1;
      SWAP (arr[k], arr[l + 1]); 
      if (arr[l] > arr[ir])
      {
        SWAP (arr[l], arr[ir]);
      }
      if (arr[l + 1] > arr[ir])
      {
        SWAP (arr[l + 1], arr[ir]);
      }
      if (arr[l] > arr[l + 1])
      {
        SWAP (arr[l], arr[l + 1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      for (;;)
      {

        do
          i++;
        while (arr[i] < a);

        do
          j--;
        while (arr[j] > a);
        if (j < i)
          break;
        SWAP (arr[i], arr[j]);
      }
      arr[l + 1] = arr[j];
      arr[j] = a;
      jstack += 2;
      if (jstack > NSTACK)
        nrerror ("NSTACK too small in sort.");
      if (ir - i + 1 >= j - l)
      {
        istack[jstack] = ir;
        istack[jstack - 1] = i;
        ir = j - 1;
      }
      else
      {
        istack[jstack] = j - 1;
        istack[jstack - 1] = l;
        l = i;
      }
    }
  }
  free_lvector (istack, 1, NSTACK);
}


#undef SWAP
#undef M
#undef NSTACK
#undef SWAP

/***********************************/
/* GAMMA DISTRIBUTION FUNCTIONS */
/***********************************/
double
gammln (double xx)
{
  double x, y, tmp, ser;
  static double cof[6] =
    { 76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
  };
  int j;
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log (tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++)
    ser += cof[j] / ++y;
  return -tmp + log (2.5066282746310005 * ser / x);
}


/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

#define ITMAX 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
void
gcf (double *gammcf, double a, double x, double *gln)
{
  int i;
  double an, b, c, d, del, h;
  *gln = logfact[(int) a - 1];
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  for (i = 1; i <= ITMAX; i++)
  {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs (d) < FPMIN)
      d = FPMIN;
    c = b + an / c;
    if (fabs (c) < FPMIN)
      c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs (del - 1.0) < EPS)
      break;
  }
  if (i > ITMAX)
    IM_err (IMERR_LOWERGAMMA, " too many iterations within gcf()");
  *gammcf = exp (-x + a * log (x) - (*gln)) * h;
}


// gcflog is the same as gcf but returns the logarithm of gammcf
// saves a little time and stops some underflows
void
gcflog (double *gammcflog, double a, double x, double *gln)
{
  int i;
  double an, b, c, d, del, h;
  *gln = logfact[(int) a - 1];
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  for (i = 1; i <= ITMAX; i++)
  {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs (d) < FPMIN)
      d = FPMIN;
    c = b + an / c;
    if (fabs (c) < FPMIN)
      c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs (del - 1.0) < EPS)
      break;
  }
  if (i > ITMAX)
    IM_err (IMERR_UPPERGAMMA, " too many iterations within gcflog()");

  //*gammcf=exp(-x+a*log(x)-(*gln))*h;
  *gammcflog = (-x + a * log (x) - (*gln)) + log (h);
}                               //gcflog


#undef FPMIN

/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */
void
gser (double *gamser, int a, double x, double *gln)
{
  int n;
  double sum, del, ap;
  *gln = logfact[a - 1];
  if (x <= 0.0)
  {
    if (x < 0.0)
      IM_err (IMERR_UPPERGAMMA,
              " gser() called with negative x value, x: %lf ", x);
    *gamser = 0.0;
    return;
  }
  else
  {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= ITMAX; n++)
    {
      ++ap;
      del *= x / ap;
      sum += del;
      if (fabs (del) < fabs (sum) * EPS)
      {
        *gamser = sum * exp (-x + a * log (x) - (*gln));
        return;
      }
    }
    IM_err (IMERR_UPPERGAMMA, " too many iterations within gser() ");
    return;
  }
}


// gserlog is the same as gser but returns the logarith of gamser
// saves a little time and stops some underflows
void
gserlog (double *gamserlog, int a, double x, double *gln)
{
  int n;
  double sum, del, ap;
  *gln = logfact[a - 1];
  if (x <= 0.0)
  {
    if (x < 0.0)
      IM_err (IMERR_LOWERGAMMA,
              " function gserlog() called from lowergamma() with negative x value, x: %lf",
              x);
    *gamserlog = 0.0;
    return;
  }
  else
  {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= ITMAX; n++)
    {
      ++ap;
      del *= x / ap;
      sum += del;
      if (fabs (del) < fabs (sum) * EPS)
      {
        //*gamser=sum*exp(-x+a*log(x)-(*gln));
        *gamserlog = log (sum) + (-x + a * log (x) - (*gln));
        return;
      }
    }
    IM_err (IMERR_LOWERGAMMA, " too many iterations within gserlog() ");
    return;
  }
}


#undef ITMAX
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

#define MAXIT 100
#define EULER 0.5772156649
#define FPMIN 1.0e-30
double
expint (int n, double x, int *islog)
{
  int i, ii, nm1;
  double a, b, c, d, del, fact, h, psi, ans;
  *islog = 0;
  nm1 = n - 1;
  if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
  {
    IM_err (IMERR_UPPERGAMMA,
            " expint() called from uppergamma() with bad value(s). n %d  x %lf",
            n, x);
  }
  else
  {
    if (n == 0)
    {
      ans = exp (-x) / x;
    }
    else
    {
      if (x == 0.0)
      {
        ans = 1.0 / nm1;
      }
      else
      {
        if (x > 1.0)
        {
          b = x + n;
          c = 1.0 / FPMIN;
          d = 1.0 / b;
          h = d;
          for (i = 1; i <= MAXIT; i++)
          {
            a = -i * (nm1 + i);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            del = c * d;
            h *= del;
            if (fabs (del - 1.0) < EPS)
            {
              *islog = 1;

              //ans = h*exp(-x);
              ans = log (h) - x;
              return ans;
            }
          }

          IM_err (IMERR_UPPERGAMMA,
                  " too many iterations in expint() called from uppergamma()");
        }
        else
        {
          ans = (nm1 != 0 ? 1.0 / nm1 : -log (x) - EULER);
          fact = 1.0;
          for (i = 1; i <= MAXIT; i++)
          {
            fact *= -x / i;
            if (i != nm1)
            {
              del = -fact / (i - nm1);
            }
            else
            {
              psi = -EULER;
              for (ii = 1; ii <= nm1; ii++)
                psi += 1.0 / ii;
              del = fact * (-log (x) + psi);
            }
            ans += del;
            if (fabs (del) < fabs (ans) * EPS)
              return ans;
          }
          IM_err (IMERR_UPPERGAMMA,
                  " too many iterations in expint() called from uppergamma()");
        }
      }
    }
  }
  return ans;
}                               //expint


#undef MAXIT
#undef EPS
#undef FPMIN
#undef EULER
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */
double
uppergamma (int a, double x)
/*  Returns the log of what Mathematica calls the incomplete gamma function.   */
{
  int logindicator;
  double gamser, gammcf, gln, p;
  double temp;
  if (x < 0.0 || a < 0.0)
    IM_err (IMERR_UPPERGAMMA, "  step %d, uppergamma arguments: a  %d, x %lf",
            step, a, x);
  if (a == 0)
  {
    temp = expint (1, x, &logindicator);
    if (!logindicator)
      p = log (temp);

    else
      p = temp;
  }
  else
  {
    if (x < (a + 1.0))
    {
      gser (&gamser, a, x, &gln);
      p = gln + log (1.0 - gamser);
    }
    else
    {
      //gcf(&gammcf,a,x,&gln);
      //p = gln + log(gammcf);
      gcflog (&gammcf, (double) a, x, &gln);
      p = gln + gammcf;
    }
  }
  if (p < -1e200)
    p = -1e200;
  return p;
}                               //uppergamma

double
lowergamma (int a, double x)
/*  modified numrec functions to return the log of  what mathematica would call the complement of 
 the incomplete gamma function.  */
{
  double gamser, gammcf, gln, p;
  if (x < 0.0 || a < 0.0)
    IM_err (IMERR_LOWERGAMMA, "  step %d, lowergamma arguments: a  %d, x %lf",
            step, a, x);
  if (a == 0)
  {
    IM_err (IMERR_LOWERGAMMA, "  step %d, lowergamma arguments: a  %d, x %lf",
            step, a, x);
  }
  if (x < (a + 1.0))
  {

    //gser(&gamser,a,x,&gln);
    //p = gln + log(gamser);
    gserlog (&gamser, a, x, &gln);
    p = gln + gamser;
  }
  else
  {
    gcf (&gammcf, (double) a, x, &gln);
    p = gln + log (1 - gammcf);
  }
  if (p < -1e200)
    p = -1e200;
  return p;
}


/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

/***********************************/
/* TEXT, CHARACTER RELATED FUNCTIONS */
/***********************************/

/* Delete the substring of length "len" at index "pos" from "s". Delete less if out-of-range. */
void
strdelete (char *s, int pos, int len)
{
  int slen;
  if (--pos < 0)
    return;
  slen = (int) (strlen (s) - pos);
  if (slen <= 0)
    return;
  s += pos;
  if (slen <= len)
  {
    *s = 0;
    return;
  }
  while ((*s = s[len]))
    s++;
}


/* insert source into dest at pos */
void
strinsert (char *dest, char *source, int pos)
{
  char temp[POPTREESTRINGLENGTHMAX];
  temp[0] = '\0';
  if (pos > 0)
    strncpy (temp, dest, (size_t) pos);
  temp[pos] = '\0';
  strcat (temp, source);
  if ((int) strlen (dest) > pos)
    strcat (temp, &dest[pos - 1]);
  strcpy (dest, temp);
}


/* remove all characters, c, from a string */
void
strremove (char *s, char c)
{
  char *cspot;
  while ((cspot = strchr (s, c)) != NULL)
  {
    *cspot = '\0';
    strcat (s, cspot + 1);
  }
}

/* truncate a string at the first instance of c */
void
strtrunc (char *s, char c)
{
  char *cspot;
  if (strchr (s, c))
  {
    cspot = strchr (s, c);
    *cspot = '\0';
  }
}


/* find next whitespace, after next non-whitespace */
char *
nextwhite (char *c)
{
  int nonw;
  nonw = !isspace (*c);
  while (*c != '\0')
  {
    while ((nonw == 1 && !isspace (*c)) || (nonw == 0 && isspace (*c)))
    {
      c++;
      if (nonw == 0 && !isspace (*c))
        nonw = 1;
    }
    return c;
  }
  return c;
}                               /* nextwhite */


/* finds the next non-space character after the next space */
char *
nextnonspace (char *textline)
{
  char *cc;
  if (textline == NULL)
    return NULL;
  cc = textline;
  while (*cc != ' ' && *cc != '\0')
    cc++;
  while (*cc == ' ')
    cc++;
  if (*cc == '\0')
    return NULL;
  else
    return cc;
}

/*
void
skip_a_line (FILE * f)
{
  int c;
  c = 'x';
  while (c != '\n')
  {
    c = fgetc (f);
  }
  return;
}
*/

int
comma_exists (char *s)
{
  int l;
  int i;
  int v;

  v = 0;
  l = strlen (s);
  for (i = 0; i < l; i++)
  {
    if (s[i] == ',' || s[i] == ';')
    {
      v = 1;
    }
  }
  return v;
}

/* remove spaces and count parentheses */
void
checktreestring (char *t)
{
  int closepcheck, openpcheck, i;

  for (closepcheck = 0, openpcheck = 0, i = 0; (unsigned) i < strlen (t); i++)  // check if the number of open parenthese match the close parentheses
  {
    closepcheck += t[i] == ')';
    openpcheck += t[i] == '(';
    if ((t[i] == ' ') || (t[i] == '\t') || (t[i] == '\n') || (t[i] == '\r'))    // remove spaces
    {
      t[i] = '\0';
      strcat (t, &(t[i + 1]));
      i--;
    }
  }
  if (closepcheck != openpcheck)
  {
    //  errr (-1, -1, 17, closepcheck, openpcheck);
    IM_err (IMERR_POPTREESTRINGFAIL,
            " wrong number of parentheses, string %s,  open '(' %d  close ')' , step %d",
            t, openpcheck, closepcheck, step);
  }
}                               /* checktreestring */

/**********************************************************************
 * wildcard.c
 *
 * written by Julian Robichaux, http://www.nsftools.com
 **********************************************************************/

#ifndef FALSE
#define FALSE 0
#endif /*  */

#ifndef TRUE
#define TRUE !FALSE
#endif /*  */

/*  See if a string matches a wildcard specification that uses * or ?
    (like "*.txt"), and return TRUE or FALSE, depending on the result.
    There's also a TRUE/FALSE parameter you use to indicate whether
    the match should be case-sensitive or not.  */
int
IsWildcardMatch (char *wildcardString, char *stringToCheck, int caseSensitive)
{
  char wcChar;
  char strChar;

  // use the starMatchesZero variable to determine whether an asterisk
  // matches zero or more characters (TRUE) or one or more characters
  // (FALSE)
  int starMatchesZero = TRUE;
  while ((strChar = *stringToCheck) && (wcChar = *wildcardString))
  {
    // we only want to advance the pointers if we successfully assigned
    // both of our char variables, so we'll do it here rather than in the
    // loop condition itself
    *stringToCheck++;
    *wildcardString++;
    // if this isn't a case-sensitive match, make both chars uppercase
    // (thanks to David John Fielder (Konan) at http://innuendo.ev.ca
    // for pointing out an error here in the original code)
    if (!caseSensitive)
    {
      wcChar = (char) toupper (wcChar);
      strChar = (char) toupper (strChar);
    }

    // check the wcChar against our wildcard list
    switch (wcChar)
    {
      // an asterisk matches zero or more characters
    case '*':
      // do a recursive call against the rest of the string,
      // until we've either found a match or the string has
      // ended
      if (starMatchesZero)
        *stringToCheck--;
      while (*stringToCheck)
      {
        if (IsWildcardMatch (wildcardString, stringToCheck++, caseSensitive))
          return TRUE;
      }
      break;
      // a question mark matches any single character
    case '?':
      break;
      // if we fell through, we want an exact match
    default:
      if (wcChar != strChar)
        return FALSE;
      break;
    }
  }

  // if we have any asterisks left at the end of the wildcard string, we can
  // advance past them if starMatchesZero is TRUE (so "blah*" will match "blah")
  while ((*wildcardString) && (starMatchesZero))
  {
    if (*wildcardString == '*')
      wildcardString++;
    else
      break;
  }

  // if we got to the end but there's still stuff left in either of our strings,
  // return false; otherwise, we have a match
  if ((*stringToCheck) || (*wildcardString))
    return FALSE;
  else
    return TRUE;
}


#undef TRUE
#undef FALSE

__forceinline void
checkmigt (int n, struct edge *gtree)
{
  int cmm;
  struct migstruct **mig;

  /* if we call for more memory for global miginfo, we just return
   * because there must be enough memory for that. */
  if (n + 2 > gtree->cmm)
  {
    assert (gtree != NULL);
    mig = &gtree->mig;
    cmm = gtree->cmm;

    cmm += MIGINC;
    if (cmm > ABSMIGMAX)
    {
      IM_err(IMERR_MIGARRAYTOOBIG," too many migrations in array %d",cmm);
    }
    *mig = realloc (*mig, cmm * sizeof (struct migstruct));
    gtree->cmm = cmm;
  }
  return;
}

__forceinline void
copymig (struct migstruct *m1, struct migstruct *m2)
{
  m1->mp = m2->mp;
  m1->mt = m2->mt;
}

/* realloc memory for migration arrays, as needed */
__forceinline void
checkmig (int i, struct migstruct **mig, int *nmig)
{
  assert (*nmig >= MIGINC);
  if (i + 2 > *nmig)
  {
    *nmig += MIGINC;
    if (*nmig > ABSMIGMAX)
    {
      IM_err(IMERR_MIGARRAYTOOBIG," too many migrations in array %d ",*nmig);
    }
    *mig = realloc (*mig, *nmig * sizeof (struct migstruct));
  }
}



/************************/
/* Other Misc Functions */
/************************/

/** 
 * logfact is absurdly large, but this is because it must be usable 
 * for all of the migration events in a * data set, which might be very large
 * for a large data set with large upper bounds on migration rates.
 * We need log factorial of n + 1/2 where n a is nonnegative integer.
 * loghalffact (n) := log (n + 1/2)! where n starts from 0 to 50 * ABSMIGMAX.
 */
#ifndef M_LNPI
#  define M_LNPI 1.14472988584940017414342735135
#endif
#ifndef M_LN2
#  define M_LN2 0.69314718055994530941723212146
#endif
void
setlogfact (void)
{
  int i;
  UByteP a;
  // int bitj;
  int ic;
  logfact[0] = 0;
  for (i = 1; i <= 100 * ABSMIGMAX; i++)
    logfact[i] = ((double) logfact[i - 1]) + log ((double) i);
  for (i = 0; i < 50 * ABSMIGMAX; i++)
  {
    loghalffact[i] = M_LNPI / 2
      + logfact[2 * i + 2] - logfact[i + 1] - (2 * i + 2) * M_LN2;
  }
  a = (UByteP) malloc (sizeof (unsigned char));
  for (ic = 0; ic < 256; ic++)
  {
    BitZero (a, 1);
    a[0] = (unsigned char) ic;
    BITNUMBERTRUE[ic] = 0;
    for (i = 0; i < 8; i++)
    {
      if (BitIsTrue (a, i))
      {
        BITNUMBERTRUE[ic] += 1;
      }
    }
  }
  XFREE (a);
  a = NULL;
}

void
ieevent (struct eevent *a)      // initialize structure
{
  a->n = 0;
  a->s = 0;
  a->ss = 0;
  a->s2 = 0;
}


double
rnd (double x, int n)
/* sort of rounds  - quite crude*/
/* only use for values in which x*10^n  does not exceed a long int */
{
  int i, k;
  long xi;
  k = 0;
  while (fabs (x) < 1)
  {
    x *= 10;
    k++;
  }
  while (fabs (x) >= 10)
  {
    x /= 10;
    k--;
  }
  for (i = 0; i < n; i++)
    x *= 10;
  xi = (long) ((x) + 0.5);
  x = (double) xi;
  if (k < 0)
  {
    for (i = 0; i > k; i--)
      x *= 10;
  }
  if (k > 0)
  {
    for (i = 0; i < k; i++)
      x /= 10;
  }
  for (i = 0; i < n; i++)
    x /= 10;
  return x;
}


#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
double
bessi (int n, double x)
{
  int j;
  double bi, bim, bip, tox, ans;
  assert (x >= 0);
  n = abs (n);
  if (x > 700)
    return (MYDBL_MAX);
  if (n == 0)
    return bessi0 (x);
  if (n == 1)
    return bessi1 (x);
  if (x == 0.0)
  {
    return 0.0;
  }
  else
  {
    tox = 2.0 / fabs (x);
    bip = ans = 0.0;
    bi = 1.0;
    for (j = 2 * (n + (int) sqrt (ACC * n)); j > 0; j--)
    {
      bim = bip + j * tox * bi;
      bip = bi;
      bi = bim;
      if (fabs (bi) > BIGNO)
      {
        ans *= BIGNI;
        bi *= BIGNI;
        bip *= BIGNI;
      }
      if (j == n)
        ans = bip;
    }
    ans *= bessi0 (x) / bi;
    return x < 0.0 && (n & 1) ? -ans : ans;
  }
}


#undef ACC
#undef BIGNO
#undef BIGNI

/* eexp is used, generally together with struct extendnum to calculated exponentials of sums of numbers 
that may have very large exponents */ 
#define LOG_10_2  0.30102999566398119521
__forceinline void
eexp (double x, double *m, int *z)
{
  double u, zr, temp;
  int n;
  double us[10] = { 1.0, 0.5, 0.16666666666666666666666666667,
  0.04166666666666666666666666667, 0.00833333333333333333333333333,
  0.001388888888888888888888888889, 0.000198412698412698412698412698,
  0.000024801587301587301587301587301, 2.75573192239858906525573192239859e-6,
  2.75573192239858906525573192239e-7};
  n = (int) floor (x / LOG2);
  zr = LOG_10_2 * (double) n;
  *z = (int) zr;
  zr -= (double) *z;
  u = x - (((double) n) * LOG2);
  temp =
    1 + u * (us[0] +
             u * (us[1] +
                  u * (us[2] +
                       u * (us[3] +
                            u * (us[4] +
                                 u * (us[5] +
                                      u * (us[6] +
                                           u * (us[7] +
                                                u * (us[8] +
                                                     u * (us[9]))))))))));
  *m = temp * pow (10.0, zr);
  if (fabs (*m) > 10)
  {
    *m = *m / 10;
    *z = *z + 1;
  }

  if (fabs (*m) < 1)
  {
    *m = *m * 10;
    *z = *z - 1;
  }
}                               //eexp

double ****
alloc4Ddouble (int d1, int d2, int d3, int d4)
{
  double ****a;
  int i, j, k;

  if ((a = (double ****) malloc (d1 * sizeof (double ***))) == NULL)
    IM_err (IMERR_MEM, "  alloc4Ddouble main malloc did not work.  step %d",
            step);
  for (i = 0; i < d1; i++)
  {
    if ((a[i] = (double ***) malloc (d2 * sizeof (double **))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc4Ddouble first loop malloc did not work. d1 %d,  d2 %d, step %d",
              d1, d2, step);
    for (j = 0; j < d2; j++)
    {
      if ((a[i][j] = (double **) malloc (d3 * sizeof (double *))) == NULL)
        IM_err (IMERR_MEM,
                "  alloc4Ddouble second loop malloc did not work. d2 %d,  d3 %d, step %d",
                d2, d3, step);
      for (k = 0; k < d3; k++)
        if ((a[i][j][k] = (double *) malloc (d4 * sizeof (double))) == NULL)
          IM_err (IMERR_MEM,
                  "  alloc4Ddouble third loop malloc did not work. d3 %d,  d4 %d, step %d",
                  d3, d4, step);
    }
  }
  return a;
}

int ****
alloc4Dint (int d1, int d2, int d3, int d4)
{
  int ****a;
  int i, j, k;

  if ((a = (int ****) malloc (d1 * sizeof (int ***))) == NULL)
    IM_err (IMERR_MEM, "  alloc4Dint main malloc did not work.  step %d",
            step);
  for (i = 0; i < d1; i++)
  {
    if ((a[i] = (int ***) malloc (d2 * sizeof (int **))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc4Dint first loop malloc did not work. d1 %d,  d2 %d, step %d",
              d1, d2, step);
    for (j = 0; j < d2; j++)
    {
      if ((a[i][j] = (int **) malloc (d3 * sizeof (int *))) == NULL)
        IM_err (IMERR_MEM,
                "  alloc4Dint second loop malloc did not work. d2 %d,  d3 %d, step %d",
                d2, d3, step);
      for (k = 0; k < d3; k++)
        if ((a[i][j][k] = (int *) malloc (d4 * sizeof (int))) == NULL)
          IM_err (IMERR_MEM,
                  "  alloc4Dint third loop malloc did not work. d3 %d,  d4 %d, step %d",
                  d3, d4, step);
    }
  }
  return a;
}

void
free4D (void ****a, int d1, int d2, int d3)
{
  int i, j, k;
  for (i = 0; i < d1; i++)
  {
    for (j = 0; j < d2; j++)
    {
      for (k = 0; k < d3; k++)
      {
        XFREE (a[i][j][k]);
        a[i][j][k] = NULL;
      }
      XFREE (a[i][j]);
      a[i][j] = NULL;
    }
    XFREE (a[i]);
    a[i] = NULL;
  }
  XFREE (a);
  a = NULL;
}

double ***
alloc3Ddouble (int layers, int rows, int cols)
{
  double ***a;
  int i, j;
  if ((a = malloc (layers * sizeof (*a))) == NULL)
    IM_err (IMERR_MEM, "  alloc3Ddouble main malloc did not work.  step %d",
            step);
  for (i = 0; i < layers; i++)
  {
    if ((a[i] = malloc (rows * sizeof (*a[i]))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc3Ddouble first loop malloc did not work. layers %d,  rows %d, step %d",
              layers, rows, step);
    for (j = 0; j < rows; j++)
      if ((a[i][j] = malloc (cols * sizeof (*a[i][j]))) == NULL)
        IM_err (IMERR_MEM,
                "  alloc3Ddouble second loop malloc did not work. rows %d,  cols %d, step %d",
                rows, cols, step);
  }
  return a;
}

int ***
alloc3Dint (int layers, int rows, int cols)
{
  int ***a;
  int i, j;
  if ((a = malloc (layers * sizeof (*a))) == NULL)
    IM_err (IMERR_MEM, "  alloc3Dint main malloc did not work.  step %d",
            step);
  for (i = 0; i < layers; i++)
  {
    if ((a[i] = malloc (rows * sizeof (*a[i]))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc3Dint first loop malloc did not work. layers %d,  rows %d, step %d",
              layers, rows, step);
    for (j = 0; j < rows; j++)
      if ((a[i][j] = malloc (cols * sizeof (*a[i][j]))) == NULL)
        IM_err (IMERR_MEM,
                "  alloc3Dint second loop malloc did not work. rows %d,  cols %d, step %d",
                rows, cols, step);
  }
  return a;
}

void
free3D (void ***a, int layers, int rows)
{
  int i, j;
  for (i = 0; i < layers; i++)
    for (j = 0; j < rows; j++)
      XFREE (a[i][j]);
  XFREE (a[i]);
  XFREE (a);
}

int **
alloc2Dint (int rows, int cols)
{
  int **a;
  int i;
  if ((a = malloc (rows * sizeof (*a))) == NULL)
    IM_err (IMERR_MEM, "  alloc2Dint main malloc did not work.  step %d",
            step);
  for (i = 0; i < rows; i++)
    if ((a[i] = malloc (cols * sizeof (*a[i]))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc2Dint loop malloc did not work. rows %d,  cols %d, step %d",
              rows, cols, step);
  return a;
}

double **
alloc2Ddouble (int rows, int cols)
{
  double **a;
  int i;
  if ((a = malloc (rows * sizeof (*a))) == NULL)
    IM_err (IMERR_MEM, "  alloc2Ddouble main malloc did not work.  step %d",
            step);
  for (i = 0; i < rows; i++)
    if ((a[i] = malloc (cols * sizeof (*a[i]))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc2Ddouble loop malloc did not work. rows %d,  cols %d, step %d",
              rows, cols, step);
  return a;
}

void
free2D (void **a, int rows)
{
  int i;
  for (i = 0; i < rows; i++)
    XFREE (a[i]);
  XFREE (a);
}

/********* PROFILING FUNCTIONS - SIMPLE  ***************/
/*  some simple profiling to check function timings 
be sure to define IMPROFILE  
and to use ifdef and endif whereever 

to invoke,  
to include a function:
add a clock_t variable in IMamp.h  
e.g. gextern  clock_t  clock_myfunction

initialize the varible and give it a string by entering a line in the init_clock()  function
e.g.  init_clockinfo(&clock_myfunction,"myfunction");

to record times,  use startclock() and addclock() 

if the function to be timed has only one exit point,  then the start and stop calls can be included in the function
otherwise it is easiest to put the startclock() and addclock() calls flanking the call to that function
e.g.
#ifdef IMPROFILE
	startclock(&clock_myfunction);
#endif 
		myfunction();
#ifdef IMPROFILE
	addclock(&clock_myfunction);
#endif 

*/

#ifdef IMPROFILE

void
startclock (struct clockinfo *c)
{
  c->b = clock ();
}

void
addclock (struct clockinfo *c)
{
  c->c += clock () - c->b;
  c->calls++;
}

void
printclock (FILE * outfile, struct clockinfo c)
{
  FP " %-30s   %8ld       %8d     %10.2lf      %5.2lf      %.7lf\n",
    c.name, c.calls, c.c, (double) c.c / CLOCKS_PER_SEC,
    (double) c.c / c.calls,
    (double) ((double) c.c / CLOCKS_PER_SEC) / c.calls);
}

void init_clockinfo (struct clockinfo *c, char name[])
{
  c->c = 0;
  c->calls = 0;
  strcpy (&c->name[0], &name[0]);
}

void init_clock (void)
{
  init_clockinfo (&clock_updategenealogy, "updategenealogy");
  init_clockinfo (&clock_updategenealogy_covar, "updategenealogy_covar");
  init_clockinfo (&clock_changet_NW, "changet_NW");
  init_clockinfo (&clock_changet_RY1, "changet_RY1");
  init_clockinfo (&clock_changet_RY2, "changet_RY2");
  init_clockinfo (&clock_indexx, "indexx");
  init_clockinfo (&clock_treeweight, "treeweight");
  init_clockinfo (&clock_integrate_tree_prob, "integrate_tree_prob");
  init_clockinfo (&clock_likelihoodIS, "likelihoodIS");
  init_clockinfo (&clock_record, "record");
  init_clockinfo (&clock_swapchains, "swapchains");
  init_clockinfo (&clock_changeu, "changeu");
}


void print_timings (FILE * outfile)
{
  FP "\n");
  FP "===============\n");
  FP " FUNCTION TIMES \n");
  FP "===============\n");
  FP "function                             calls         cycles        seconds  cycles/call  seconds/call\n");
  printclock (outfile, clock_updategenealogy);
  printclock (outfile, clock_updategenealogy_covar);
  printclock (outfile, clock_changet_NW);
  printclock (outfile, clock_changet_RY1);
  printclock (outfile, clock_changet_RY2);
  printclock (outfile, clock_indexx);
  printclock (outfile, clock_treeweight);
  printclock (outfile, clock_integrate_tree_prob);
  printclock (outfile, clock_likelihoodIS);
  printclock (outfile, clock_record);
  printclock (outfile, clock_swapchains);
  printclock (outfile, clock_changeu);

  FP "\n");
}
#endif

/* Function logsum approximately computes the value of v of the following
 * equation:
 * exp(v) = exp(a_1) + exp(a_2) + ...
 * The first arguement is the number of elements of the sum, and the rest of the
 * argements are the n elements. For example, we can compute 
 * exp(3.0) + exp(1.5) + exp(2.5), which is exp(3.604131) by using function call
 * logsum (3, 3.0, 1.5, 2.5). The returned value should be 3.604131.
 */
#define LOG_DBL_MAX    7.0978271289338397e+02
double logsum (int n, ...)
{
  register int i;
  double a;
  double b;
  va_list ap;

  assert (n > 1);
  va_start (ap, n);
  a = va_arg (ap, double);
  for (i = 1; i < n; i++)
  {
    b = va_arg (ap, double);

    if (a > b)
    {
      if (a - b < LOG_DBL_MAX)
        a = b + log (exp (a - b) + 1.0);
    }
    else
    {
      if (b - a < LOG_DBL_MAX)
        a = a + log (exp (b - a) + 1.0);
      else
        a = b;
    }
  }
  va_end (ap);
  return a;
}

double logsuma (int n, double *d) 
{
  register int i;
  double a;
  double b;

  assert (n > 1);
  a = d[0];
  for (i = 1; i < n; i++)
  {
    b = d[i];

    if (a > b)
    {
      if (a - b < LOG_DBL_MAX)
        a = b + log (exp (a - b) + 1.0);
    }
    else
    {
      if (b - a < LOG_DBL_MAX)
        a = a + log (exp (b - a) + 1.0);
      else
        a = b;
    }
  }
  return a;
}

#undef LOG_DBL_MAX

/* SASNGCHUL: Wed Apr 15 2009
 * I borrowed some codes from MICSAT for computing cumulative normal distribution. 
 * This is for implementing their Branch-Swapping algorithm. I have not decided this 
 * would be included or not. */
/* Frequently used numerical constants: */
#define OneUponSqrt2Pi .39894228040143267794
#define twopi 6.283195307179587
#define LnSqrt2Pi -0.9189385332046727417803296 
#define SQRT2 1.414213562373095049
#define SQRTPI 1.772453850905516027

/* ---------------------------------------------------------------------------

   UNIVARIATE NORMAL PROBABILITY

   ---------------------------------------------------------------------------*/

#define UPPERLIMIT 100.0 
/* I won't return either of univariate normal density or 
	probability when x < -UPPERLIMIT  or x > UPPERLIMIT. */

#define P10 242.66795523053175
#define P11 21.979261618294152
#define P12 6.9963834886191355
#define P13 -.035609843701815385
#define Q10 215.05887586986120
#define Q11 91.164905404514901
#define Q12 15.082797630407787
#define Q13 1.0

#define P20 300.4592610201616005
#define P21 451.9189537118729422
#define P22 339.3208167343436870
#define P23 152.9892850469404039
#define P24 43.16222722205673530
#define P25 7.211758250883093659
#define P26 .5641955174789739711
#define P27 -.0000001368648573827167067
#define Q20 300.4592609569832933
#define Q21 790.9509253278980272
#define Q22 931.3540948506096211
#define Q23 638.9802644656311665
#define Q24 277.5854447439876434
#define Q25 77.00015293522947295
#define Q26 12.78272731962942351
#define Q27 1.0

#define P30 -.00299610707703542174
#define P31 -.0494730910623250734
#define P32 -.226956593539686930
#define P33 -.278661308609647788
#define P34 -.0223192459734184686
#define Q30 .0106209230528467918
#define Q31 .191308926107829841
#define Q32 1.05167510706793207
#define Q33 1.98733201817135256
#define Q34 1.0

double cumnorm(double ox, double mean, double sd)
{
	int sn;
	double R1, R2,/* R3,*/ y, y2, y3, y4, y5, y6, y7;
	double erf, erfc, z, z2, z3, z4;
	double phi,x;

        x = (ox-mean)/sd;

	if (x < -UPPERLIMIT) return 0.0;
	if (x > UPPERLIMIT) return 1.0;

	y = x / SQRT2;
	if (y < 0) {
		y = -y;
		sn = -1;
	}
	else
		sn = 1;

	y2 = y * y;
	y4 = y2 * y2;
	y6 = y4 * y2;

	if(y < 0.46875) {
		R1 = P10 + P11 * y2 + P12 * y4 + P13 * y6;
		R2 = Q10 + Q11 * y2 + Q12 * y4 + Q13 * y6;
		erf = y * R1 / R2;
		if (sn == 1)
			phi = 0.5 + 0.5*erf;
		else 
			phi = 0.5 - 0.5*erf;
	}
	else 
		if (y < 4.0) {
			y3 = y2 * y;
			y5 = y4 * y;
			y7 = y6 * y;
			R1 = P20 + P21 * y + P22 * y2 + P23 * y3 + 
			    P24 * y4 + P25 * y5 + P26 * y6 + P27 * y7;
			R2 = Q20 + Q21 * y + Q22 * y2 + Q23 * y3 + 
			    Q24 * y4 + Q25 * y5 + Q26 * y6 + Q27 * y7;
			erfc = exp(-y2) * R1 / R2;
			if (sn == 1)
				phi = 1.0 - 0.5*erfc;
			else
				phi = 0.5*erfc;
		}
		else {
			z = y4;
			z2 = z * z;
			z3 = z2 * z;
			z4 = z2 * z2;
			R1 = P30 + P31 * z + P32 * z2 + P33 * z3 + P34 * z4;
			R2 = Q30 + Q31 * z + Q32 * z2 + Q33 * z3 + Q34 * z4;
			erfc = (exp(-y2)/y) * (1.0 / SQRTPI + R1 / (R2 * y2));
			if (sn == 1)
				phi = 1.0 - 0.5*erfc;
			else 
				phi = 0.5*erfc;
		}

	return phi;
}


