/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"
#include "update_gtree_common.h"
#include "xtrapbits.h"

/* This file:
  Includes local functions with three locals for declaration:
  declared locally - see static declarations below
  declared in update_gree_common.h
	used by  this file,  update_gtree.c and update_t_NW.c
  declared globally in imamp.h */

#define M_LNPI 1.14472988584940017414342735135
#ifndef M_LN2
#define M_LN2 0.69314718055994530941723212146
#endif /* M_LN2 */
extern double logfact[100 * ABSMIGMAX + 1];
extern double loghalffact[50 * ABSMIGMAX];
extern int BITNUMBERTRUE[256];


#define gsl_sf_lnfact(n) logfact[(n)]
#define IMA_math_lnfact_half(v,n,quotient,remainder) \
  (quotient) = (n) / 2; \
  (remainder) = (n) % 2; \
  if ((remainder) == 0) \
    { \
      (v) = logfact[(quotient)]; \
    } \
  else \
    { \
      (v) = M_LNPI / 2 \
          + logfact[2 * (quotient) + 2] \
          - logfact[(quotient) + 1] \
          - (2 * (quotient) + 2) * M_LN2; \
    }

// variables declared in update_gtree
extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;

/* declarded in update_gtree.h */
/* extern int mrootdrop; not used here */
/* extern double lmrootdrop; not used here */
/* extern struct genealogy_weights holdgweight_updategenealogy; not used here */
/* extern struct genealogy_weights holdallgweight_updategenealogy; not used here */
/* extern struct probcalc holdallpcalc_updategenealogy; not used here */
extern int rootmove;            /* declared in update_gtree.c */
static int holddownA[MAXLINKED];
static int medgedrop;           /* NO USAGE: SHOULD BE DELETED */
static double lmedgedrop;       /* NO USAGE: SHOULD BE DELETED */
static double holdsisdlikeA[MAXLINKED];
struct genealogy holdgtree;

//variables local to this file 
#define ESTARTLENGTH  20        // fairly arbitrary starting size of e and eci
#define EADDLENGTH  10          // fairly aribtrary amount to add when needed

static int elength;
static struct gtreeevent *e; /* Can we change this to other name? */
static unsigned long *eci;
static struct edge *copyedge;
static double *nnminus1;

/* prototype of local functions */
static void IMA_initmemory_edgemiginfo ();
static double integrate_coalescent_term (int cc, double fc, double hcc,
                                         double max, double min);
static double integrate_migration_term (int cm, double fm, double max,
                                        double min);
static double integrate_migration_term_expo_prior (int cm, double fm,
                                                   double exmean);
static double integrate_splitrate_term (double f, double max, double min);
static int simmpath (int ci, struct edgemiginfo *edgem, int period, int numm,
                     int lastm, double timein, double upt, int pop,
                     int constrainpop);

/********* LOCAL FUNCTIONS **************/

void
IMA_initmemory_edgemiginfo ()
{
  struct edgemiginfo *em;

  if (modeloptions[NOMIGRATION] == 0)
  {
    em = &newedgemig;
    em->mtimeavail = malloc (npops * sizeof (double));
    em->mp = malloc (npops * sizeof (int));
    em->sisid = -1;
    em->li = -1;
    em = &newsismig;
    em->mtimeavail = malloc (npops * sizeof (double));
    em->mp = malloc (npops * sizeof (int));
    em->sisid = -1;
    em->li = -1;
    em = &oldedgemig;
    em->mtimeavail = malloc (npops * sizeof (double));
    em->mp = malloc (npops * sizeof (int));
    em->sisid = -1;
    em->li = -1;
    em = &oldsismig;
    em->mtimeavail = malloc (npops * sizeof (double));
    em->mp = malloc (npops * sizeof (int));
    em->sisid = -1;
    em->li = -1;
  }
  else
  {
    em = &newedgemig;
    em->mtimeavail = malloc (npops * sizeof (double));
    em->mp = malloc (npops * sizeof (int));
    em->sisid = -1;
    em->li = -1;
    em = &newsismig;
    em->mtimeavail = malloc (npops * sizeof (double));
    em->mp = malloc (npops * sizeof (int));
    em->sisid = -1;
    em->li = -1;
    em = &oldedgemig;
    em->mtimeavail = malloc (npops * sizeof (double));
    em->mp = malloc (npops * sizeof (int));
    em->sisid = -1;
    em->li = -1;
    em = &oldsismig;
    em->mtimeavail = malloc (npops * sizeof (double));
    em->mp = malloc (npops * sizeof (int));
    em->sisid = -1;
    em->li = -1;
  }

  return;
}

double
integrate_coalescent_term (int cc, double fc, double hcc, double max,
                           double min)
{                               /* note fc includes inheritance scalar, i.e. it is fc/h  - see treeweight  
                                   the hcc term is cc*log(h) */

  double p;
  double a, b, c, d;
  if (cc > 0)
  {
    assert (fc > 0);
    if (min == 0)
    {
      p = uppergamma ((int) cc - 1, 2 * fc / max) + LOG2 - hcc + (1 - cc) * log (fc);
    }
    else
    {
#ifdef MORESTABLE
      a = uppergamma ((int) cc - 1, 2 * fc / max);
      b = uppergamma ((int) cc - 1, 2 * fc / min);
      LogDiff (p, a, b);
      p += (LOG2 - hcc + (1 - cc) * log (fc));
#else
      p = log (exp (uppergamma ((int) cc - 1, 2 * fc / max)) 
               - exp (uppergamma ((int) cc - 1, 2 * fc / min))) 
          + LOG2 - hcc + (1 - cc) * log (fc);
#endif /* MORESTABLE */
    }
  }
  else
  {
    assert (fabs (hcc) < 1e-10);
    if (2 * fc / max > 0)       /* cc == 0 */
    {
      if (min == 0)
      {
//  fc values coming out very large.  this deals with it,  but why so large ????
#ifdef MORESTABLE
        a = log (max) - 2.0 * fc / max;
        b = LOG2 + log (fc) + uppergamma (0, 2.0 * fc / max);
        LogDiff (p, a, b);
#else
        p = log (max * exp (-2 * fc / max) - 2 * fc * exp (uppergamma (0, 2 * fc / max)));
#endif /* MORESTABLE */
      }
      else
      {
#ifdef MORESTABLE
        a = uppergamma (0, 2 * fc / max);
        b = uppergamma (0, 2 * fc / min);
        LogDiff (c, a, b);
        c += LOG2 + log (fc);
        a = log (max) - 2.0 * fc / max;
        b = log (min) - 2.0 * fc / min;
        LogDiff (d, a, b);
        LogDiff (p, d, c);
#else
        p = log (max * exp (-2 * fc / max) - min * exp (-2 * fc / min)
                 - 2 * fc * (exp (uppergamma (0, 2 * fc / max)) -
                             exp (uppergamma (0, 2 * fc / min))));
#endif /* MORESTABLE */
      }
    }
    else
    {
      p = log (max - min);
    }
  }

  return p;
}                               /* integrate_coalescent_term */

double
integrate_migration_term (int cm, double fm, double max, double min)
{
  double p;
  double a, b, c;
  if (cm > 0)
  {
    assert (fm > 0);
    if (min == 0)
    {
      p = (-1 - cm) * log (fm) + lowergamma ((int) cm + 1, fm * max);
    }
    else
    {

#ifdef MORESTABLE
      a = uppergamma (cm + 1, fm * min);
      b = uppergamma (cm + 1, fm * max);
      LogDiff (c, a, b);
      p = (-1 - cm) * log (fm) + c;
#else
      p = (-1 - cm) * log (fm) 
          + log (exp (uppergamma (cm + 1, fm * min)) 
                 - exp (uppergamma (cm + 1, fm * max)));
#endif /* MORESTABLE */
    }
  }
  else
  {
    if (fm > MPRIORMIN)                 /* cm == 0   use a cutoff, because sometimes the value that comes in is very low, when it should be zero */
    {
      if (min == 0)
      {
        if (max == MPRIORMIN)
        {
          p = 0;
        }
        else
        {
#ifdef MORESTABLE
          a = 0.0;
          b = -fm * max;
          LogDiff (c, a, b);
          p = c - log (fm);
#else
          p = log ((1 - exp (-fm * max)) / fm);
#endif /* MORESTABLE */
        }
      }
      else
      {
#ifdef MORESTABLE
        a = -fm * min;
        b = -fm * max;
        LogDiff (c, a, b);
        p = c - log (fm);
#else
        p = log ((exp (-fm * min) - exp (-fm * max)) / fm);
#endif /* MORESTABLE */
      }
    }
    else
    {
      p = log (max - min);
    }
  }

//      assert(p > -1e200 && p < 1e200);
  return p;
}                               /* integrate_migration_term */

double
integrate_migration_term_expo_prior (int cm, double fm, double exmean)
{
  double p;
    p = -log (exmean) +( -(cm + 1) * log (fm + 1.0/exmean)) + logfact[cm];
  return p;
}                               /* integrate_migration_term_expo_prior */

double
integrate_splitrate_term (double f, double max, double min)
{
  if (min > 0)
  {
    return (2 - npops) * log (f) + logfact[npops] +
      log (exp (uppergamma (npops - 2, f / max)) -
           exp (uppergamma (npops - 2, f / min)));
  }
  else
  {
    return (2 - npops) * log (f) + logfact[npops] + uppergamma (npops - 2,
                                                                f / max);
  }
}                               //integrate_splitrate_term




/* simulate the migration path along the moved edge
this gets called starting from the top, period by period  
'period' is the current period,  numm is the number of migration events so far on this edge 
if it must end up in a particular population,  then constrainpop is that population */
int
simmpath (int ci, struct edgemiginfo *edgem, int period, int numm, int lastm,
          double timein, double upt, int pop, int constrainpop)
{
  int i, lastpop, startm;
  int li;
  int ei;
  struct edge *gtree;

  assert (numm > 0);
  startm = lastm + 1;
  lastm = lastm + numm;

  li = edgem->li;
  ei = edgem->edgeid;
  gtree = C[ci]->G[li].gtree;
  checkmig (lastm + 1, &gtree[ei].mig, &gtree[ei].cmm);

  for (i = startm; i <= lastm; i++)
    edgem->mig[i].mt = upt + uniform () * timein;
  edgem->mig[i].mt = -1;
  if (numm > 1)
    hpsortmig (&edgem->mig[startm] - 1, numm);
  lastpop = pop;
  if (constrainpop < 0)
  {
    for (i = startm; i <= lastm; i++)
    {
      edgem->mig[i].mp =
        picktopop (lastpop, C[ci]->plist[period], npops - period);
      lastpop = edgem->mig[i].mp;
    }
  }
  else if (numm >= 2)
  {
    i = startm;
    while (i < lastm - 1)
    {
      edgem->mig[i].mp =
        picktopop (lastpop, C[ci]->plist[period], npops - period);
      lastpop = edgem->mig[i].mp;
      i++;
    }
    edgem->mig[lastm - 1].mp =
      picktopop2 (lastpop, C[ci]->plist[period], npops - period,
                  constrainpop);
    edgem->mig[lastm].mp = constrainpop;
  }
  else
  {
    if (numm == 1)
      edgem->mig[lastm].mp = constrainpop;
  }
  return lastm;
}                               /* simmpath */



/********** FUNCTIONS DECLARED IN UPDATE_GTREE_COMMON.c***********/


/* calcmrate()
the count of migration events divided by the length of the relevant portion of the tree yields the 
migration rate on the tree, per unit time 

the migration rate prior is in units of migration events per mutation event 
to convert this to a relevant prior for the locus,  need to get the mutation rate scalar for this locus

8/08  removed use of mutation rate scalar in this
had thought it made sense for a given locus,  but now think not - essentially a bug
*/

/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre tree - e.g. help to avoid promoting
correlations and improve mixing  */
/* update_genealogy_covar()  still uses the overall migration rate for the tree */
#define MINMIG 1                //0.1
double
calcmrate (int mc, double mt)
{
  //double r;
  assert (mc >= 0);
  assert (mt > 0.0);

  if (mc == 0)
  {
    if (mt < 1)
      return 0.1;
    else
      return 0.1 / mt;
  }
  else
  {
    if (mt < 1)
      return 1.0;
    else
      return ((double) mc) / mt;
  }
}                               //calcmrate


/* SANGCHUL: Wed Dec 10 10:42:24 EST 2008
 * We could have used a simple memset function if members of structure
 * edgemiginfo are not dynamically allocated. We have to replace memset
 * function with function IMA_reset_edgemiginfo in order to reset structure
 * edgemiginfo.
 * member li must be initialized once. It should be not be changed in function
 * IMA_reset_edgemiginfo. Global edgemiginfo should have -1 value of li, which
 * is used in function checkmig. */
void
IMA_reset_edgemiginfo (struct edgemiginfo *em)
{
  em->edgeid = -1;
  em->sisid = -1;
  em->b = -1;
  em->e = -1;
  em->pop = -1;
  em->temppop = -1;
  em->fpop = -1;
  em->upt = -1.0;
  em->dnt = -1.0;
  assert (em->mtimeavail != NULL);
  memset (em->mtimeavail, 0, npops * sizeof (double));
  em->mtall = 0.0;              /* why this was -1.0 */
  assert (em->mp != NULL);
  memset (em->mp, 0, npops * sizeof (int));
  em->mpall = 0;                /* why this was -1 */
  assert (em->mig != NULL);
  em->mig[0].mt = -1.0;
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// LEVINE: P(D|G) computation
//
// internal node index of N: node - ngenes
// BitTrue (N[node - ngenes], i): add i to set element N[node - ngenes] 
static void
likelihoodDG_N (int ci, int li, UByteP * N)
{
  int ngenes;
  int parent;
  int i;
  struct edge *gtree;
  assert (N != NULL);
  gtree = C[ci]->G[li].gtree;
  ngenes = L[li].numgenes;
  for (i = 0; i < ngenes; i++)
  {
    parent = gtree[i].down;
    while (!(parent < 0))       // while parent is not root
    {
      BitTrue (N[parent - ngenes], i);
      parent = gtree[parent].down;
    }
  }
  return;
}

// up[0] is left, and up[1] is right child node.
// If a child node is tip, then we have a single-element set. Otherwise, we take
// the set element of the child node from set N.
static void
likelihoodDG_Nlr (int ci, int li, UByteP * Nl, UByteP * Nr, UByteP * N)
{
  int ngenes;
  int ngnodes;
  int upl;
  int upr;
  int i;
  struct edge *gtree;
  int size;
  assert (N != NULL);
  assert (Nl != NULL);
  assert (Nr != NULL);
  gtree = C[ci]->G[li].gtree;
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  for (i = ngenes; i < ngnodes; i++)
  {
    upl = gtree[i].up[0];
    if (gtree[upl].up[0] < 0)
    {
      BitTrue (Nl[i - ngenes], upl);
    }
    else
    {
      BitAssign (Nl[i - ngenes], N[upl - ngenes], size);
    }

    upr = gtree[i].up[1];
    if (gtree[upr].up[1] < 0)
    {
      BitTrue (Nr[i - ngenes], upr);
    }
    else
    {
      BitAssign (Nr[i - ngenes], N[upr - ngenes], size);
    }
  }
  return;
}

/* We assign Nl, and Nr to each Ui and Vi, respectively. */
/* We remove any tip that does not belong to population i. */
static void
likelihoodDG_UV (int ci, int li, UByteP * Pi, UByteP ** UV,
                 UByteP ** Ui, UByteP ** Vi, UByteP * Nl, UByteP * Nr)
{
  int i;
  int j;
  int k;
  int eu;
  int ev;
  int ngenes;
  int ngnodes;
  int node;
  struct edge *gtree;
  int size;
  assert (UV != NULL);
  assert (Ui != NULL);
  assert (Vi != NULL);
  assert (Nl != NULL);
  assert (Nr != NULL);
  gtree = C[ci]->G[li].gtree;
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  for (i = 0; i < npops; i++)
  {
    /* copy sets Nl and Nr to Ui and Vi. */
    for (j = ngenes; j < ngnodes; j++)
    {
      node = j - ngenes;
      /* remove all tips that do not belong to population i */
      BitIntersection (Ui[i][node], Nl[node], Pi[i], size, k);
      BitIntersection (Vi[i][node], Nr[node], Pi[i], size, k);
      /*
       * BitAssign (Ui[i][node], Nl[node], size, k);
       * BitAssign (Vi[i][node], Nr[node], size, k);
       * for (k = 0; k < ngenes; k++)
       *  {
       *    if (gtree[k].pop != i) if gene k does not belong to population i 
       *      {
       *        BitFalse(Ui[i][node], k);
       *        BitFalse(Vi[i][node], k);
       *      }
       *  }
       **/
      /* if either Ui or Vi is empty, then we remove that internal node. */
      BitEmpty (eu, Ui[i][node], size, k);
      BitEmpty (ev, Vi[i][node], size, k);
      if (eu == 0 || ev == 0)
      {
        /* Sets for i, node are removed from Ui and Vi. */
        BitFalse (UV[i][node], 0);
      }
      else
      {
        BitTrue (UV[i][node], 0);
      }
    }
  }
  return;
}

/* Ni = union of U and V */
static void
likelihoodDG_Ni (int ci, int li, UByteP ** Ni, UByteP ** UV,
                 UByteP ** Ui, UByteP ** Vi)
{
  int i;
  int j;
  int k;
  int ngenes;
  int ngnodes;
  int node;
  int size;
  assert (Ni != NULL);
  assert (Ui != NULL);
  assert (Vi != NULL);
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  for (i = 0; i < npops; i++)
  {
    /* take the union of Ui and Vi to have Ni */
    for (j = ngenes; j < ngnodes; j++)
    {
      node = j - ngenes;
      if (BitIsTrue (UV[i][node], 0))
      {
        BitUnion (Ni[i][node], Ui[i][node], Vi[i][node], size, k);
      }
    }
  }
  return;
}

/* Remove any element that is not paired with an element of the same set. */
/* We need to know which gene is paired with which other. */
static void
likelihoodDG_R (int ci, int li, UByteP ** Ri, UByteP ** Ni, UByteP ** UV)
{
  int i;
  int j;
  int k;
  int k_paired;
  int ngenes;
  int ngnodes;
  int node;
  int size;
  int *indices;
  assert (Ri != NULL);
  assert (Ni != NULL);
  indices = NULL;
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  for (i = 0; i < npops; i++)
  {
    for (j = ngenes; j < ngnodes; j++)
    {
      node = j - ngenes;
      /*BitAssign (Ri[i][node], Ni[i][node], size, k); */
      /* Is each element of Ri paired with an element of the Ri? */
      /* Is a pair of tips in set Ni? */
      /* for (each pair of population i) FIXME??? */
      /* or a list of elements in Ni */
      if (BitIsTrue (UV[i][node], 0))
      {
        for (k = 0; k < ngenes; k++)
        {
          if (BitIsTrue (Ni[i][node], k))
          {
            /* CHECK: element that is paired with that element */
            /*k_paired = diploid_genecopy (&gAGenes, li, k); */
            k_paired = L[li].pairs[k];
            if (BitIsTrue (Ni[i][node], k_paired))
            {
              BitTrue (Ri[i][node], k);
              BitTrue (Ri[i][node], k_paired);
            }
          }
        }
      }
    }
  }
  return;
}

/* One element in Ui and one in Vi are paired are fused to a single element, */
/* which is added to set Zi. */
static void
likelihoodDG_Z (int ci, int li, UByteP ** Zi, UByteP ** UV,
                UByteP ** Ui, UByteP ** Vi)
{
  int i;
  int j;
  int k;
  int k_paired;
  int ngenes;
  int ngnodes;
  int node;
  int size;
  int *indices;
  assert (Zi != NULL);
  assert (Ui != NULL);
  assert (Vi != NULL);
  indices = NULL;
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  for (i = 0; i < npops; i++)
  {
    for (j = ngenes; j < ngnodes; j++)
    {
      node = j - ngenes;
      if (BitIsTrue (UV[i][node], 0))
      {
        for (k = 0; k < ngenes; k++)    /* each element in Ui[i][node]) */
        {
          if (BitIsTrue (Ui[i][node], k))
          {
            /* CHECK: element that is paired with that element */
            /*k_paired = diploid_genecopy (&gAGenes, li, k); */
            k_paired = L[li].pairs[k];
            if (BitIsTrue (Vi[i][node], k_paired))
            {
              BitTrue (Zi[i][node], k);
              BitTrue (Zi[i][node], k_paired);
            }
          }
        }
      }
    }
  }
  return;
}

/* URi is the intersection of Ui and Ri */
static void
likelihoodDG_UR (int ci, int li, UByteP ** URi,
                 UByteP ** UV, UByteP ** Ui, UByteP ** Ri)
{
  int i;
  int j;
  int k;
  int ngenes;
  int ngnodes;
  int node;
  int size;
  assert (URi != NULL);
  assert (Ui != NULL);
  assert (Ri != NULL);
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  for (i = 0; i < npops; i++)
  {
    /* take the union of Ui and Vi to have Ni */
    for (j = ngenes; j < ngnodes; j++)
    {
      node = j - ngenes;
      if (BitIsTrue (UV[i][node], 0))
      {
        BitIntersection (URi[i][node], Ui[i][node], Ri[i][node], size, k);
      }
    }
  }
  return;
}

#ifndef IMA_math_lnfact_half
// If n is an integer, we call function of log n!.
// If n is represented as k + 1/2, we use the equality,
// \left(n+\frac{1}{2}\right)!=\sqrt{\pi}\prod_{k=0}^{n}\frac{2k+1}{2}
static double
IMA_math_lnfact_half (int n)
{
  double v;
  int quotient;
  int remainder;
  assert (!(n < 0));
  quotient = n / 2;
  remainder = n % 2;
  if (remainder == 0)
  {
    v = gsl_sf_lnfact (quotient);
  }
  else
  {
    v = M_LNPI / 2
      + gsl_sf_lnfact (2 * quotient + 2)
      - gsl_sf_lnfact (quotient + 1) - (2 * quotient + 2) * M_LN2;
  }
  return v;
}
#endif /* IMA_math_lnfact_half */

// Levine and Hypergeometric parts are computed.
// We return a logarithmic scaled value.
static double
likelihoodDG_compute (int n, int *ri, int *zi, int *ui, int *vi, int *uri)
{
  double v;
  double d;
  int i;
  int quotient;
  int remainder1;
  v = 0;
  for (i = 0; i < n; i++)
  {
    // Levine part
    v += zi[i] * M_LN2;
    IMA_math_lnfact_half (d, uri[i] + ri[i] - uri[i], quotient,
                          remainder1);
    v += d;
    v -= gsl_sf_lnfact (zi[i]);
    IMA_math_lnfact_half (d, uri[i] - zi[i], quotient, remainder1);
    v -= d;
    IMA_math_lnfact_half (d, ri[i] - uri[i] - zi[i], quotient,
                          remainder1);
    v -= d;
    v += gsl_sf_lnfact (uri[i]);
    v += gsl_sf_lnfact (ri[i] - uri[i]);
    v -= gsl_sf_lnfact (uri[i] + ri[i] - uri[i]);
    assert (ui[i] + vi[i] - ri[i] >= 0);
    assert (ui[i] - uri[i] >= 0);
    assert (vi[i] - ri[i] + uri[i] >= 0);
    assert (ri[i] >= 0);
    assert (uri[i] >= 0);
    assert (ri[i] - uri[i] >= 0);
    assert (ui[i] + vi[i] >= 0);
    assert (ui[i] >= 0);
    assert (vi[i] >= 0);
    // Hypergeometric part
    v += gsl_sf_lnfact (ui[i] + vi[i] - ri[i]);
    v -= gsl_sf_lnfact (ui[i] - uri[i]);
    v -= gsl_sf_lnfact (vi[i] - ri[i] + uri[i]);
    v += gsl_sf_lnfact (ri[i]);
    v -= gsl_sf_lnfact (uri[i]);
    v -= gsl_sf_lnfact (ri[i] - uri[i]);
    v -= gsl_sf_lnfact (ui[i] + vi[i]);
    v += gsl_sf_lnfact (ui[i]);
    v += gsl_sf_lnfact (vi[i]);
  }

  return v;
}

// We count numbers of elements of 5 sets: R, Z, U, V, and UR.
static double
likelihoodDG_count (int ci, int li, UByteP ** Ri,
                    UByteP ** Zi, UByteP ** UV,
                    UByteP ** Ui, UByteP ** Vi, UByteP ** URi)
{
  double v;
  int c;
  int *ri;
  int *zi;
  int *ui;
  int *vi;
  int *uri;
  int i;
  int j;
  int k;
  int l;
  int size;
  int nset;
  int ngenes;
  int ngnodes;
  int node;
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  ri = NULL;
  zi = NULL;
  ui = NULL;
  vi = NULL;
  uri = NULL;
  v = 0.0;
  for (i = 0; i < npops; i++)
  {
    nset = 0;
    for (j = ngenes; j < ngnodes; j++)
    {
      node = j - ngenes;
      if (BitIsTrue (UV[i][node], 0))
      {
        nset++;
      }
    }

    // All counts will be saved to ri, zi, ui, vi, and uri.
    ri = (int *) malloc (nset * sizeof (int));
    zi = (int *) malloc (nset * sizeof (int));
    ui = (int *) malloc (nset * sizeof (int));
    vi = (int *) malloc (nset * sizeof (int));
    uri = (int *) malloc (nset * sizeof (int));
    k = 0;
    for (j = ngenes; j < ngnodes; j++)
    {
      node = j - ngenes;
      if (BitIsTrue (UV[i][node], 0))
      {
        BitNumberTrue (c, Ri[i][node], size, l);
        ri[k] = c;
        BitNumberTrue (c, Zi[i][node], size, l);
        assert (c % 2 == 0);
        zi[k] = c / 2;
        BitNumberTrue (c, Ui[i][node], size, l);
        ui[k] = c;
        BitNumberTrue (c, Vi[i][node], size, l);
        vi[k] = c;
        BitNumberTrue (c, URi[i][node], size, l);
        uri[k] = c;
        k++;
      }
    }

#ifdef COMMENT
    for (j = 0; j < nset; j++)
    {
      fprintf (stdout, "Pops %d-[%d]: %d\t%d\t%d\t%d\t%d\n", i, j,
               ri[j], zi[j], ui[j], vi[j], uri[j]);
    }
#endif /* COMMENT */
    // Levine and Hypergeometric calculation
    v += likelihoodDG_compute (nset, ri, zi, ui, vi, uri);
    XFREE (ri);
    XFREE (zi);
    XFREE (ui);
    XFREE (vi);
    XFREE (uri);
  }

  return v;
}

#ifdef COMMENT
static void
print_powerset (const char *s, UByteP * A, int size, int n)
{
  int i;
  int j;
  for (i = 0; i < n; i++)
  {
    fprintf (stdout, "Set [%s] %d: ", s, i);
    for (j = 0; j < size * 8; j++)
    {
      if ((BitIsTrue (A[i], j)))
      {
        fprintf (stdout, "%d ", j);
      }
    }
    fprintf (stdout, "\n");
  }
  return;
}

static void
print_setpowerset (const char *s, UByteP ** A, int size, int m, int n)
{
  int i;
  int j;
  int k;
  for (i = 0; i < m; i++)
  {
    fprintf (stdout, "Set [%s] %d: ", s, i);
    for (j = 0; j < n; j++)
    {
      fprintf (stdout, "<<< [%d]: ", j);
      for (k = 0; k < size * 8; k++)
      {
        if ((BitIsTrue (A[i][j], k)))
        {
          fprintf (stdout, "%d ", k);
        }
      }
      fprintf (stdout, "[%d] >>> ", j);
    }
    fprintf (stdout, "\n");
  }
  return;
}
#endif /* COMMENT */

// Please, look up Jody's poplation assignment notes and find section Levine in
// SANGCHUL's note of imamp.nw.
// 
// We compute likelihood of diploid information given a genealogy. We have three
// basic sets, N, Nl, and Nr. We have 6 sets for each population: Ui, Vi, Ni,
// Ri, Zi, and URi. There are K contemporary populations, and each population
// has ng number of genes. We ignore the case where we have genes that are 
// paired with no other gene. We may need to consider the case. It may be okay 
// to ignore any tip nodes that are not paired with any gene. 
double
likelihoodDG (int ci, int li)
{
  double v;
  UByteP *Pi = NULL;
  UByteP *N = NULL;
  UByteP *Nl = NULL;
  UByteP *Nr = NULL;
  UByteP **UV = NULL;
  UByteP **Ui = NULL;
  UByteP **Vi = NULL;
  UByteP **Ni = NULL;
  UByteP **Ri = NULL;
  UByteP **Zi = NULL;
  UByteP **URi = NULL;
  int ngenes;
  int ngnodes;
  int size;
  int i;
  int j;
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  Pi = C[ci]->G[li].Pi;
  N = C[ci]->G[li].N;
  Nl = C[ci]->G[li].Nl;
  Nr = C[ci]->G[li].Nr;
  UV = C[ci]->G[li].UV;
  Ui = C[ci]->G[li].Ui;
  Vi = C[ci]->G[li].Vi;
  Ni = C[ci]->G[li].Ni;
  Ri = C[ci]->G[li].Ri;
  Zi = C[ci]->G[li].Zi;
  URi = C[ci]->G[li].URi;
  /* memory allocation */
  /* 
   * BitPowerSetNew(N, size, ngnodes, i);
   * BitPowerSetNew(Nl, size, ngnodes, i);
   * BitPowerSetNew(Nr, size, ngnodes, i);
   * BitSetPowerSetNew(UV, size, npops, ngnodes, i, j);
   * BitSetPowerSetNew(Ui, size, npops, ngnodes, i, j);
   * BitSetPowerSetNew(Vi, size, npops, ngnodes, i, j);
   * BitSetPowerSetNew(Ni, size, npops, ngnodes, i, j);
   * BitSetPowerSetNew(Ri, size, npops, ngnodes, i, j);
   * BitSetPowerSetNew(Zi, size, npops, ngnodes, i, j);
   * BitSetPowerSetNew(URi, size, npops, ngnodes, i, j); 
   **/
  BitPowerSetZero (N, size, ngnodes, i);
  BitPowerSetZero (Nl, size, ngnodes, i);
  BitPowerSetZero (Nr, size, ngnodes, i);
  BitSetPowerSetZero (UV, size, npops, ngnodes, i, j);
  BitSetPowerSetZero (Ui, size, npops, ngnodes, i, j);
  BitSetPowerSetZero (Vi, size, npops, ngnodes, i, j);
  BitSetPowerSetZero (Ni, size, npops, ngnodes, i, j);
  BitSetPowerSetZero (Ri, size, npops, ngnodes, i, j);
  BitSetPowerSetZero (Zi, size, npops, ngnodes, i, j);
  BitSetPowerSetZero (URi, size, npops, ngnodes, i, j);
  // sets N, Nl, Nr, Ui, Vi, Ni, Ri, Zi, and URi
  likelihoodDG_N (ci, li, N);
  likelihoodDG_Nlr (ci, li, Nl, Nr, N);
  likelihoodDG_UV (ci, li, Pi, UV, Ui, Vi, Nl, Nr);
  likelihoodDG_Ni (ci, li, Ni, UV, Ui, Vi);
  likelihoodDG_R (ci, li, Ri, Ni, UV);
  likelihoodDG_Z (ci, li, Zi, UV, Ui, Vi);
  likelihoodDG_UR (ci, li, URi, UV, Ui, Ri);
#ifdef COMMENT
  print_powerset ("N", N, size, ngnodes);
  print_powerset ("Nl", Nl, size, ngnodes);
  print_powerset ("Nr", Nr, size, ngnodes);
  print_setpowerset ("UV", UV, size, npops, ngnodes);
  print_setpowerset ("Ui", Ui, size, npops, ngnodes);
  print_setpowerset ("Vi", Vi, size, npops, ngnodes);
  print_setpowerset ("Ni", Ni, size, npops, ngnodes);
  print_setpowerset ("Ri", Ri, size, npops, ngnodes);
  print_setpowerset ("Zi", Zi, size, npops, ngnodes);
  print_setpowerset ("URi", URi, size, npops, ngnodes);
#endif /* COMMENT */
  // compute Levine and Hypergeometric equations
  v = likelihoodDG_count (ci, li, Ri, Zi, UV, Ui, Vi, URi);
  /* memory deallocation
   * BitSetPowerSetDelete(URi, npops, ngnodes, i, j);
   * BitSetPowerSetDelete(Zi, npops, ngnodes, i, j);
   * BitSetPowerSetDelete(Ri, npops, ngnodes, i, j);
   * BitSetPowerSetDelete(Ni, npops, ngnodes, i, j);
   * BitSetPowerSetDelete(Vi, npops, ngnodes, i, j);
   * BitSetPowerSetDelete(Ui, npops, ngnodes, i, j);
   * BitSetPowerSetDelete(UV, npops, ngnodes, i, j);
   * BitPowerSetDelete(Nr, ngnodes, i);
   * BitPowerSetDelete(Nl, ngnodes, i);
   * BitPowerSetDelete(N, ngnodes, i);
   **/
  return v;
}

void
likelihoodDG_init (int ci, int li)
{
  int ngenes;
  int ngnodes;
  int size;
  int i;
  int j;
  struct edge *gtree = NULL;
  ngenes = L[li].numgenes;
  size = ngenes / 8 + 1;
  ngnodes = 2 * ngenes - 1;
  gtree = C[ci]->G[li].gtree;
  /* memory allocation */
  BitPowerSetNew (C[ci]->G[li].Pi, size, npops, i);
  BitPowerSetNew (C[ci]->G[li].N, size, ngnodes, i);
  BitPowerSetNew (C[ci]->G[li].Nl, size, ngnodes, i);
  BitPowerSetNew (C[ci]->G[li].Nr, size, ngnodes, i);
  BitSetPowerSetNew (C[ci]->G[li].UV, size, npops, ngnodes, i, j);
  BitSetPowerSetNew (C[ci]->G[li].Ui, size, npops, ngnodes, i, j);
  BitSetPowerSetNew (C[ci]->G[li].Vi, size, npops, ngnodes, i, j);
  BitSetPowerSetNew (C[ci]->G[li].Ni, size, npops, ngnodes, i, j);
  BitSetPowerSetNew (C[ci]->G[li].Ri, size, npops, ngnodes, i, j);
  BitSetPowerSetNew (C[ci]->G[li].Zi, size, npops, ngnodes, i, j);
  BitSetPowerSetNew (C[ci]->G[li].URi, size, npops, ngnodes, i, j);
  for (i = 0; i < npops; i++)
  {
    for (j = 0; j < ngenes; j++)
    {
      if (gtree[j].pop == i)    /* if gene j does not belong to population i */
      {
        BitTrue (C[ci]->G[li].Pi[i], j);
      }
    }
  }

  return;
}

void
likelihoodDG_fin ()
{
  int ngenes;
  int ngnodes;
  int size;
  int ci;
  int li;
  int i;
  int j;
  for (ci = 0; ci < numchains; ci++)
  {
    for (li = 0; li < nloci; li++)
    {
      ngenes = L[li].numgenes;
      size = ngenes / 8 + 1;
      ngnodes = 2 * ngenes - 1;
      /* memory deallocation */
      BitSetPowerSetDelete (C[ci]->G[li].URi, npops, ngnodes, i, j);
      BitSetPowerSetDelete (C[ci]->G[li].Zi, npops, ngnodes, i, j);
      BitSetPowerSetDelete (C[ci]->G[li].Ri, npops, ngnodes, i, j);
      BitSetPowerSetDelete (C[ci]->G[li].Ni, npops, ngnodes, i, j);
      BitSetPowerSetDelete (C[ci]->G[li].Vi, npops, ngnodes, i, j);
      BitSetPowerSetDelete (C[ci]->G[li].Ui, npops, ngnodes, i, j);
      BitSetPowerSetDelete (C[ci]->G[li].UV, npops, ngnodes, i, j);
      BitPowerSetDelete (C[ci]->G[li].Nr, ngnodes, i);
      BitPowerSetDelete (C[ci]->G[li].Nl, ngnodes, i);
      BitPowerSetDelete (C[ci]->G[li].N, ngnodes, i);
      BitPowerSetDelete (C[ci]->G[li].Pi, npops, i);
    }
  }
  return;
}


void
init_gtreecommon (void)         // initialize copy edge
{
  int i, li, largestsamp;
  copyedge = malloc (3 * (sizeof (struct edge)));
  for (i = 0; i < 3; i++)
    copyedge[i].mig = malloc (ABSMIGMAX * sizeof (struct migstruct));   // very wasteful of space - should use dynamic memory and checkmig()
  if (somestepwise)
    for (i = 0; i < 3; i++)
    {
      copyedge[i].A = calloc (MAXLINKED, sizeof (int));
      copyedge[i].dlikeA = calloc (MAXLINKED, sizeof (double));
    }
  for (i = 0, li = 0; li < nloci; li++)
    if (i < L[li].numgenes)
      i = L[li].numgenes;
  largestsamp = i;
  nnminus1 = malloc ((largestsamp + 1) * sizeof (double));
  nnminus1[0] = 0;
  for (i = 1; i <= largestsamp; i++)
  {
    nnminus1[i] = (double) (i) * ((double) i - 1);
  }
  IMA_initmemory_edgemiginfo ();
}                               /* init_gtreecommon */


void
free_gtreecommon (void)
{
  int i;
  for (i = 0; i < 3; i++)
  {
    XFREE (copyedge[i].mig);
    if (somestepwise)
    {
      XFREE (copyedge[i].A);
      XFREE (copyedge[i].dlikeA);
    }
  }
  XFREE (copyedge);
  XFREE (nnminus1);
  return;
}                               //free_gtreecommon

/* We remove warning: shadow of global variable holdgtree. */
void
init_holdgtree (struct genealogy *G, int numgenes)
{
  int i;
  int numlines = 2 * numgenes - 1;
  G->gtree = calloc ((size_t) numlines, (sizeof (struct edge)));
  for (i = 0; i < numlines; i++)
  {
    G->gtree[i].mig = malloc (MIGINC * sizeof (struct migstruct));
    G->gtree[i].mig[0].mt = -1;
    G->gtree[i].cmm = MIGINC;
    G->gtree[i].up[0] = -1;
    G->gtree[i].up[1] = -1;
    G->gtree[i].down = -1;
    G->gtree[i].time = 0;
    G->gtree[i].mut = -1;
    G->gtree[i].pop = -1;
    if (somestepwise)
    {
      G->gtree[i].A = calloc (MAXLINKED, sizeof (int));
      G->gtree[i].dlikeA = calloc (MAXLINKED, sizeof (double));
    }
  }
  init_genealogy_weights (&(G->gweight));
}                               //init_holdgtree

void
free_holdgtree (struct genealogy *G, int numgenes)
{
  int i, numlines = 2 * numgenes - 1;
  for (i = 0; i < numlines; i++)
  {
    XFREE (G->gtree[i].mig);
    if (somestepwise)
    {
      XFREE (G->gtree[i].A);
      XFREE (G->gtree[i].dlikeA);
    }
  }
  free_genealogy_weights (&G->gweight);
  XFREE (G->gtree);
}                               //free_holdgtree

/* FIXME: DELETE VARIABLES: lmedgedrop, medgedrop!? */
void
storeoldedges (int ci, int li, int edge, int sisedge, int downedge)
{
  int i;
  double uptime;
  struct edge *gtree = C[ci]->G[li].gtree;
  copyedge[0].down = gtree[edge].down;
  i = -1;

  do
  {
    i++;
    if (!(ABSMIGMAX > i))
    {
      IM_err (IMERR_TOOMANYMIG, "step %d: locus [%d] edge [%d] mig %d > %d",
              step, li, edge, ABSMIGMAX, i);
    }
    copyedge[0].mig[i] = gtree[edge].mig[i];
  } while (copyedge[0].mig[i].mt > -0.5);
  copyedge[0].cmm = gtree[edge].cmm;
  medgedrop = i;
  if (edge < L[li].numgenes)
  {
    uptime = 0;
  }
  else
  {
    uptime = gtree[gtree[edge].up[0]].time;
  }
  if (uptime < C[ci]->tvals[lastperiodnumber - 1])
  {
    if (gtree[edge].time < C[ci]->tvals[lastperiodnumber - 1])
    {
      lmedgedrop = gtree[edge].time - uptime;
    }
    else
    {
      lmedgedrop = C[ci]->tvals[lastperiodnumber - 1] - uptime;
    }
  }
  else
  {
    lmedgedrop = 0;
    assert (medgedrop == 0);
  }
  copyedge[0].time = gtree[edge].time;
  copyedge[0].pop = gtree[edge].pop;
  copyedge[1].down = gtree[sisedge].down;
  i = -1;
  do
  {
    i++;
    if (!(ABSMIGMAX > i))
    {
      IM_err (IMERR_TOOMANYMIG, "step %d: locus [%d] edge [%d] mig %d > %d",
              step, li, sisedge, ABSMIGMAX, i);
    }
    copyedge[1].mig[i] = gtree[sisedge].mig[i];
  } while (copyedge[1].mig[i].mt > -0.5);
  copyedge[1].cmm = gtree[sisedge].cmm;
  copyedge[1].time = gtree[sisedge].time;
  copyedge[1].pop = gtree[sisedge].pop;
  copyedge[2].down = gtree[downedge].down;
  i = -1;
  do
  {
    i++;
    if (!(ABSMIGMAX > i))
    {
      IM_err (IMERR_TOOMANYMIG, "step %d: locus [%d] edge [%d] mig %d > %d",
              step, li, downedge, ABSMIGMAX, i);
    }
    copyedge[2].mig[i] = gtree[downedge].mig[i];
  } while (copyedge[2].mig[i].mt > -0.5);
  copyedge[2].cmm = gtree[downedge].cmm;
  copyedge[2].time = gtree[downedge].time;
  copyedge[2].pop = gtree[downedge].pop;
  copyedge[0].up[0] = gtree[edge].up[0];
  copyedge[0].up[1] = gtree[edge].up[1];
  copyedge[1].up[0] = gtree[sisedge].up[0];
  copyedge[1].up[1] = gtree[sisedge].up[1];
  copyedge[2].up[0] = gtree[downedge].up[0];
  copyedge[2].up[1] = gtree[downedge].up[1];
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    storeAinfo (ci, li, gtree, edge, sisedge, downedge);
}                               /* storoldeges */


/* set the gtree back to the way it was */
void
restoreedges (int ci, int li, int edge, int sisedge, int downedge,
              int newsisedge)
/*all this can be optimized some*/
{
  int i, j, ai, down;
  struct edge *gtree = C[ci]->G[li].gtree;
  if (newsisedge != sisedge)
  {
    down = gtree[downedge].down;
    if (down != -1)
    {
      if (gtree[down].up[0] == downedge)
        gtree[down].up[0] = newsisedge;
      else
        gtree[down].up[1] = newsisedge;
    }
    else
    {
      C[ci]->G[li].root = newsisedge;
      C[ci]->G[li].roottime = gtree[gtree[newsisedge].up[0]].time;
      assert (C[ci]->G[li].roottime <= TIMEMAX);
    }
    gtree[newsisedge].down = down;
    if (down != -1)
    {
      i = 0;
      while (gtree[newsisedge].mig[i].mt > -0.5)
        i++;
      j = -1;

      do
      {
        j++;
        checkmig (i + j, &(gtree[newsisedge].mig), &(gtree[newsisedge].cmm));
        gtree[newsisedge].mig[i + j] = gtree[downedge].mig[j];
      } while (gtree[downedge].mig[j].mt > -0.5);
    }
    else
    {
      gtree[newsisedge].mig[0].mt = -1;
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
          gtree[newsisedge].dlikeA[ai] = 0;
    }
    gtree[newsisedge].time = gtree[downedge].time;
  }
  gtree[edge].down = copyedge[0].down;
  i = -1;

  do
  {
    i++;
    checkmig (i, &(gtree[edge].mig), &(gtree[edge].cmm));
    gtree[edge].mig[i] = copyedge[0].mig[i];
  } while (gtree[edge].mig[i].mt > -0.5);
  gtree[edge].time = copyedge[0].time;
  gtree[edge].pop = copyedge[0].pop;
  down = gtree[sisedge].down;
  gtree[sisedge].down = copyedge[1].down;
  if (down != -1)
  {
    if (gtree[down].up[0] == sisedge)
      gtree[down].up[0] = downedge;

    else
      gtree[down].up[1] = downedge;
  }
  i = -1;

  do
  {
    i++;
    checkmig (i, &(gtree[sisedge].mig), &(gtree[sisedge].cmm));
    gtree[sisedge].mig[i] = copyedge[1].mig[i];
  } while (gtree[sisedge].mig[i].mt > -0.5);
  gtree[sisedge].time = copyedge[1].time;
  gtree[sisedge].pop = copyedge[1].pop;
  gtree[downedge].down = copyedge[2].down;
  i = -1;

  do
  {
    i++;
    checkmig (i, &(gtree[downedge].mig), &(gtree[downedge].cmm));
    gtree[downedge].mig[i] = copyedge[2].mig[i];
  } while (gtree[downedge].mig[i].mt > -0.5);
  gtree[downedge].time = copyedge[2].time;
  gtree[downedge].pop = copyedge[2].pop;
  gtree[downedge].up[0] = copyedge[2].up[0];
  gtree[downedge].up[1] = copyedge[2].up[1];
  if (gtree[downedge].down == -1)
  {
    C[ci]->G[li].roottime = gtree[gtree[downedge].up[0]].time;
    assert (C[ci]->G[li].roottime <= TIMEMAX);
    C[ci]->G[li].root = downedge;
  }
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
    {
      gtree[edge].A[ai] = copyedge[0].A[ai];
      gtree[sisedge].A[ai] = copyedge[1].A[ai];
      gtree[downedge].A[ai] = copyedge[2].A[ai];
      gtree[edge].dlikeA[ai] = copyedge[0].dlikeA[ai];
      gtree[sisedge].dlikeA[ai] = copyedge[1].dlikeA[ai];
      gtree[downedge].dlikeA[ai] = copyedge[2].dlikeA[ai];
      if (holdsisdlikeA[ai] != 0)
        gtree[newsisedge].dlikeA[ai] = holdsisdlikeA[ai];
    }
}                               /* restoreedges  */

// changes made 5/27/08  so that this could be called directly from update_genealogy_covar() 

/* info for edge goes in copyedge[0],
   info for sisedge goes in copyedge[1]
   info for downedge goes in copyedge[2] */
void
storeAinfo (int ci, int li, struct edge *gtree, int edge, int sisedge,
            int downedge)
{
  int ai;

  for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
  {
    copyedge[0].A[ai] = gtree[edge].A[ai];
    copyedge[0].dlikeA[ai] = gtree[edge].dlikeA[ai];
    copyedge[1].A[ai] = gtree[sisedge].A[ai];
    copyedge[1].dlikeA[ai] = gtree[sisedge].dlikeA[ai];
    copyedge[2].A[ai] = gtree[downedge].A[ai];
    //        if (copyedge[2].down != -1)  // commented this 5/27/08
    if (gtree[downedge].down != -1)   // added this 5/27/08 
    {
      copyedge[2].dlikeA[ai] = gtree[downedge].dlikeA[ai];
      //holddownA[ai] = gtree[copyedge[2].down].A[ai];
      holddownA[ai] = gtree[gtree[downedge].down].A[ai];      // changed this 5/27/08
    }
    else
    {
      copyedge[2].down = -1;  // inserted this 5/27/08
      holddownA[ai] = -1;
      copyedge[2].dlikeA[ai] = 0;
    }
  }
  return;
}

/* calculate the probability of the migration events summarised in edgem and sisem  */
/* closely follows mathematica notebook  'genealogy_updating_11_16_07.nb' */
double
getmprob (int ci, double mrate, struct edgemiginfo *edgem,
          struct edgemiginfo *sisem)
{
  double tempp, n, d;
  double r, logmrate, mt;
  int i, cm, pop, topop, popc, mp, pop1, pop2, cm1, cm2;

  tempp = 0.0;
  logmrate = log (mrate);
  if (sisem == NULL || sisem->mtall <= 0)       // only deal with edgem 
  {
    cm = 0;
    for (i = edgem->b; i <= edgem->e - 1; i++)
    {
      assert (edgem->mtimeavail[i] > 0.0);
      r = edgem->mtimeavail[i] * mrate * (npops - (i + 1));
      tempp += edgem->mp[i] * logmrate - r;
      cm += edgem->mp[i];
    }

    if (edgem->e < lastperiodnumber)
    {
      /* if assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1  do not enter here */
      if (edgem->e == lastperiodnumber - 1 
          && (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0 || npops == 2))     // only 2 pops in last period
      {
        assert (edgem->mtimeavail[edgem->e] > 0);
        if (ODD (edgem->mp[edgem->e]))
          tempp += edgem->mp[edgem->e] * logmrate - mylogsinh (mrate * edgem->mtimeavail[edgem->e]);
        else
          tempp += edgem->mp[edgem->e] * logmrate - mylogcosh (mrate * edgem->mtimeavail[edgem->e]);
      }
      else                      // 3 or more pops in the last period 
      {
        if (cm == 0)
          pop = edgem->pop;
        else
          pop = edgem->mig[cm - 1].mp;  // the last migration before this last period 
        while (C[ci]->poptree[pop].e <= edgem->e)
          pop = C[ci]->poptree[pop].down;
        topop = edgem->fpop;
        popc = (npops - edgem->e - 1);
        r = edgem->mtimeavail[edgem->e] * mrate * popc;
        if (pop == topop)
        {
          assert (edgem->mp[edgem->e] != 1);
          if (edgem->mp[edgem->e] == 0)
          {
            n = -r;
            d = log (1 - r * exp (-r));
          }
          else                  // >=2 
          {
            d =
              log ((1 - r * exp (-r)) * pathcondition (1,
                                                       edgem->
                                                       mp[edgem->e], popc));
            n = edgem->mp[edgem->e] * logmrate - r;
          }
        }
        else
        {
          assert (edgem->mp[edgem->e] > 0);
          if (edgem->mp[edgem->e] == 1)
          {
            d = log ((1 - exp (-r)) / popc);
            n = logmrate - r;
          }
          else                  // >= 2
          {
            d =
              log ((1 - exp (-r)) * pathcondition (0,
                                                   edgem->mp[edgem->
                                                             e], popc));
            n = edgem->mp[edgem->e] * logmrate - r;
          }
        }
        tempp += n - d;
      }
    }
  }
  else
  {
    if (edgem->mtall <= 0)      // only deal with sisem, code is mostly a copy of that above
    {
      cm = 0;
      for (i = sisem->b; i <= sisem->e - 1; i++)
      {
        r = sisem->mtimeavail[i] * mrate * (npops - (i + 1));
        tempp += sisem->mp[i] * logmrate - r;
        cm += sisem->mp[i];
      }
      if (sisem->e < lastperiodnumber)
      {
        if (sisem->e == lastperiodnumber - 1 && (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0 || npops == 2))   // only 2 pops in last period
        {
          if (ODD (sisem->mp[sisem->e]))
            tempp +=
              sisem->mp[sisem->e] * logmrate -
              mylogsinh (mrate * sisem->mtimeavail[sisem->e]);
          else
            tempp +=
              sisem->mp[sisem->e] * logmrate -
              mylogcosh (mrate * sisem->mtimeavail[sisem->e]);
        }
        else                    // 3 or more pops in the last period 
        {
          if (cm == 0)
            pop = sisem->pop;
          else
            pop = sisem->mig[cm - 1].mp;        // the last migration before this last period 
          while (C[ci]->poptree[pop].e <= sisem->e)
            pop = C[ci]->poptree[pop].down;
          topop = sisem->fpop;
          popc = (npops - sisem->e - 1);
          r = sisem->mtimeavail[sisem->e] * mrate * popc;
          if (pop == topop)
          {
            if (sisem->mp[sisem->e] == 0)
            {
              n = -r;
              d = log (1 - r * exp (-r));
            }
            else                // >=2 
            {
              d =
                log ((1 - r * exp (-r)) * pathcondition (1,
                                                         sisem->
                                                         mp[sisem->e], popc));
              n = sisem->mp[sisem->e] * logmrate - r;
            }
          }
          else
          {
            if (sisem->mp[sisem->e] == 1)
            {
              d = log ((1 - exp (-r)) / popc);
              n = logmrate - r;
            }
            else                // >= 2
            {
              d =
                log ((1 - exp (-r)) * pathcondition (0,
                                                     sisem->
                                                     mp[sisem->e], popc));
              n = sisem->mp[sisem->e] * logmrate - r;
            }
          }
          tempp += n - d;
        }
      }
    }
    else                        // have to deal with both edgem and sisem
    {
      cm1 = cm2 = 0;
      assert (edgem->e == sisem->e);
      for (i = IMIN (edgem->b, sisem->b); i <= edgem->e - 1; i++)
      {
        if (i >= edgem->b)
        {
          r = edgem->mtimeavail[i] * mrate * (npops - (i + 1));
          tempp += edgem->mp[i] * logmrate - r;
          cm1 += edgem->mp[i];
        }
        if (i >= sisem->b)
        {
          r = sisem->mtimeavail[i] * mrate * (npops - (i + 1));
          tempp += sisem->mp[i] * logmrate - r;
          cm2 += sisem->mp[i];
        }
      }
      if (edgem->e < lastperiodnumber)
      {
        if (edgem->e == lastperiodnumber - 1 && (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0 || npops == 2))   // only 2 pops in last period
        {
          mt = edgem->mtimeavail[edgem->e] + sisem->mtimeavail[sisem->e];
          mp = edgem->mp[edgem->e] + sisem->mp[sisem->e];
          if (ODD (mp))
            tempp += mp * logmrate - mylogsinh (mrate * mt);

          else
            tempp += mp * logmrate - mylogcosh (mrate * mt);
        }
        else                    //>=2 pops
        {
          if (cm1 == 0)
            pop1 = edgem->pop;
          else
            pop1 = edgem->mig[cm1 - 1].mp;      // the last migration before this last period 

          while (C[ci]->poptree[pop1].e <= edgem->e)
            pop1 = C[ci]->poptree[pop1].down;

          if (cm2 == 0)
            pop2 = sisem->pop;
          else
            pop2 = sisem->mig[cm2 - 1].mp;      // the last migration before this last period 

          while (C[ci]->poptree[pop2].e <= sisem->e)
            pop2 = C[ci]->poptree[pop2].down;

          mt = edgem->mtimeavail[edgem->e] + sisem->mtimeavail[sisem->e];
          mp = edgem->mp[edgem->e] + sisem->mp[sisem->e];
          popc = npops - edgem->e - 1;
          r = mt * mrate * popc;
          if (pop1 != pop2)
          {
            d = log ((1 - exp (-r)) * pathcondition (0, mp, popc));
            n = mp * logmrate - r;
          }
          else
          {
            if (mp == 0)
            {
              d = log (1 - r * exp (-r));
              n = -r;
            }
            else                //>= 2
            {
              assert (mp != 1);
              d = log ((1 - r * exp (-r)) * pathcondition (1, mp, popc));
              n = mp * logmrate - r;
            }
          }
          tempp += n - d;
        }
      }
    }
  }
  assert (tempp > -1e200 && tempp < 1e200);
  return tempp;
}                               /* getmprob */

/* fillmiginfoperiods() */
/* called by fillmiginfo  and by addmigration.  Fills in info on time periods spanned by an edge*/
/* b and e are the time periods during with the beginning and ends of the edge occur, respectively */
void
fillmiginfoperiods (int ci, struct edgemiginfo *em)
{
  int i;
  em->b = 0;
  while (em->upt > C[ci]->tvals[em->b])
    em->b++;
  em->e = em->b;
  while (em->dnt > C[ci]->tvals[em->e])
    em->e++;
  if (em->e == em->b)
  {
    if (em->b == lastperiodnumber)
      em->mtimeavail[em->b] = 0;
    else
      em->mtimeavail[em->b] = em->dnt - em->upt;
  }
  else
  {
    em->mtimeavail[em->b] = C[ci]->tvals[em->b] - em->upt;
    if (em->e == lastperiodnumber)
      em->mtimeavail[em->e] = 0;
    else
      em->mtimeavail[em->e] = em->dnt - C[ci]->tvals[em->e - 1];
    for (i = em->b + 1; i < em->e; i++)
      em->mtimeavail[i] = C[ci]->tvals[i] - C[ci]->tvals[i - 1];
  }
  if (em->b < lastperiodnumber)
  {
    em->mtall = DMIN (C[ci]->tvals[lastperiodnumber - 1], em->dnt) - em->upt;
  }
  else
  {
    em->mtall = 0;
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    assert (em->b == 0 && em->e == 0);
  }
  return;
}                               /*fillmiginfoperiods */

/* fillmiginfo() 
for updating edges and migration it is necessary to gather the info about the original edge 
and the original sister edge into one place to make sure it is easily available */
void
fillmiginfo (int ci, int li, struct edge *gtree, int edge, int sisedge)
{
  int i, j;
  double uptime;

/* REMOVED!
  memset (&oldedgemig, 0, sizeof (struct edgemiginfo));
  memset (&oldsismig, 0, sizeof (struct edgemiginfo));
*/
  IMA_reset_edgemiginfo (&oldedgemig);
  IMA_reset_edgemiginfo (&oldsismig);

  oldedgemig.edgeid = edge;
  oldsismig.edgeid = sisedge;
  oldedgemig.li = li;
  oldsismig.li = li;
  if (edge < L[li].numgenes)
    uptime = 0;
  else
    uptime = gtree[gtree[edge].up[0]].time;
  oldedgemig.pop = gtree[edge].pop;
  oldedgemig.dnt = gtree[edge].time;
  oldedgemig.upt = uptime;
  oldedgemig.fpop = gtree[gtree[edge].down].pop;
  fillmiginfoperiods (ci, &oldedgemig);
  oldedgemig.mig[0].mt = -1;
  for (j = 0; j <= oldedgemig.b; j++)
    oldedgemig.mp[j] = 0;
  oldedgemig.mpall = 0;
  i = 0;
  j = oldedgemig.b;
  while (gtree[edge].mig[i].mt > -0.5)
  {
    while (gtree[edge].mig[i].mt > C[ci]->tvals[j])
    {
      j++;
      oldedgemig.mp[j] = 0;
    }
    oldedgemig.mp[j]++;
    oldedgemig.mpall++;
    oldedgemig.mig[i] = gtree[edge].mig[i];
    i++;
  }
  oldedgemig.mig[i].mt = -1;
  if (sisedge >= 0)
  {
    // assert (gtree[sisedge].down == C[ci]->G[li].root); not true if this is called from t updating 
    if (sisedge < L[li].numgenes)
      uptime = 0;
    else
      uptime = gtree[gtree[sisedge].up[0]].time;
    oldsismig.pop = gtree[sisedge].pop;
    oldsismig.dnt = gtree[sisedge].time;
    oldsismig.upt = uptime;
    oldsismig.fpop = gtree[gtree[sisedge].down].pop;
    fillmiginfoperiods (ci, &oldsismig);
    oldsismig.mig[0].mt = -1;
    for (j = 0; j <= oldsismig.b; j++)
      oldsismig.mp[j] = 0;
    oldsismig.mpall = 0;
    i = 0;
    j = oldsismig.b;
    while (gtree[sisedge].mig[i].mt > -0.5)
    {
      while (gtree[sisedge].mig[i].mt > C[ci]->tvals[j])
      {
        j++;
        oldsismig.mp[j] = 0;
      }
      oldsismig.mp[j]++;
      oldsismig.mpall++;
      oldsismig.mig[i] = gtree[sisedge].mig[i];
      i++;
    }
    oldsismig.mig[i].mt = -1;
  }
}                               /* fillmiginfo */


// copy the migration info in newedgemig and newsismig  to the genealogy
void
copynewmig_to_gtree (int ci, int li)
{
  int i, pop, downperiod;
  struct edge *gtree = C[ci]->G[li].gtree;
  i = 0;
  while (newedgemig.mig[i].mt > 0)
  {
    checkmig (i, 
              &(gtree[newedgemig.edgeid].mig), 
              &(gtree[newedgemig.edgeid].cmm));
    gtree[newedgemig.edgeid].mig[i] = newedgemig.mig[i];
    i++;
  }
  gtree[newedgemig.edgeid].mig[i].mt = -1;
  if (newsismig.edgeid >= 0)
  {
    i = 0;

    //pop = newsismig.pop;
    while (newsismig.mig[i].mt > 0)
    {
      checkmig (i, &(gtree[newsismig.edgeid].mig),
                &(gtree[newsismig.edgeid].cmm));
      gtree[newsismig.edgeid].mig[i] = newsismig.mig[i];
      i++;
    }
    gtree[newsismig.edgeid].mig[i].mt = -1;

    // set top of down edge - seems to be necessary when root moves 
    if (i > 0)
      pop = gtree[newsismig.edgeid].mig[i - 1].mp;
    else
      pop = gtree[newsismig.edgeid].pop;
    downperiod = findperiod (ci, gtree[newsismig.edgeid].time);
    while (C[ci]->poptree[pop].e <= downperiod && C[ci]->poptree[pop].e != -1)
      pop = C[ci]->poptree[pop].down;
    gtree[gtree[newsismig.edgeid].down].pop = pop;
  }
}                               /* copynewmig_to_gtree */


/* store a few basic tree statistics */
void
storetreestats (int ci, int li, int mode)
{
  static double holdlength, holdtlength;
  static double holdroottime;
  static int holdroot;
  static int holdmig;
  if (mode == 0)
  {
    holdlength = C[ci]->G[li].length;
    holdtlength = C[ci]->G[li].tlength;
    holdroottime = C[ci]->G[li].roottime;
    holdroot = C[ci]->G[li].root;
    holdmig = C[ci]->G[li].mignum;
  }
  else
  {
    C[ci]->G[li].length = holdlength;
    C[ci]->G[li].tlength = holdtlength;
    C[ci]->G[li].mignum = holdmig;
    C[ci]->G[li].roottime = holdroottime;
    C[ci]->G[li].root = holdroot;
  }
}                               /* storetreestats */


/* return the population to which migration happens, chosen at random from those in plist*/
int
picktopop (int nowpop, int plist[], int numpops)
{
  int topop;

  do
  {
    topop = randposint (numpops);
  } while (plist[topop] == nowpop);

  assert (topop < numpops);
  return plist[topop];
}                               // picktopop


/* return the population to which migration happens, chosen at random from those in plist*/
/* this is same as picktopop() except that this will specifically avoid population notother */
int
picktopop2 (int nowpop, int plist[], int numpops, int notother)
{
  int topop;
  do
  {
    topop = randposint (numpops);
  } while (plist[topop] == nowpop || plist[topop] == notother);

  assert (topop < numpops);
  return plist[topop];
}                               //picktopop2



/* adds migration for single edges,  updates temppop as needed */
int
mwork (int ci, struct edgemiginfo *edgem, int lastmigperiod, double mrate)
{
  int i, lastm, mpall;
  double r;
  double timestartperiod;
  if (lastmigperiod < edgem->b)
    return 0;
  lastm = -1;
  mpall = 0;
  timestartperiod = edgem->upt;
  for (i = edgem->b, edgem->mpall = 0; i <= lastmigperiod; i++)
  {
    if (i == lastperiodnumber)
      break;
    while (C[ci]->poptree[edgem->temppop].e <= i
           && C[ci]->poptree[edgem->temppop].e != -1)
      edgem->temppop = C[ci]->poptree[edgem->temppop].down;

    /* the lengths is mt[i],  the rate per instant is mrate, and there are npops-i populations to go to */
    r = edgem->mtimeavail[i] * mrate * (npops - (i + 1));
    if (i < edgem->e)           // can do migration to any available pop
    {
      if ((edgem->mp[i] = poisson (r, -1)) > 0)
        lastm =
          simmpath (ci, edgem, i, edgem->mp[i], lastm,
                    edgem->mtimeavail[i], timestartperiod, edgem->temppop,
                    -1);
    }
    else
    {
      if (npops - i == 2)
      {
        if (edgem->temppop == edgem->fpop)      //even
          edgem->mp[i] = poisson (r, 0);
        else                    // odd
          edgem->mp[i] = poisson (r, 1);
        if (edgem->mp[i] > 0)
          lastm =
            simmpath (ci, edgem, i, edgem->mp[i], lastm,
                      edgem->mtimeavail[i], timestartperiod,
                      edgem->temppop, -1);
      }
      else
      {
        if (edgem->temppop == edgem->fpop)      //cannot be just one migration event
          edgem->mp[i] = poisson (r, 3);
        else                    // cannot be zero migration events 
          edgem->mp[i] = poisson (r, 2);
        if (edgem->mp[i] > 0)
          lastm =
            simmpath (ci, edgem, i, edgem->mp[i], lastm,
                      edgem->mtimeavail[i], timestartperiod,
                      edgem->temppop, edgem->fpop);
      }
    }
    if (edgem->mp[i] > 0)
    {
      mpall += edgem->mp[i];
      edgem->temppop = edgem->mig[lastm].mp;
    }
    timestartperiod = C[ci]->tvals[i];
  }
  if (edgem->mtimeavail[i] > 0)
  {
    if (C[ci]->poptree[edgem->temppop].e == i)
      edgem->temppop = C[ci]->poptree[edgem->temppop].down;
  }
  return mpall;
}                               /* mwork */

/* pathcondition() use to calculate the probability of a migration path condition 
e.g. that the path must end in a certain population, when there are popsless1 available populations
and there are moves migrations to do 
if number of populations is large or number of moves is large the probabability is 1/(npops-1)
uses recursion */
double
pathcondition (int issame, int moves, int popsless1)
{
  if ((popsless1 == 2 && moves >= 12) || (popsless1 > 2 && moves >= 6))
    return 1.0 / (popsless1 + 1);
  if (issame)
  {
    if (moves > 0)
      return pathcondition (0, moves - 1, popsless1);
    else
      return 1;
  }
  else
  {
    if (moves > 1)
      return pathcondition (0, moves - 2, popsless1) / popsless1 +
        pathcondition (0, moves - 1, popsless1) * (1.0 - 1.0 / popsless1);
    else
    {
      if (moves == 1)
        return 1.0 / popsless1;

      else
        return 0;
    }
  }
}                               /*pathcondition */


/************ GLOBAL FUNCTIONS **************/



/* return the period in the population tree of a particular time point */
int
findperiod (int ci, double t)
{
  int k = 0;
  assert (t >= 0);
  while (k < lastperiodnumber && C[ci]->tvals[k] <= t)
    k++;
  return k;
}

/* return the population that the edge is in at time.  If time is a population splitting time, it returns
the population state in the period above that split*/
int
nowedgepop (int ci, struct edge *gtree, double ptime)
{
  int pop, j;

  pop = gtree->pop;
  j = 0;
  while (gtree->mig[j].mt > -1 && gtree->mig[j].mt < ptime)
  {
    pop = gtree->mig[j].mp;
    j++;
  }
  while (ptime > C[ci]->poptree[pop].time && pop != -1)
  {
    pop = C[ci]->poptree[pop].down;
  }
  assert (pop >= 0 && pop < numtreepops);
  return pop;
}

void
init_treeweight (void)
{
  elength = ESTARTLENGTH;
  e = malloc (elength * sizeof (struct gtreeevent));
  eci = malloc (elength * sizeof (unsigned long));
}                               // init_treeweight

void
free_treeweight (void)
{
  XFREE (e);
  XFREE (eci);
}                               // free_treeweight


/* updated this on 11/3/08  to avoid most of the mallocing that was being done
created static dynamic arrays  e and eci  and just realloc them if they are not big enough
Also dropped the use of memset - just not needed */

void
treeweight (int ci, int li)
/* events are migrations, m = 1, or coalescent, c = 0 or population splitting, t = -1
e contains information for events. event i is contained in e+i and includes
the time till the next event, time
the event that happens at that time cmt (0 or 1 or -1)
the population in which the event happens 
if that event is migration, the pop is the one from which the migrant is leaving
the number of individuals in each population, just before the event happens.

accumulates info in G->gweight  
information is accumulated in each period k,  even if a population spans multiple periods. 
The info on coalescent events will have to be summed over periods for those populations that 
span multiple periods.  will have to be summed for population s
*/
{
  int i, ii, j, jj, k, n[2 * MAXPOPS - 1], nsum, ncount;
  int ip, jp;
  int nowpop;
  double t, timeinterval, fmtemp, lasttime, sumtime = 0;
  int ec;
  double h2term, hlog;
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  struct genealogy_weights *gweight = &(G->gweight);
  int ng = L[li].numgenes;
  double lastsplitt, timeadd;

#ifdef IMPROFILE
  startclock (&clock_treeweight);
#endif
// assert(_CrtCheckMemory());
  ncount = ng - 1;              // there must be ng-1 coalescent events

  if (modeloptions[NOMIGRATION] == 0)
  {
    G->mignum = 0;
    for (i = 0; i < 2 * ng - 1; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > -0.5)
        j++;
      G->mignum += j;             // count migration events 
    }
    ncount += G->mignum;
  }


  /* add an 'event' for every population tree split time that is younger than the root of the genealogy */
  ncount += findperiod (ci, G->roottime);

  if (ncount + 2 >= elength)
  {
    e = realloc (e, (ncount + EADDLENGTH) * sizeof (struct gtreeevent));
    eci = realloc (eci, (ncount + EADDLENGTH) * sizeof (unsigned long));
    elength = ncount + EADDLENGTH;
  }

  //e = malloc (ncount * sizeof (struct gtreeevent));    //maybe can remove this for most times through this loop - only change if ncount is bigger than before
  //memset (e, 0, ncount * sizeof (struct gtreeevent));   // initialize to zero's quickly 

  ec = 0;
  for (i = 0; i < 2 * ng - 1; i++)
  {
    nowpop = gtree[i].pop;
    assert (nowpop < numtreepops);
    if (i < ng)
    {
      t = 0;
    }
    else
    {
      // every internal branch begins in a population by a coalescent event, we know the time and the pop
      t = gtree[gtree[i].up[0]].time;
      assert (i >= ng);
      (e + ec)->time = t;
      (e + ec)->cmt = 0;
      (e + ec)->periodi = findperiod (ci, t);
      /// !!! added this next line on 8/21/08   why was this missing  ???
      // shouldn't this not be needed,  shouldn't nowpop be correct before this point ?? 
      /*  2/16/09  fixed the bug in changet_NW  that made this section necessary 
         while (C[ci]->poptree[nowpop].e <= (e + ec)->periodi
         && C[ci]->poptree[nowpop].e != -1)
         {
         nowpop = C[ci]->poptree[nowpop].down;
         } */
      (e + ec)->pop = nowpop;
      ec += 1;
    }
    j = 0;
    while (gtree[i].mig[j].mt > -0.5)
    {
      //assert(gtree[i].mig[j].mt >= t); 
      t = gtree[i].mig[j].mt;
      (e + ec)->time = t;
      (e + ec)->periodi = findperiod (ci, t);

      /* could change population by passing into the next period */
      while (C[ci]->poptree[nowpop].e <= (e + ec)->periodi
             && C[ci]->poptree[nowpop].e != -1)
      {
        nowpop = C[ci]->poptree[nowpop].down;
      }
      assert (ISELEMENT (nowpop, C[0]->periodset[(e + ec)->periodi]));
      (e + ec)->pop = nowpop;
      (e + ec)->topop = gtree[i].mig[j].mp;
      assert (nowpop != (e + ec)->topop);
      nowpop = (e + ec)->topop;
      assert ((e + ec)->topop >= 0 && (e + ec)->topop < numtreepops);
      (e + ec)->cmt = 1;        // a migration event
      ec += 1;
      j++;
    }
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
  {
    k = findperiod (ci, G->roottime);
    for (i = 0; i < k; i++)
    {
      (e + ec)->time = C[ci]->tvals[i];
      (e + ec)->periodi = i;
      (e + ec)->cmt = -1;
      ec += 1;
    }
  }
  assert (ncount == ec); /* FIXME: island model is checked here? */
  // assert(_CrtCheckMemory());

#ifdef IMPROFILE
  startclock (&clock_indexx);
#endif
  //eci = malloc (ncount * sizeof (unsigned long));
  indexx ((unsigned long) ec, e - 1, eci - 1);  // sorts an index of locations in ec by time;-1 because NR routines use arrays that begin at 1
#ifdef IMPROFILE
  addclock (&clock_indexx);
#endif

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (i = 0; i < npops; i++)
      n[i] = 0;
    for (i = 0; i < ng; i++)
    {
      n[gtree[i].pop]++;
    }
  }
  else
  {
    for (i = 0; i < npops; i++)
      n[i] = L[li].samppop[i];  // use regular num indexing  (not sets)
  }

  for (i = npops; i < 2 * npops - 1; i++)
    n[i] = 0;
  nsum = L[li].numgenes;
  lasttime = 0;
  G->length = G->tlength = 0;   // does not really belong here, but this is an easy place to measure tlength, the total length of genealogy
  h2term = 1 / (2 * L[li].hval);
  lastsplitt = C[ci]->tvals[lastperiodnumber - 1];
  k = 0;                        //period
  for (j = 0; j < ec; j++)
  {
    i = eci[j] - 1;             // -1 because NR routines use arrays that begin at 1
    timeinterval = (e + i)->time - lasttime;
    assert (timeinterval >= 0);

    timeadd = nsum * timeinterval;
    G->length += timeadd;
    if ((e + i)->time < lastsplitt)
    {
      G->tlength += timeadd;
    }
    else
    {
      if (lasttime < lastsplitt)
        G->tlength += nsum * (lastsplitt - lasttime);
    }
    lasttime = (e + i)->time;
    for (ii = 0; ii < npops - k; ii++)
    {
      if (modeloptions[SINGLEPOPULATION] == 0)
      {
        ip = C[ci]->plist[k][ii];
      }
      else
      {
        ip = 0;
      }
      gweight->fc[k][ii] += nnminus1[n[ip]] * timeinterval * h2term;
      assert (gweight->fc[k][ii] < DBL_MAX);
      if (!modeloptions[NOMIGRATION] && k < lastperiodnumber)
      {
        fmtemp = n[ip] * timeinterval;
        for (jj = 0; jj < npops - k; jj++)
          if (jj != ii)
          {
            gweight->fm[k][ii][jj] += fmtemp;
          }
      }
    }
    switch ((e + i)->cmt)
    {
    case 0:
      if (modeloptions[SINGLEPOPULATION] == 0)
      {
        assert (ISELEMENT ((e + i)->pop, C[ci]->periodset[findperiod (ci, lasttime)]));
        ip = (e + i)->pop;
        assert (n[ip] > 1);
        ii = 0;
        // this while() may be slow, could have some kind of lookup 
        while (C[ci]->plist[k][ii] != ip)
          ii++;
      }
      else
      {
        ip = 0;
        ii = 0;
      }
      gweight->cc[k][ii]++;
      n[ip]--;
      nsum--;
      break;
    case 1:
      assert (modeloptions[SINGLEPOPULATION] == 0);
      ip = (e + i)->pop;
      jp = (e + i)->topop;
      ii = jj = 0;
      // this while() may be slow, could have some kind of lookup 
      while (C[ci]->plist[k][ii] != ip)
        ii++;
      assert (ii < npops - k);
      // this while() may be slow, could have some kind of lookup 
      while (C[ci]->plist[k][jj] != jp)
        jj++;
      assert (jj < npops - k);
      gweight->mc[k][ii][jj]++;
      assert (n[ip] > 0);
      n[ip]--;
      n[jp]++;
      break;
    case -1:
      assert (modeloptions[SINGLEPOPULATION] == 0);
      assert (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0);
      assert ((e + i)->time >= C[ci]->tvals[k]);
      assert (nsum > 1);
      assert (k < numsplittimes);
      sumtime = C[ci]->tvals[k];
      k++;
      n[C[ci]->addpop[k]] =
        n[C[ci]->droppops[k][0]] + n[C[ci]->droppops[k][1]];
      n[C[ci]->droppops[k][0]] = n[C[ci]->droppops[k][1]] = 0;
      break;
    }
  }
  // assert(fabs((e + i)->time -C[ci]->G[li].roottime) < 1e-8);
  hlog = log (L[li].hval);
  for (k = 0; k < numsplittimes + 1; k++)
    for (i = 0; i < npops - k; i++)
      gweight->hcc[k][i] += hlog * gweight->cc[k][i];
  
  assert (nsum == 1);
  //assert(fabs(G->length - templ) < 1e-10);
#ifdef IMPROFILE
  addclock (&clock_treeweight);
#endif
}                               /*treeweight */

#undef ESTARTLENGTH
#undef EADDLENGTH



/* integrate_tree_prob - determine integrated probability of the genealogy  */
/* this function does not know about loci,  all the information needed is contained in 
  *gweight,  regardless of whether it is for one genealogy (called from updategenealogy) or all genealogies (called
  from changet ) 
  regardless of whether it is all genealogies or a single genealogy the caculations results are put into C[ci]->allpcalc*/

void
integrate_tree_prob (int ci,
                     struct genealogy_weights *gweight,
                     struct genealogy_weights *holdgweight,
                     struct probcalc *pcalc,
                     struct probcalc *holdpcalc, double *holdt)
{
  double psum;
  int i, j;
  int c, holdc;
  double hc, f, holdf;
  int checkm;
#ifdef IMPROFILE
  startclock (&clock_integrate_tree_prob);
#endif
  /* check to see if there are any migration events where there should not be under the model */
  if (modeloptions[NOMIGRATION] == 0 && nomigrationchecklist.n > 0)
  {
    i = 0;
    checkm = 1;
    while (checkm && i < nomigrationchecklist.n)
    {
      checkm = gweight->mc[nomigrationchecklist.p[i]]
                          [nomigrationchecklist.r[i]]
                          [nomigrationchecklist.c[i]] == 0;
      i++;
    }
  }
  else
  {
    checkm = 1;
  }

  if (!checkm)
  {
    psum = -MYDBL_MAX;          // very large neg value to ensure rejection of update if migration events occur where they should not under the model
  // not entirely sure if this is the best way to reject updates 
  }
  else
  {
    psum = 0;
    for (i = 0; i < numpopsizeparams; i++)
    {
      c = holdc = 0;
      f = holdf = hc = 0.0;
      for (j = 0; j < itheta[i].wp.n; j++)
      {
        c += gweight->cc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        holdc += holdgweight->cc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        f += gweight->fc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        holdf += holdgweight->fc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        hc += gweight->hcc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
      }
      if (c == holdc && f == holdf)
      {
        pcalc->qintegrate[i] = holdpcalc->qintegrate[i];
      }
      else
      {
        pcalc->qintegrate[i] =
          integrate_coalescent_term (c, f, hc, itheta[i].pr.max,
                                     itheta[i].pr.min);
      }

      psum += pcalc->qintegrate[i];
      //     assert(psum > -1e200 && psum < 1e200);
    }
    if (modeloptions[SPLITTINGRATEPARAMETER])
    {
      c = npops;
      hc = (double) npops;
      f = holdf = 0.0;
      for (i = 0; i < numsplittimes; i++)
      {
        holdf += holdt[i] * (npops - i);
        f += C[ci]->tvals[i] * (npops - i);
      }
      if (f == holdf)
        pcalc->sintegrate = holdpcalc->sintegrate;
      else
        pcalc->sintegrate =
          integrate_splitrate_term (f, isplit.pr.max, isplit.pr.min);
      psum += pcalc->sintegrate;
      //   assert(psum > -1e200 && psum < 1e200);
    }

    if (!modeloptions[NOMIGRATION])
    {
      for (i = 0; i < nummigrateparams; i++)
      {
        c = holdc = 0;
        f = holdf = 0.0;
        for (j = 0; j < imig[i].wp.n; j++)
        {
          c += gweight->mc[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
          holdc +=
            holdgweight->mc[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].
                                                              wp.c[j]];
          f += gweight->fm[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
          holdf +=
            holdgweight->fm[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].
                                                              wp.c[j]];
        }
        if (c == holdc && f == holdf)
        {
          pcalc->mintegrate[i] = holdpcalc->mintegrate[i];
        }
        else
        {

          if (modeloptions[EXPOMIGRATIONPRIOR])
            pcalc->mintegrate[i] =
              integrate_migration_term_expo_prior (c, f, imig[i].pr.mean);
          else
            pcalc->mintegrate[i] =
              integrate_migration_term (c, f, imig[i].pr.max, imig[i].pr.min);
        }
        psum += pcalc->mintegrate[i];
        //     assert(psum > -1e200 && psum < 1e200);
      }
    }
    //assert(psum > -1e200 && psum < 1e200);
  }
  //assert(psum > -1e200 && psum < 1e200);
  pcalc->probg = psum;
#ifdef IMPROFILE
  addclock (&clock_integrate_tree_prob);
#endif
}                               /* integrate_tree_prob */

// same as integrate_tree_prob but does not have checks to see if values match current ones */
void
initialize_integrate_tree_prob (int ci,
                                struct genealogy_weights *gweight,
                                struct probcalc *pcalc)
{
  double psum;
  int i, j;
  int c;
  double hc, f;
  int checkm;
  /* check to see if there are any migration events where there should not be under the model */
  if (modeloptions[NOMIGRATION] == 0 && nomigrationchecklist.n > 0)
  {
    i = 0;
    checkm = 1;
    while (checkm && i < nomigrationchecklist.n)
    {
      checkm = gweight->mc[nomigrationchecklist.p[i]]
                          [nomigrationchecklist.r[i]]
                          [nomigrationchecklist.c[i]] == 0;
      i++;
    }
  }
  else
  {
    checkm = 1;
  }

  if (!checkm)
  {
    psum = -MYDBL_MAX;          // very large neg value to ensure rejection of update if migration events occur where they should not under the model
  }
  else
  {
    psum = 0;
    for (i = 0; i < numpopsizeparams; i++)
    {
      c = 0;
      f = hc = 0.0;
      for (j = 0; j < itheta[i].wp.n; j++)
      {
        c += gweight->cc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        f += gweight->fc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        hc += gweight->hcc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
      }
      pcalc->qintegrate[i] =
        integrate_coalescent_term (c, f, hc, itheta[i].pr.max,
                                   itheta[i].pr.min);
      psum += pcalc->qintegrate[i];
      assert (psum > -1e200 && psum < 1e200);
    }
    if (modeloptions[SPLITTINGRATEPARAMETER])
    {
      c = npops;
      f = 0.0;
      for (i = 0; i < numsplittimes; i++)
      {
        f += C[ci]->tvals[i] * (npops - i);
      }
      pcalc->sintegrate =
        integrate_splitrate_term (f, isplit.pr.max, isplit.pr.min);
      psum += pcalc->sintegrate;
    }
    if (!modeloptions[NOMIGRATION])
    {
      for (i = 0; i < nummigrateparams; i++)
      {
        c = 0;
        f = 0.0;
        for (j = 0; j < imig[i].wp.n; j++)
        {
          c += gweight->mc[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
          f += gweight->fm[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
        }
        if (modeloptions[EXPOMIGRATIONPRIOR])
        {
          pcalc->mintegrate[i] =
            integrate_migration_term_expo_prior (c, f, imig[i].pr.mean);
        }
        else
        {
          pcalc->mintegrate[i] =
            integrate_migration_term (c, f, imig[i].pr.max, imig[i].pr.min);
        }
        psum += pcalc->mintegrate[i];
      }
    }
    assert (psum > -1e200 && psum < 1e200);
  }
  pcalc->probg = psum;
}                               /* initialize_integrate_tree_prob */



/* for HKY model */
void
copyfraclike (int ci, int li)
{
  int i, j, k;
  struct edge *gtree = C[ci]->G[li].gtree;
  int ng = L[li].numgenes;
  for (i = ng; i < 2 * ng - 1; i++)
  {
    if (gtree[i].hkyi.newfrac[0][0] != -1)
    {
      for (j = 0; j < L[li].numsites; j++)
        for (k = 0; k < 4; k++)
          gtree[i].hkyi.frac[j][k] = gtree[i].hkyi.newfrac[j][k];
    }
  }
}

/* for HKY model */
void
storescalefactors (int ci, int li)
{
  int i, k;
  struct edge *gtree = C[ci]->G[li].gtree;
  for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
    for (k = 0; k < L[li].numsites; k++)
      gtree[i].hkyi.oldscalefactor[k] = gtree[i].hkyi.scalefactor[k];
}


/* for HKY model */
void
restorescalefactors (int ci, int li)
{
  int i, k;
  struct edge *gtree = C[ci]->G[li].gtree;
  for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
  {
    for (k = 0; k < L[li].numsites; k++)
    {
      gtree[i].hkyi.scalefactor[k] = gtree[i].hkyi.oldscalefactor[k];
    }
  }
}


/*  finishSWupdateA() returns the likelihood associated with the allele states on the parts of the gtree that have changed
   it returns the difference in likelihoods associated with these parts.  thus the overall likelihood of data can just be
   updated by the value returned from this function. 

	when  finishSWupdate() is entered the gtree is not complete, as the allele state of the new downedge is not known
    must pick an allele state for the new downedge:
    the new downedge has the same number as the old downedge, so it is just called downedge hereafter

    at the start the current allele value for downedge is oldA - this however comes from a different part 
       of the gtree from the old downedge that was erased
   
    consider - the point at the base of edge is connected to 3 other nodes (2 up and 1 down) unless
    edge is the root in which case there are only 2 ups
    
    start w/ a old value of A at the top of downedge = oldA

    determine the mean (weighted by the length of the branch to the node) of the difference between the A values
    of the connecting nodes and oldA 
    Treat this mean value as the parameter of a geometric distribution.
    pick an rv that is in the necessary direction (+ or -) and add this to oldA  this becomes newA. 
    
    the allele value for downedge becomes newA. the new genealogy is now complete.
    
    now consider the reverse update from new to old A

    again, calculate the weighted mean difference (for the old genealogy) between adjacent allele states and newA. 
    pick a geometric rv 

    Consider the update as this - When the gtree jumped from one genealogy to another, this caused oldA to be set to newA

        
    Aterm = log(Pr(oldstate in node | T* -> T)/Pr(new state in node | T -> T*)) is the ratio of two geometric probabilities. 
	Aterm is part of the Hastings term, the proposal ratio.
    This is passed to the calling function along with the new likelihood. */
double
finishSWupdateA (int ci, int li, int ai, int edge,
                 int downedge, int sisedge, int newsisedge,
                 double u, double *Aterm)
{
  int i, j, newA, oldA;
  int d, dA, e[3];
  double geotermnew, geotermold;
  double upt, time[3];
  int wsumdiff;
  double oldlikeadj = 0, likeadj = 0;
  struct edge *gtree = C[ci]->G[li].gtree;
  oldA = gtree[downedge].A[ai];
  if (newsisedge != sisedge)
    holdsisdlikeA[ai] = gtree[newsisedge].dlikeA[ai];
  else
    holdsisdlikeA[ai] = 0;

  /* NOTE: Function storeoldedge must be called to set copyedge. Sang Chul has 
   * used finishSWupdateA without calling function storeoldedge. */
  oldlikeadj =
    copyedge[0].dlikeA[ai] + copyedge[1].dlikeA[ai] +
    copyedge[2].dlikeA[ai] + holdsisdlikeA[ai];
  e[0] = edge;
  e[1] = newsisedge;
  e[2] = gtree[e[0]].down;
  for (i = 0; i < 3; i++)
  {
    if (e[i] >= L[li].numgenes)
    {
      upt = gtree[gtree[e[i]].up[0]].time;
      assert (upt == gtree[gtree[e[i]].up[1]].time);
    }
    else
    {
      upt = 0;
    }

    if (gtree[e[i]].down != -1)
    {
      time[i] = gtree[e[i]].time - upt;
      assert (time[i] >= 0.0);
    }
  }
  assert (gtree[e[0]].time == gtree[e[1]].time);
  wsumdiff = 0;
  for (i = 0, j = 0; i < 3; i++)
  {
    if (gtree[e[i]].down != -1)
    {
      if (i < 2)
        d = gtree[e[i]].A[ai] - oldA;
      else
        d = gtree[gtree[e[i]].down].A[ai] - oldA;
      wsumdiff += abs (d);
      j++;
    }
  }
  geotermnew = j / ((double) (wsumdiff + j));
  assert (geotermnew > 0);
  if (geotermnew > 0.95)
    geotermnew = 0.95;
  dA = geometric (geotermnew) - 1;
  if (bitran () /*uniform() < 0.5 */ )
    dA = -dA;
  if (dA >= 0)
    newA = IMIN (L[li].maxA[ai], oldA + dA);
  else
    newA = IMAX (L[li].minA[ai], oldA + dA);
  gtree[downedge].A[ai] = newA;
  assert (newA >= L[li].minA[ai]);
  dA = newA - oldA;             /* difference between new value, and what it would be based simply on weighted mean */
  if (gtree[sisedge].down != -1 && sisedge != newsisedge)
  {
    if (sisedge >= L[li].numgenes)
    {
      upt = gtree[gtree[sisedge].up[0]].time;
    }
    else
    {
      upt = 0;
    }

    /* set dlikeA for the old sisedge (which is now continuous with the old downedge */
    d = gtree[sisedge].A[ai] - gtree[gtree[sisedge].down].A[ai];
    assert (d <= L[li].maxA[ai]);
    gtree[sisedge].dlikeA[ai] = -(gtree[sisedge].time - upt) * u + log (bessi (d, (gtree[sisedge].time - upt) * u));
  }
  else
  {
    gtree[sisedge].dlikeA[ai] = 0;
  }

  likeadj = gtree[sisedge].dlikeA[ai];
  for (i = 0, j = 0; i < 3; i++)
  {
    if (gtree[e[i]].down != -1)
    {
      if (i < 2)
        d = gtree[e[i]].A[ai] - newA;
      else
        d = gtree[gtree[e[i]].down].A[ai] - newA;
      gtree[e[i]].dlikeA[ai] = -(time[i] * u) + log (bessi (d, time[i] * u));
      likeadj += gtree[e[i]].dlikeA[ai];
    }
  }
  wsumdiff = 0;
  j = 0;
  for (i = 0; i < 2; i++)
  {
    d = copyedge[i].A[ai] - newA;
    wsumdiff += abs (d);
    j++;
  }
  if (copyedge[2].down != -1)
  {
    d = holddownA[ai] - newA;
    wsumdiff += abs (d);
    j++;
  }
  geotermold = j / ((double) (wsumdiff + j));
  if (geotermold > 0.95)
    geotermold = 0.95;

  /* FIXME: Tue Jan 27 17:35:50 EST 2009
   * SANGCHUL thought that we could return 0 without computing anything. Now
   * that I think that we may have to compute parts of the above to propose a
   * new state of alleles at internal nodes. We may place this before computing
   * likelihood because we do not need compute likelihood when we run without
   * data. In the same reason, we may want to rearrange the code of function
   * likelihoodSW. */
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
  {
    *Aterm = 0;
    return 0;
  }
  else
  {
    *Aterm =
      (abs (dA) * log (1 - geotermold) + log (geotermold)) -
      (abs (dA) * log (1 - geotermnew) + log (geotermnew));
    return (likeadj - oldlikeadj);
  }
}                               /* finishSWupdateA */


double
updateA (int ci, int li, int ai, double u, int *count)
/* implements a geometric distribution based on the weighted mean difference in allele sizes surrounding the value in the 
node to be updated. branch lengths are used for the weighting */
/* would it be good to randomize the order in which nodes are updated ?  */
{
  int i, upA[2], downA, newA, oldA;
  int d;
  int up[2], upup[2], down;
  double pdg, like, dlikeup[2], dlikedown;
  double tup[2], tdown, weightsum, metropolishastingsterm;
  double wsumdiff, geotermnew, geotermold, Aterm;
  struct edge *gtree = C[ci]->G[li].gtree;
  int ng = L[li].numgenes;
  *count = 0;
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
  {
    *count = (int) (updateAfrac * (ng - 1));    // don't do anything as it can have no effect 
    return 0;
  }
  for (i = ng; i < 2 * ng - 1; i++)
  {
    /* don't bother with all nodes  just do uprop of them */
    if (uniform () < updateAfrac)
    {
      L[li].A_rec[ai].upinf->tries++;
      up[0] = gtree[i].up[0];
      up[1] = gtree[i].up[1];
      upup[0] = gtree[up[0]].up[0];
      upup[1] = gtree[up[1]].up[0];
      down = gtree[i].down;
      if (upup[0] == -1)
        tup[0] = gtree[up[0]].time;
      else
        tup[0] = gtree[up[0]].time - gtree[upup[0]].time;
      if (upup[1] == -1)
        tup[1] = gtree[up[1]].time;
      else
        tup[1] = gtree[up[1]].time - gtree[upup[1]].time;
      oldA = gtree[i].A[ai];
      upA[0] = gtree[gtree[i].up[0]].A[ai];
      upA[1] = gtree[gtree[i].up[1]].A[ai];
      pdg = gtree[gtree[i].up[0]].dlikeA[ai] + gtree[gtree[i].up[1]].dlikeA[ai];

      wsumdiff = abs (oldA - upA[0]) / tup[0] + abs (oldA - upA[1]) / tup[1];
      weightsum = 1 / tup[0] + 1 / tup[1];
      if (down != -1)
      {
        tdown = gtree[i].time - gtree[up[0]].time;
        downA = gtree[gtree[i].down].A[ai];
        wsumdiff += abs (oldA - downA) / tdown;
        weightsum += 1 / tdown;
        pdg += gtree[i].dlikeA[ai];
      }
      else
      {
        downA = -1;
      }
      geotermnew = weightsum / (wsumdiff + weightsum);
      assert (geotermnew > 0);
      if (geotermnew > 0.95)
        geotermnew = 0.95;
      d = geometric (geotermnew) - 1;
      assert (d >= 0);
      if (bitran () /*uniform() < 0.5 */ )
        newA = IMAX (L[li].minA[ai], oldA - d);
      else
        newA = IMIN (L[li].maxA[ai], oldA + d);

      if (newA != oldA)
      {
        wsumdiff = abs (newA - upA[0]) / tup[0] + abs (newA - upA[1]) / tup[1];
        if (down != -1)
        {
          wsumdiff += abs (newA - downA) / tdown;
        }
        geotermold = weightsum / (wsumdiff + weightsum);
        assert (geotermold > 0);
        if (geotermold > 0.95)
          geotermold = 0.95;
        d = newA - upA[0];
        like = dlikeup[0] = -(tup[0] * u) + log (bessi (d, tup[0] * u));
        d = newA - upA[1];
        like += dlikeup[1] = -(tup[1] * u) + log (bessi (d, tup[1] * u));
        if (down != -1)
        {
          d = newA - downA;
          like += dlikedown = -(tdown * u) + log (bessi (d, tdown * u));
        }
        Aterm = abs (newA - oldA) * log ((1 - geotermold) / (1 - geotermnew)) + log (geotermold / geotermnew);
        metropolishastingsterm = exp (beta[ci] * (like - pdg) + Aterm);
        if (metropolishastingsterm >= 1.0
            || metropolishastingsterm > uniform ())
        {
          gtree[i].A[ai] = newA;
          if (down != -1)
            gtree[i].dlikeA[ai] = dlikedown;
          gtree[up[0]].dlikeA[ai] = dlikeup[0];
          gtree[up[1]].dlikeA[ai] = dlikeup[1];
          //*count = *count + 1;
          *count = *count + 1;
          L[li].A_rec[ai].upinf->accp++;
        }
      }
    }
  }
  like = 0;
  for (i = 0; i < 2 * ng - 1; i++)
  {
    if (gtree[i].down != -1)
    {
      like += gtree[i].dlikeA[ai];
    }
  }
  return like;
}                               /*updateA */

/* getnewt() these window widths just seem to work,  no method to them */
/*whichupdate:
0 : NW
1:  RY1
2:  RY2
*/

#define MAXWIN 10
double
getnewt (int timeperiod, double t_u_prior, double t_d_prior, double oldt,
         int whichupdate)
{
  int i = 0;
  double twin, newt;
#ifndef EXPOTPROPOSE
  double U;
#endif
  if (timeperiod == (lastperiodnumber - 1) && modeloptions[SPLITTINGRATEPARAMETER])     // allow for wider range between t_d_prior and t_u_prior
  {
    switch (whichupdate)
    {
    case 0:
      if (npops == 2)
        twin = (t_d_prior - t_u_prior) / (25 * log ((double) nloci + 1.0) * npops);        // update rate just two low in splitting rate models with 2 pops, so raised the denominator 
      else
        twin = (t_d_prior - t_u_prior) / (12 * log ((double) nloci + 1.0) * npops);
      break;
    case 1:
      twin = (t_d_prior - t_u_prior) / (6 * log ((double) nloci + 1.0) * npops);
      break;
    case 2:
      twin = (t_d_prior - t_u_prior) / (20 * log ((double) nloci + 1.0) * npops);
      break;
    }
  }
  else
  {
    switch (whichupdate)
    {
    case 0:
      twin =
        (t_d_prior -
         t_u_prior) / (6 * log ((double) nloci + 1.0) * (npops - timeperiod));
      break;
    case 1:
      twin =
        (t_d_prior - t_u_prior) / (log ((double) nloci + 1.0) * (npops - timeperiod));
      break;
    case 2:
      twin =
        (t_d_prior -
         t_u_prior) / (10 * log ((double) nloci + 1.0) * (npops - timeperiod));
      break;
    }
  }

#ifdef EXPOTPROPOSE
  if (bitran ())
    newt = oldt - expo (1 / twin);
  else
    newt = oldt + expo (1 / twin);
#else // regular uniform distribution
  U = uniform ();
  newt = (oldt - twin / 2) + U * twin;
#endif

  /* this is a real kluge  1/15/09
     we want to have newt picked from an exponential and yet still fall between the t_d_prior and t_u_prior
     and also the probably of proposing newt given oldt equal to the probability of proposing oldt given newt 
     how best to handle values outside of the boundn of t_d_prior and t_u_prior ??  
     this folding back routine can be chaotic and fail, so we set a maximum on how many times through and then pick a new value


     For some reason it is not sufficient to just pick a new uniform value when newt is outsite of the bounds. Not sure why this is. 
   */
  while ((newt >= t_d_prior) || (newt <= t_u_prior))
  {
    if (newt >= t_d_prior)
      newt = 2.0 * t_d_prior - newt;
    else
      newt = 2.0 * t_u_prior - newt;
    i++;
    if (i == MAXWIN)
    {
      newt = t_u_prior + uniform () * (t_d_prior - t_u_prior);
#ifdef _DEBUG
      printf ("******* TRAPPED IN getnewt(), CALLED BY %d ********\n",
              whichupdate);
#endif
    }

  }
  assert (newt < t_d_prior && newt > t_u_prior);
  return newt;
}                               //getnewt()
