/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

/* 

This implements an update of multiple genealogy branches.  Each branch gets updated in a way similat ro that in update_gtree()
see  correllated_updates_notes_7_08.nb

tested with and without this update a fair bit in 8/08
does not seem to offer any big improvement over single branch updating, but keep it in
there may be some data sets that benefit.
*/

#undef GLOBVARS

#include "imamp.h"
#include "update_gtree_common.h"
#include "updateassignment.h"
/* #include "update_gtree.h"
 */

extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;

/* declarded in update_gtree.h
int holddownA[MAXLINKED];
int medgedrop;
int mrootdrop;
int rootmove;
double lmedgedrop;
double lmrootdrop;
double holdsisdlikeA[MAXLINKED];
struct genealogy_weights holdgweight_updategenealogy;
struct genealogy_weights holdallgweight_updategenealogy;
struct probcalc holdallpcalc_updategenealogy; 
struct genealogy holdgtree;
*/
extern int rootmove;            /* declared in update_gtree.c */
extern struct genealogy holdgtree;      /* declared in update_gtree_common.c */


// local to this file

static int *num_gtree_covar_updates;
static struct genealogy_weights holdallgweight;
static struct probcalc holdallpcalc;
static int *inpoplist;
static double *slidedist;
static int largestsamp;

static void copy_gtree (int li, struct genealogy *G, int mode);
static double get_popbranchlengths (int cpop, int *numnotroot,
                                    struct edge *gtree, struct genealogy *G,
                                    int ci, int li);
static void fillpicklists (int li, struct edge *gtree, struct genealogy *G,
                           int *incount, int cpop, int *inpoplist);
static void addmigration_covar (int ci, int li, int *oldmigcount,
                                double *oldtlength, double *form,
                                double *form_inv, double *backm,
                                double *backm_inv, double migmultiplier);
#define MAX_NUM_GTREE_COVAR_UPDATES 5
#define MIN_NUM_GTREE_COVAR_UPDATES 1
#define NUM_GTREE_COVAR_UPDATES_DIVISOR 3.0
#define GTREE_COVAR_MIGRATION_MULTIPLIER  3.0
#define GTREE_COVAR_SLIDEPARAM_MULTIPLIER  2.0
#define SLIDESTDVMAX 2          // set this shorter than is used in updategenealogy()


// copy_gtree  if mode==1  copy *L into holdgtree  if mode==0 copy holdgtree into *L
void
copy_gtree (int li, struct genealogy *G, int mode)
{
  struct edge *togtree, *fromgtree;
  struct genealogy *fromG, *toG;
  int i, j, ai;

  if (mode)
  {
    toG = &holdgtree;
    fromG = G;
  }
  else
  {
    toG = G;
    fromG = &holdgtree;
  }
  togtree = toG->gtree;
  fromgtree = fromG->gtree;
  toG->length = fromG->length;
  toG->mignum = fromG->mignum;
  toG->root = fromG->root;
  toG->roottime = fromG->roottime;
  toG->tlength = fromG->tlength;
  for (i = 0; i < L[li].numlines; i++)
  {
    togtree[i].down = fromgtree[i].down;
    j = -1;
    do
    {
      j++;
#ifdef ASSIGNMENTWMIG
      checkmig (j, &(togtree[i]));
#else
      checkmig (j, &(togtree[i].mig), &(togtree[i].cmm));
#endif /* ASSIGNMENTWMIG */
      togtree[i].mig[j] = fromgtree[i].mig[j];
    } while (fromgtree[i].mig[j].mt > -0.5);

    togtree[i].time = fromgtree[i].time;
    togtree[i].pop = fromgtree[i].pop;
    togtree[i].up[0] = fromgtree[i].up[0];
    togtree[i].up[1] = fromgtree[i].up[1];
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
      {
        togtree[i].A[ai] = fromgtree[i].A[ai];
        togtree[i].dlikeA[ai] = fromgtree[i].dlikeA[ai];
      }
  }
}                               /* copy_gtree */


double
get_popbranchlengths (int cpop, int *numnotroot, struct edge *gtree,
                      struct genealogy *G, int ci, int li)
{
  int i;
  double upt;
  double tempbranchlength = 0;
  int popcount = 0, rootcount = 0;
  int period = C[ci]->poptree[cpop].b;

  if (cpop < npops)
  {
    popcount = L[li].samppop[cpop];
    for (i = 0; i < L[li].numgenes; i++)
      if (gtree[i].pop == cpop)
        tempbranchlength += gtree[i].time;
  }
  for (i = L[li].numgenes; i < L[li].numlines; i++)
    if (i != G->root && gtree[i].pop == cpop)
    {
      popcount++;
      upt =
        (period ==
         0) ? gtree[gtree[i].up[0]].time : DMAX (gtree[gtree[i].up[0]].time,
                                                 C[ci]->tvals[period - 1]);
      assert (gtree[i].time > upt);
      tempbranchlength += gtree[i].time - upt;
      rootcount += (gtree[i].pop == C[ci]->rootpop);
    }
  tempbranchlength /= (double) popcount;
  *numnotroot = L[li].numlines - rootcount - 1;
  return tempbranchlength;
}                               /* get_popbranchlength */

// fill lists of node numbers
void
fillpicklists (int li, struct edge *gtree, struct genealogy *G, int *incount,
               int cpop, int *inpoplist)
{
  int i;
  for (i = 0, *incount = 0; i < L[li].numlines; i++)
  {
    if (i != G->root)
    {
      if (gtree[i].pop == cpop)
      {
        inpoplist[*incount] = i;
        (*incount)++;
      }
    }
  }
}                               /* fillpicklists */


void
addmigration_covar (int ci, int li, int *oldmigcount, double *oldtlength,
                    double *form, double *form_inv, double *backm,
                    double *backm_inv, double migmultiplier)
{
  int newsis, edge;
  double mparamf, mparamb;
  struct edge *gtree = C[ci]->G[li].gtree;
  double mtime;
  int mcount;
  assert (C[ci]->G[li].mignum >= 0 && C[ci]->G[li].tlength > 0);

  /* determine current migration rate to use for update */
  mparamf = calcmrate (C[ci]->G[li].mignum, C[ci]->G[li].tlength);

  /* store information on edge, and sister edge if needed */
  //memset (&newedgemig, 0, sizeof (struct edgemiginfo));
  IMA_reset_edgemiginfo (&newedgemig);
  newedgemig.edgeid = edge = oldedgemig.edgeid;
  newedgemig.li = li;
  oldedgemig.li = li;
  if (edge < L[li].numgenes)
    newedgemig.upt = 0;
  else
    newedgemig.upt = gtree[gtree[edge].up[0]].time;
  newedgemig.fpop = gtree[gtree[edge].down].pop;
  newedgemig.pop = newedgemig.temppop = gtree[edge].pop;
  newedgemig.dnt = gtree[edge].time;
  newedgemig.mig[0].mt = -1;
  fillmiginfoperiods (ci, &newedgemig);
  if (gtree[edge].down == C[ci]->G[li].root)    /* simulate migration on the sister branch as well */
  {
    //memset (&newsismig, 0, sizeof (struct edgemiginfo));
    IMA_reset_edgemiginfo (&newsismig);
    newedgemig.fpop = -1;       //pop unknown, as edge must be determined by migration 
    if (gtree[gtree[edge].down].up[0] == edge)
      newsis = gtree[gtree[edge].down].up[1];
    else
      newsis = gtree[gtree[edge].down].up[0];

    if (newsis < L[li].numgenes)
      newsismig.upt = 0;
    else
      newsismig.upt = gtree[gtree[newsis].up[0]].time;

    newsismig.edgeid = newsis;
    newsismig.li = li;
    oldsismig.li = li;
    newsismig.fpop = -1;        //pop unknown, as edge must be determined by migration 
    newsismig.pop = newsismig.temppop = gtree[newsis].pop;
    newsismig.dnt = gtree[newsis].time;
    newsismig.mig[0].mt = -1;
    fillmiginfoperiods (ci, &newsismig);
  }
  else
  {
    newedgemig.fpop = gtree[gtree[edge].down].pop;
    newsismig.edgeid = -1;
    newsismig.mtall = 0;        // no second edge to deal with
  }
  getm (ci, mparamf * migmultiplier, &newedgemig, &newsismig);
  if (gtree[newedgemig.edgeid].down == C[ci]->G[li].root)
    gtree[C[ci]->G[li].root].pop = newedgemig.fpop;
  if (newedgemig.mtall == 0 && newsismig.mtall == 0)    // both new edges occur after the last split time
  {
    *form = 0;
    *form_inv = 0;
  }
  else
  {
    *form = getmprob (ci, mparamf * migmultiplier, &newedgemig, &newsismig);
    *form_inv =
      getmprob (ci, mparamf * (1.0 / migmultiplier), &newedgemig, &newsismig);
  }
  /* calculate probability of reverse update    */
  /* 7/11/08  - fixed a nasty bug here, 
     was causing wrong values for mtime and mcount, which in turn was causing wrong values
     for reverse update in calculation of hastings term.
     not sure just how it might have shaped the what it was doing to results */
  mtime = *oldtlength - oldedgemig.mtall - oldsismig.mtall + newedgemig.mtall;
  mcount = *oldmigcount - oldedgemig.mpall - oldsismig.mpall + newedgemig.mpall;
  if (newsismig.edgeid >= 0)
  {
    mtime += newsismig.mtall;
    mcount += newsismig.mpall;
  }
/* find the migation rate for the backward update */
  if (oldedgemig.mtall == 0 && oldsismig.mtall == 0)    // both new edges occur after the last split time
  {
    *backm = 0;
    *backm_inv = 0;
  }
  else
  {
    mparamb = calcmrate (mcount, mtime);        // is it correct to use the inverse multiplier for the reverse update ????
    *backm = getmprob (ci, mparamb * migmultiplier, &oldedgemig, &oldsismig);
    *backm_inv =
      getmprob (ci, mparamb * (1.0 / migmultiplier), &oldedgemig, &oldsismig);
  }
  *oldmigcount = mcount;
  *oldtlength = mtime;
}                               /* addmigration_covar */


/**********GLOBAL FUNCTIONS***/
void
init_updategenealogy_covar (void)
{
  int li;
  init_genealogy_weights (&holdallgweight);
  init_probcalc (&holdallpcalc);
  for (largestsamp = 0, li = 0; li < nloci; li++)
    if (largestsamp < L[li].numgenes)
      largestsamp = L[li].numgenes;
  init_holdgtree (&holdgtree, largestsamp);
  inpoplist = malloc (2 * largestsamp * sizeof (int));
  slidedist = malloc (MAX_NUM_GTREE_COVAR_UPDATES * sizeof (double));
  num_gtree_covar_updates = malloc (nloci * sizeof (int));
  for (li = 0; li < nloci; li++)
  {
    num_gtree_covar_updates[li] = (int) IMAX (MIN_NUM_GTREE_COVAR_UPDATES,
                                              (int)
                                              IMIN
                                              (MAX_NUM_GTREE_COVAR_UPDATES,
                                               (double) L[li].numlines /
                                               (NUM_GTREE_COVAR_UPDATES_DIVISOR
                                                * (double) numtreepops)));
  }
}                               // init_updategenealogy


void
free_updategenealogy_covar (void)
{
  free_genealogy_weights (&holdallgweight);
  free_probcalc (&holdallpcalc);
  free_holdgtree (&holdgtree, largestsamp);
  XFREE (inpoplist);
  XFREE (slidedist);
  XFREE (num_gtree_covar_updates);
}                               //free_updategenealogy

int
updategenealogy_covar (int ci, int li)
{
  int i, j;
  int ai, mpart;
  int edge, oldsis, newsis, freededge, accp;
  double newpdg, newpdg_a[MAXLINKED];
  double newplg;
  double pickededgeweight, migweight, metropolishastingsterm, U;
  double tpw;
  double Aterm[MAXLINKED], Atermsum;
  double tlengthpart;
  double slideparam, tempslidedist;
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  int rejectIS;
  double like;
  int incount[MAX_NUM_GTREE_COVAR_UPDATES];
  int cpop;
  double temp;
  int tmrcachange, topolchange;
  int netrootmove;
  double migmultiplier, slidemultiplier;
  static double prob_mig_forward_w_multiplier[MAX_NUM_GTREE_COVAR_UPDATES],
    prob_mig_forward_w_multiplier_inv[MAX_NUM_GTREE_COVAR_UPDATES],
    prob_mig_reverse_w_multiplier[MAX_NUM_GTREE_COVAR_UPDATES],
    prob_mig_reverse_w_multiplier_inv[MAX_NUM_GTREE_COVAR_UPDATES];
  double prob_mig_forward_w_multiplier_product = 0,
    prob_mig_forward_w_multiplier_inv_product = 0,
    prob_mig_reverse_w_multiplier_product = 0,
    prob_mig_reverse_w_multiplier_inv_product = 0;
  int startedge;
  int numnotrootpop0, numnotrootpop1;
  double tempslidemult, tempslideinv;
  double traplowexp;
  double holdt[MAXPERIODS];

  // initialize and make copies structures that hold quantities for calculating prob of genealogy
  copy_gtree (li, G, 1);
  copy_treeinfo (&holdgtree.gweight, &G->gweight);

  copy_treeinfo (&holdallgweight, &C[ci]->allgweight);
  copy_probcalc (&holdallpcalc, &C[ci]->allpcalc);
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];

  do
  {
    startedge = randposint (L[li].numlines);
  } while (gtree[startedge].down == -1);  // include rootpop || gtree[startedge].pop == C[ci]->rootpop);  //(gtree[startedge].down == -1);

  cpop = gtree[startedge].pop;
  slideparam = get_popbranchlengths (cpop, &numnotrootpop0, gtree, G, ci, li);
  migmultiplier = (bitran ())? GTREE_COVAR_MIGRATION_MULTIPLIER : 1.0 / GTREE_COVAR_MIGRATION_MULTIPLIER;
  slidemultiplier = (bitran ())? GTREE_COVAR_SLIDEPARAM_MULTIPLIER : 1.0 / GTREE_COVAR_SLIDEPARAM_MULTIPLIER;
  // Atermsum only used for Stepwise mutation model
  Atermsum = 0;
  migweight = 0;
  tlengthpart = G->tlength;
  mpart = G->mignum;
  netrootmove = 0;
  for (i = 0; i < num_gtree_covar_updates[li]; i++)     // main gtree updating loop 
  {
    fillpicklists (li, gtree, G, incount + i, cpop, inpoplist);
    if (i > 0)
    {
      assert (incount[i] > 0);
      edge = inpoplist[randposint (incount[i])];
      assert (gtree[edge].pop == cpop);
    }
    else
    {
      edge = startedge;
    }
    freededge = gtree[edge].down;
    if ((oldsis = gtree[freededge].up[0]) == edge)
      oldsis = gtree[freededge].up[1];

    /* copy information on the edge,  and if it connects to the root, then the sister edge as well */
    if (freededge == G->root)
      fillmiginfo (ci, li, gtree, edge, oldsis);
    else
      fillmiginfo (ci, li, gtree, edge, -1);
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      storeAinfo (ci, li, gtree, edge, oldsis, freededge);

    // remove any migrations  from the slidingedge 
    gtree[edge].mig[0].mt = -1;

    // pick the slidedistance
    slidedist[i] = normdev (0.0, slideparam * slidemultiplier);
    // join the sister and the down branches at the point where edge used to connect, this frees up the down branch 
    joinsisdown (ci, li, oldsis, &tmrcachange);

    // do the slide and identify the new sister branch and where the new connection point for the edge is 
    newsis = oldsis;
    tempslidedist = slidedist[i];
    if (mprior == 0)
      slider_nomigration (ci, li, edge, &newsis, &(gtree[edge].time),
                          &tempslidedist);
    else
      slider (ci, li, edge, &newsis, &(gtree[edge].time), &tempslidedist);

// now separate the new sister branch into a shorter sis branch and a down branch 
    splitsisdown (ci, li, edge, freededge, newsis);
    netrootmove |= rootmove;

    // add migration events  (must do this even if cpop ==C[ci]->rootpop because the sister branch might have had migration
    if (!modeloptions[NOMIGRATION])
    {
      addmigration_covar (ci, li, &mpart, &tlengthpart,
                          prob_mig_forward_w_multiplier + i,
                          prob_mig_forward_w_multiplier_inv + i,
                          prob_mig_reverse_w_multiplier + i,
                          prob_mig_reverse_w_multiplier_inv + i,
                          migmultiplier);
      prob_mig_reverse_w_multiplier_product +=
        prob_mig_reverse_w_multiplier[i];
      prob_mig_reverse_w_multiplier_inv_product +=
        prob_mig_reverse_w_multiplier_inv[i];
      prob_mig_forward_w_multiplier_product +=
        prob_mig_forward_w_multiplier[i];
      prob_mig_forward_w_multiplier_inv_product +=
        prob_mig_forward_w_multiplier_inv[i];
    }

    // copy the migration info in newedgemig and newsismig  to the genealogy
    copynewmig_to_gtree (ci, li);

    // update ancestral allele values for models with STRs - need to put A and dlikeA values into copyedge so that finishSWupdateA() can work 

    switch (L[li].model)
    {
    case STEPWISE:
      {
        if (i == 0)
          for (ai = 0; ai < L[li].nlinked; ai++)
            newpdg_a[ai] = 0;
        for (ai = 0; ai < L[li].nlinked; ai++)
        {
          newpdg_a[ai] +=
            finishSWupdateA (ci, li, ai, edge, freededge, oldsis, newsis,
                             G->uvals[ai], &Aterm[ai]);
          //checklikelihoodSW(ci, li,ai, G->u[ai].mcinf.val);  
          Atermsum += Aterm[ai];
        }
        break;
      }
    case JOINT_IS_SW:
      if (i == 0)
        for (ai = 1; ai < L[li].nlinked; ai++)
          newpdg_a[ai] = 0;
      for (ai = 1; ai < L[li].nlinked; ai++)
      {
        newpdg_a[ai] +=
          finishSWupdateA (ci, li, ai, edge, freededge, oldsis, newsis,
                           G->uvals[ai], &Aterm[ai]);
        //checklikelihoodSW(ci, li,ai, G->u[ai].mcinf.val);  
        Atermsum += Aterm[ai];
      }
      break;
    }
  }                             // updating loop

  //calculate Hastings for forward changes 

  // loop forward through the updates,  for each update there was an incount and a notincount
  for (j = 0, tempslidemult = tempslideinv = 1 /*, tempprod = 0 */ ;
       j < num_gtree_covar_updates[li]; j++)
  {
    if (j > 0)
      temp = 1 / (double) incount[j];
    else
      temp = 1;                 // include rootpop     1/(double) numnotrootpop0;//temp = 1;
    tempslidemult *=
      normprob (0.0, slideparam * slidemultiplier, slidedist[j]) * temp;
    tempslideinv *=
      normprob (0.0, slideparam * (1.0 / slidemultiplier), slidedist[j]) * temp;
  }
  pickededgeweight = -log (tempslidemult + tempslideinv);
  // trap extremely small values that cause a zero to go to log  in migweight
  traplowexp =
    DMAX (DMAX
          (DMAX
           (prob_mig_reverse_w_multiplier_product,
            prob_mig_reverse_w_multiplier_inv_product),
           prob_mig_forward_w_multiplier_product),
          prob_mig_forward_w_multiplier_inv_product);
  if (traplowexp < -700)
    migweight = traplowexp;
  else
    migweight =
      log ((exp (prob_mig_reverse_w_multiplier_product) +
            exp (prob_mig_reverse_w_multiplier_inv_product)) /
           (exp (prob_mig_forward_w_multiplier_product) +
            exp (prob_mig_forward_w_multiplier_inv_product)));

  //calculate Hastings for reverse changes, loop backwards through the updates,  for each update there was an incount
  slideparam = get_popbranchlengths (cpop, &numnotrootpop1, gtree, G, ci, li);
  for (j = num_gtree_covar_updates[li], tempslidemult = tempslideinv =
       1 /*, tempprod = 0 */ ; j > 0; j--)
  {
    if (j < num_gtree_covar_updates[li])
      temp = 1 / (double) incount[j - 1];
    else
      temp = 1;                 // include rootpop  1/(double) numnotrootpop1; //temp  = 1;
    //tempprod += log(normprob (0,slideparam, -slidedist[j-1]) * temp );
    tempslidemult *= normprob (0.0, slideparam * slidemultiplier, -slidedist[j - 1]) * temp;
    tempslideinv *= normprob (0.0, slideparam * (1.0 / slidemultiplier), -slidedist[j - 1]) * temp;
  }
  //pickededgeweight += tempprod;
  pickededgeweight += log (tempslidemult + tempslideinv);

// determine all the weights needed for calculating the probability of the genealogy
  setzero_genealogy_weights (&G->gweight);
  treeweight (ci, li);
  sum_subtract_treeinfo (&C[ci]->allgweight, &G->gweight,
                         &(holdgtree.gweight));

/* calculate P(D|G)  for new genealogy */
  rejectIS = 0;                 /* in case P(D|G) for IS model is zero */
  switch (L[li].model)
  {
  case HKY:
    if (assignmentoptions[JCMODEL] == 1)
    {
      newpdg_a[0] = likelihoodJC (ci, li, G->uvals[0]);
      newpdg = newpdg_a[0];
    }
    else
    {
      newpdg = newpdg_a[0] =
        likelihoodHKY (ci, li, G->uvals[0], G->kappaval, -1, -1, -1, -1);
    }
    break;
  case INFINITESITES:
    newpdg = newpdg_a[0] = like = likelihoodIS (ci, li, G->uvals[0]);
    rejectIS = (like == REJECTINFINITESITESCONSTANT);
    break;
  case STEPWISE:
    {
      for (ai = 0, newpdg = 0; ai < L[li].nlinked; ai++)
      {
        newpdg_a[ai] += G->pdg_a[ai];
        newpdg += newpdg_a[ai];
      }
      //            checklikelihoodSW(ci, li,G->uvals[ai]);  
      break;
    }
  case JOINT_IS_SW:
    newpdg = newpdg_a[0] = likelihoodIS (ci, li, G->uvals[0]);
    rejectIS = (newpdg == REJECTINFINITESITESCONSTANT);
    for (ai = 1; ai < L[li].nlinked; ai++)
    {
      newpdg_a[ai] += G->pdg_a[ai];
      newpdg += newpdg_a[ai];
    }
    break;

  }
  if (calcoptions[CALCLEVINELIKELIHOOD])
  {
    if (topolchange)
    {
      newplg = likelihoodDG (ci, li);
    }
    else
    {
      newplg = G->plg;
    }
  }
  else
  {
    newplg = 0;
  }

  accp = 0;

/* final weight calculation */
/* tpw is the ratio of new and old prior probability of the genealogies.  It is actually the ratio of the total across all loci,  but
since only genealogy li is being changed at the present time,  the ratio works out to just be the ratio for genealogy li */
  copy_probcalc (&holdallpcalc, &C[ci]->allpcalc);
  if (rejectIS == 0)
  {
    // the metropolis term includes p(D|G) and p(G),  
    tpw = -C[ci]->allpcalc.probg;
    //integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
    integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight,
                         &C[ci]->allpcalc, &holdallpcalc, &holdt[0]);

    tpw += C[ci]->allpcalc.probg;
    //metropolishastingsterm = tpw + newpdg - G->pdg;
    metropolishastingsterm = tpw + gbeta * (newpdg - G->pdg + newplg - G->plg);
    U = uniform ();
    //  assert (beta[ci] * metropolishastingsterm + pickededgeweight + migweight + Atermsum > -1e200 && beta[ci] * metropolishastingsterm + pickededgeweight + migweight + Atermsum < 1e200);
    metropolishastingsterm = exp (beta[ci] * metropolishastingsterm + pickededgeweight + migweight + Atermsum);
//      assert (metropolishastingsterm >= 0);
    if (metropolishastingsterm > 1.0 || metropolishastingsterm > U)
    {
      /* accept the update */
      C[ci]->allpcalc.pdg += newpdg - G->pdg;
      C[ci]->allpcalc.plg += newplg - G->plg;   /* SANGCHUL: plus Levine's likelihood */
      G->pdg = newpdg;
      G->plg = newplg;
      for (ai = 0; ai < L[li].nlinked; ai++)
        G->pdg_a[ai] = newpdg_a[ai];
      if (L[li].model == HKY)
      {
        copyfraclike (ci, li);
        storescalefactors (ci, li);
      }
      accp = 1;
    }
  }

  /* reject the update */
  if (accp == 0)
  {
    // copy back all the weights and results associated with calculating the probability of the genealogy 
    copy_gtree (li, G, 0);
    copy_treeinfo (&G->gweight, &holdgtree.gweight);
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc);
    if (L[li].model == HKY)     // reset HKY terms
      restorescalefactors (ci, li);
    tmrcachange = 0;
  }
  return accp;
}                               /* updategenealogy_covar */

#undef MAX_NUM_GTREE_COVAR_UPDATES
#undef MIN_NUM_GTREE_COVAR_UPDATES
#undef NUM_GTREE_COVAR_UPDATES_DIVISOR
#undef GTREE_COVAR_MIGRATION_MULTIPLIER
#undef GTREE_COVAR_SLIDEPARAM_MULTIPLIER
#undef SLIDESTDVMAX

