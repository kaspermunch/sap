/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS

#include "imamp.h"
#include "update_gtree_common.h"
#include "updateassignment.h"


/*********** LOCAL STUFF **********/

static struct genealogy_weights holdallgweight_t_RY;
static struct genealogy_weights holdgweight_t_RY[MAXLOCI];
static struct probcalc holdallpcalc_t_RY;
static int largestsamp;
static int **skipflag;

static double beforesplit (int tnode, double t, double newt, double tu,
                           double tl, double ptime);
static double aftersplit (int tnode, double t, double newt, double tu,
                          double tl, double ptime);
static void storetreestats_all_loci (int ci, int mode);
static int edgerecurse (int ci, int li, struct genealogy *G,
                        struct edge *gtree, int edge, double tau_d,
                        double tau_u, int period, int *flag);
static void setzerorecurse (int ci, int li, struct genealogy *G,
                            struct edge *gtree, int i, int *flag);
static void setedgeflag (int ci, int li, struct genealogy *G,
                         struct edge *gtree, double tau_d, double tau_u,
                         int period, int *flag);

static double forward (double ptime, double oldt, double r, double diff,
                       double t_u_prior);
static double backward (double ptime, double newt, double r, double diff,
                        double t_u_prior);

/********* LOCAL FUNCTIONS **************/

void
storetreestats_all_loci (int ci, int mode)
{
  static double holdlength[MAXLOCI], holdtlength[MAXLOCI];
  static double holdroottime[MAXLOCI];
  static int holdroot[MAXLOCI];
  static int holdmig[MAXLOCI];
  int li;
  if (mode == 0)
  {
    for (li = 0; li < nloci; li++)
    {
      holdlength[li] = C[ci]->G[li].length;
      holdtlength[li] = C[ci]->G[li].tlength;
      holdroottime[li] = C[ci]->G[li].roottime;
      holdroot[li] = C[ci]->G[li].root;
      holdmig[li] = C[ci]->G[li].mignum;
    }
  }
  else
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].length = holdlength[li];
      C[ci]->G[li].tlength = holdtlength[li];
      C[ci]->G[li].mignum = holdmig[li];
      C[ci]->G[li].roottime = holdroottime[li];
      C[ci]->G[li].root = holdroot[li];
    }
  }
  return;
}                               // storetreestats  

double
aftersplit (int tnode, double oldt, double newt, double tau_d, double tau_u,
            double ptime)
{
  /* 2_24_09  nasty bug here.   For some reason had turned off the limiting case  that assumes tau_d is infinite
    when tnode == (lastperiodnumber - 1).   However when time is very close to oldt,  this would break down and we 
    return a new time that is less than newt when the old time was greater than oldt.
    So turn back on the usage of the limit as tau_d goes to infinity */ 
  if (tnode == lastperiodnumber - 1)
     return ptime + newt - oldt;
  else
    return tau_d - (tau_d - newt) * (tau_d - ptime) / (tau_d - oldt);
}

double
beforesplit (int tnode, double oldt, double newt, double tau_d, double tau_u,
             double ptime)
{
  if (tnode == 0)
    return ptime * newt / oldt;
  else
    return tau_u + (ptime - tau_u) * (newt - tau_u) / (oldt - tau_u);
}

/* edge recurse is a bit tricky because the possible population states of all points on the edge must be checked */
int
edgerecurse (int ci, int li, struct genealogy *G, struct edge *gtree,
             int edge, double tau_d, double tau_u, int period, int *flag)
{
  int j, pop, popcrosst;
  double uptime, t, td;

  // check the top of the edge
  uptime = (edge < L[li].numgenes) ? 0.0 : gtree[gtree[edge].up[0]].time;
  t = DMAX (tau_u, uptime);
  pop = nowedgepop (ci, &gtree[edge], t);
  if (C[ci]->poptree[pop].e <= (lastperiodnumber - 1)
      && t == C[ci]->tvals[C[ci]->poptree[pop].e - 1])
    pop = C[ci]->poptree[pop].down;
  if (pop == C[ci]->droppops[period + 1][0] ||
      pop == C[ci]->droppops[period + 1][1] ||
      pop == C[ci]->addpop[period + 1])
  {
    return 0;                   // get out - the edge is relevant 
  }
  j = 0;
  while (gtree[edge].mig[j].mt >= 0.0 && gtree[edge].mig[j].mt < t)
    j++;
  td = DMIN (gtree[edge].time, tau_d);
  popcrosst = 0;
  while ((C[ci]->poptree[pop].e <= (lastperiodnumber - 1)
          && t > C[ci]->tvals[C[ci]->poptree[pop].e - 1]) || t < td)
  {
    if (gtree[edge].mig[j].mt >= 0.0)
    {
      t = gtree[edge].mig[j].mt;
      if (popcrosst == 1)
      {
        pop = gtree[edge].mig[j].mp;
        j++;
      }
      popcrosst = 1 - popcrosst;
    }
    else
    {
      t = td;
    }
    if (pop == C[ci]->droppops[period + 1][0] ||
        pop == C[ci]->droppops[period + 1][1] ||
        pop == C[ci]->addpop[period + 1])
    {
      return 0;                 // get out - the edge is relevant 
    }
    if (C[ci]->poptree[pop].e <= (lastperiodnumber - 1)
        && t > C[ci]->tvals[C[ci]->poptree[pop].e - 1])
      pop = C[ci]->poptree[pop].down;
  }
  if (edge >= L[li].numgenes)
  {
    flag[gtree[edge].up[0]] =
      edgerecurse (ci, li, G, gtree, gtree[edge].up[0], tau_d, tau_u, period,
                   flag);
    flag[gtree[edge].up[1]] =
      edgerecurse (ci, li, G, gtree, gtree[edge].up[1], tau_d, tau_u, period,
                   flag);
    if ((flag[gtree[edge].up[0]] == -1 || flag[gtree[edge].up[0]] == 1)
        && (flag[gtree[edge].up[1]] == -1 || flag[gtree[edge].up[1]] == 1))
    {
      return 1;
    }
    else
    {
      flag[gtree[edge].up[0]] = 0;
      flag[gtree[edge].up[1]] = 0;
      return 0;
    }
  }
  else
  {
    return 1;
  }
  assert (0);                   // should not get here
}

/* reset skipflag values to 0 
necessary because after edgerecurse() it can happen that an edge that crosses the tau_d 
line will have a skipflag value of 0,  but one of its descendants will have a 
skipflag value of 1
*/
void
setzerorecurse (int ci, int li, struct genealogy *G, struct edge *gtree,
                int i, int *flag)
{
  if (flag[i] != -1)
  {
    flag[i] = 0;
    if (i >= L[li].numgenes)
    {
      setzerorecurse (ci, li, G, gtree, gtree[i].up[0], flag);
      setzerorecurse (ci, li, G, gtree, gtree[i].up[1], flag);
    }
  }
}                               /* setzerorecurse */

/* determines the value of skipflag */
void
setedgeflag (int ci, int li, struct genealogy *G, struct edge *gtree,
             double tau_d, double tau_u, int period, int *flag)
{
  int i;
  double uptime;

  for (i = 0; i < L[li].numlines; i++)
  {
    if (period == npops - 1)
    {
      flag[i] = 0;
    }
    else
    {
      uptime = (i < L[li].numgenes) ? 0 : gtree[gtree[i].up[0]].time;
      if (uptime > tau_d)
        flag[i] = -2;
      else if (gtree[i].time < tau_u)
        flag[i] = -1;
      else
        flag[i] = 0;
    }
  }
  if (period < npops - 1)
  {
    for (i = 0; i < L[li].numlines; i++)
    {
      uptime = (i < L[li].numgenes) ? 0.0 : gtree[gtree[i].up[0]].time;
      if (uptime < tau_d && gtree[i].time > tau_d)
      {
        assert (flag[i] == 0);
        flag[i] =
          edgerecurse (ci, li, G, gtree, i, tau_d, tau_u, period, flag);
        if (flag[i] == 0 && i >= L[li].numgenes)
          setzerorecurse (ci, li, G, gtree, i, flag);
      }
    }
  }
}                               /* setdgeflag */


double
forward (double ptime, double oldt, double r, double diff, double t_u_prior)
{
  if (ptime > t_u_prior)
  {
    if (ptime < oldt)
      return t_u_prior + (ptime - t_u_prior) * r;
    else
      return ptime + diff;
  }
  else
  {
    return ptime;
  }
}

double
backward (double ptime, double newt, double r, double diff, double t_u_prior)
{
  if (ptime > t_u_prior)
  {
    if (ptime < newt)
      return t_u_prior + (ptime - t_u_prior) / r;
    else
      return ptime - diff;
  }
  else
  {
    return ptime;
  }
}


/*************GLOBAL FUNCTIONS ******************/


void
init_t_RY (void)
{
  int li, j;
  init_genealogy_weights (&holdallgweight_t_RY);
  for (li = 0; li < nloci; li++)
    init_genealogy_weights (&holdgweight_t_RY[li]);
  init_probcalc (&holdallpcalc_t_RY);
  for (largestsamp = 0, j = 0; j < nloci; j++)
    if (largestsamp < L[j].numlines)
      largestsamp = L[j].numlines;
  skipflag = alloc2Dint (nloci, 2 * largestsamp - 1);
  RY2_t_upinf = calloc ((size_t) npops - 1, sizeof (struct update_rate_calc));
}                               // init_changet_RY


void
free_t_RY (void)
{
  int li;
  free_genealogy_weights (&holdallgweight_t_RY);
  for (li = 0; li < nloci; li++)
  {
    free_genealogy_weights (&holdgweight_t_RY[li]);
  }
  free_probcalc (&holdallpcalc_t_RY);
  free2D ((void **) skipflag, nloci);
  XFREE (RY2_t_upinf);
  RY2_t_upinf = NULL;
}                               // free_changet_RY

/* changet_RY2  
does a Rannala Yang style update,  but only changes times int the targeted period. Those times in higher numbered
periods are just shifted up or down by a constant */
/* a tricky issue arises in this update:
the update for time i causes changes to time i and all times j> i 
this means that it is possible to propose an update that would mean some times to be beyond the prior for t. 
This is unlikely under the model with a splitting rate parameter, because time priors are set to be quite large under this model.
But it is a big issue for models without splitting rates that have smaller values for the splitting time priors. 

So how to handle this in models with low splitting time priors ?  
the key questions are 
	- is what should the upper bound be on the range for which the new time is picked
	- how to handle times that make higher times over the limit

The answer to the first question is that it seems not to matter too much. 
For example we can use either
	t_d_prior = C[ci]->tvals[timeperiod+1];
	t_d_prior = T[lastperiodnumber-1].pr.max;
and still get the same results

For the answer to the second,  I'm just rejecting outright any update attempt that picks a time that puts 
the last time over the line:
if ( newt > oldt  && C[ci]->tvals[lastperiodnumber -1] + (newt - oldt) > T[lastperiodnumber - 1].pr.max)
	return 0
this seems to work (i.e. returns the correct prior on splitting times when no data is used)
*/

int
changet_RY2 (int ci, int timeperiod)
{
  double metropolishastingsterm, newt, oldt;
  double pdgnew[MAXLOCI + MAXLINKED], pdgnewsum, pdgoldsum, probgnewsum,
    temppdg;
  double tpw;
  int li, i, j, ai, ui;
  int ec, em;
  double U;
  double r;
  struct genealogy *G;
  struct edge *gtree;
  double tdiff, t_u, t_d_prior, t_u_prior;
  double holdt[MAXPERIODS];

  assert (timeperiod < lastperiodnumber - 1);
  /* select a new split interval time */
  t_u = (timeperiod == 0) ? 0 : C[ci]->tvals[timeperiod - 1];
  t_u_prior = DMAX (T[timeperiod].pr.min, t_u);
  t_d_prior = C[ci]->tvals[timeperiod + 1];
//      t_d_prior = T[lastperiodnumber-1].pr.max;

  oldt = C[ci]->tvals[timeperiod];
  //twin = (t_d_prior - t_u_prior) / (3 * log(nloci+1)* (npops - timeperiod)*(npops - timeperiod));  
  newt = getnewt (timeperiod, t_u_prior, t_d_prior, oldt, 2);
  assert (newt < T[timeperiod].pr.max);
  if (newt > oldt
      && C[ci]->tvals[lastperiodnumber - 1] + (newt - oldt) >
      T[lastperiodnumber - 1].pr.max)
  {
    return 0;
  }
  r = (newt - t_u_prior) / (oldt - t_u_prior);
  tdiff = newt - oldt;

  pdgnewsum = 0;
  probgnewsum = 0;
  ec = em = 0;
  storetreestats_all_loci (ci, 0);
  copy_treeinfo (&holdallgweight_t_RY, &C[ci]->allgweight);
  copy_probcalc (&holdallpcalc_t_RY, &C[ci]->allpcalc);
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];

  pdgoldsum = C[ci]->allpcalc.pdg;
  setzero_genealogy_weights (&C[ci]->allgweight);
  C[ci]->tvals[timeperiod] = newt;
  C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time = newt;
  C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;


  for (i = timeperiod + 1; i < lastperiodnumber; i++)
  {
    C[ci]->tvals[i] = forward (C[ci]->tvals[i], oldt, r, tdiff, t_u_prior);
    C[ci]->poptree[C[ci]->droppops[i + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[i + 1][1]].time = C[ci]->tvals[i];
  }
  for (i = 0; i < nurates; i++)
    pdgnew[i] = 0;
  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;

    copy_treeinfo (&holdgweight_t_RY[li], &G->gweight);
    for (i = 0; i < L[li].numlines; i++)
      if (gtree[i].down != -1)
      {
        ec += (gtree[i].time > t_u_prior && gtree[i].time < oldt);
        gtree[i].time = forward (gtree[i].time, oldt, r, tdiff, t_u_prior);
        j = 0;
        while (gtree[i].mig[j].mt > -0.5)
        {
          em += (gtree[i].mig[j].mt > t_u_prior && gtree[i].mig[j].mt < oldt);
          gtree[i].mig[j].mt =
            forward (gtree[i].mig[j].mt, oldt, r, tdiff, t_u_prior);
          j++;
        }
      }
    G->roottime = forward (G->roottime, oldt, r, tdiff, t_u_prior);
    setzero_genealogy_weights (&G->gweight);
    treeweight (ci, li);


    sum_treeinfo (&C[ci]->allgweight, &G->gweight);
    ai = 0;
    ui = L[li].uii[ai];
    switch (L[li].model)
    {
      assert (pdgnew[ui] == 0);
    case HKY:
      if (assignmentoptions[JCMODEL] == 1)
      {
        temppdg = pdgnew[ui] =
          likelihoodJC (ci, li, G->uvals[0]);
      }
      else
      {
        temppdg = pdgnew[ui] =
          likelihoodHKY (ci, li, G->uvals[0], G->kappaval, -1, -1, -1, -1);
      }
      break;
    case INFINITESITES:
      if (calcoptions[DONTCALCLIKELIHOODMUTATION])
        temppdg = pdgnew[ui] = 0;
      else
        temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      break;
    case STEPWISE:
      temppdg = 0;
      for (; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    case JOINT_IS_SW:
      if (calcoptions[DONTCALCLIKELIHOODMUTATION])
        temppdg = pdgnew[ui] = 0;
      else
        temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      for (ai = 1; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    }
    pdgnewsum += temppdg;
  }
  tpw = gbeta * (pdgnewsum - pdgoldsum);

  //integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
  integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight_t_RY,
                       &C[ci]->allpcalc, &holdallpcalc_t_RY, &holdt[0]);
  tpw += C[ci]->allpcalc.probg - holdallpcalc_t_RY.probg;
  assert (ODD (ec) == 0);
  ec /= 2;

  metropolishastingsterm = exp (beta[ci] * tpw + (ec + em) * log (r));
  //assert(metropolishastingsterm >= 0 && metropolishastingsterm < 1e200);
  U = uniform ();
  if (metropolishastingsterm >= 1.0 || metropolishastingsterm > U)
  {
    for (li = 0; li < nloci; li++)

    {
      C[ci]->G[li].pdg = 0;
      for (ai = 0; ai < L[li].nlinked; ai++)
      {
        C[ci]->G[li].pdg_a[ai] = pdgnew[L[li].uii[ai]];
        C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
      }
      //  assert (C[ci]->G[li].pdg < 0);
      if (L[li].model == HKY)
      {
        storescalefactors (ci, li);
        copyfraclike (ci, li);
      }
    }
    C[ci]->allpcalc.pdg = pdgnewsum;
    return 1;
  }
  else
  {
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_RY);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_RY);
    assert (pdgoldsum == C[ci]->allpcalc.pdg);
    C[ci]->tvals[timeperiod] = oldt;
    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = oldt;
    for (i = timeperiod + 1; i < lastperiodnumber; i++)
    {
      C[ci]->tvals[i] = backward (C[ci]->tvals[i], newt, r, tdiff, t_u_prior);
      C[ci]->poptree[C[ci]->droppops[i + 1][0]].time =
        C[ci]->poptree[C[ci]->droppops[i + 1][1]].time = C[ci]->tvals[i];
    }
    for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      storetreestats_all_loci (ci, 1);
      copy_treeinfo (&G->gweight, &holdgweight_t_RY[li]);
      for (i = 0; i < L[li].numlines; i++)
        if (gtree[i].down != -1)
        {
          gtree[i].time = backward (gtree[i].time, newt, r, tdiff, t_u_prior);
          j = 0;
          while (gtree[i].mig[j].mt > -0.5)
          {
            gtree[i].mig[j].mt =
              backward (gtree[i].mig[j].mt, newt, r, tdiff, t_u_prior);
            j++;
          }
        }

    }
    for (li = 0; li < nloci; li++)
    {
      if (L[li].model == HKY)
        restorescalefactors (ci, li);
      /* have to reset the dlikeA values in the trees for stepwise model */
      if (L[li].model == STEPWISE)
        for (ai = 0; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      if (L[li].model == JOINT_IS_SW)
        for (ai = 1; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
    }
    return 0;
  }
}                               /* changet_RY2 */


/* 
Notes on changet_RY()  implements updating of Rannala and Yang (2003)   

This application is pretty much the same as theirs - changing times on a species tree that contains a gene tree. 
The big difference is that IM includes migration.   This means that we have to count migration events and change 
migration times in the same was as we change coalescent times and count how many get changed. 

R&Y also only change coalescent times in populations affected by a changed t.  But because we have migration
there is more entanglement between populations.  It seems best to change all times within an interval that is 
affected by a changing t. 


in R&Y usage
u (upper) means older, deeper in the tree
l (lower)  means younger, more recent in the tree
but this causes confusion with jhey usage in which 
u means upper - more recent. 

so use  u - for upper (more recent)
use d for down  (older)

tau current split time
tau*  new value 
tau_d - older node time (i.e. time of next oldest node - deeper in the tree)
tau_u - more recent node time (i.e. time of next youngest node - more recent in time)
tau_d  > tau > tau_u

if tau is the first node,  then tau_u = 0. 
if tau is the last node, then tau_d = infinity

for an event at time t where tau_u < t < tau_d

if t > tau  see aftersplit()  t is older than the current split time 
t* = tau_d - (tau_d - tau*)(tau_d - t)/(tau_d-tau)  (A7)

if t <= tau  see beforesplit()
t* = tau_u + (tau* - tau_u)(t - tau_u)/(tau - tau_u) (A8)

m is the number of events moved using A7,  n is the number moved using A8

then Hastings term for the update is:
 tau_d - tau*      tau* - tau_u 
(------------)^m  (------------)^n
 tau-u - tau        tau - tau_u

 For IM,  we use the same except m and n include both includes migation and coalescent events

For edges where downtime < tau_d || uptime > tau_d  are not involved in the update
For any edge that spends zero time in either splitpop  or the ancestral pop, during the tau_u/tau_d interval
it is  possible to not update the coalescent time or migration times of 

The difficulty is that it is possible that an uptime for one edge gets moved because the times on one of its daughter edges got moved. 
This means that for checking about skipping an edge, because it is not involved in any
population associated with the splittin time update
we need to check entire groups of branches that descend from 
an edge that crosses the tau_u boundary. 

use a scoring system for edges  - see setedgeflag()
-1 to ignore because out of bounds above
-2 to ignore because out of bounds below
0 to  deal with, edge is relevant
1  to ignore because not in splitpops or ancestral pop

set up a recursive function  see edgerecurse()
for an edge that crosses the tau_d line,  check to see if that edge
and all descendent edges that to not cross tau_l  are in the 
populations affected by the splitting time update
If all of those edges are not relevant then 
they all get a skipflag value of 1
If any one of them does spend any time in any of the 
populations involved int the population split
then all of the edges have their skipflag value 
set to 0  - see also setzerorecurse()

*/

/* let u refer to the more recent time  and d to the older time  */
int
changet_RY1 (int ci, int timeperiod)    // after Rannala and Yang (2003)  - rubberband method
{
  double metropolishastingsterm, newt, oldt;
  double pdgnew[MAXLOCI + MAXLINKED], pdgnewsum, pdgoldsum, probgnewsum,
    temppdg;
  double t_u_hterm, t_d_hterm, tpw;
  int li, i, j, ecd, ecu, emd, emu, ai, ui;
  double U;
  struct genealogy *G;
  struct edge *gtree;
  double t_d, t_u, t_u_prior, t_d_prior;
  double holdt[MAXPERIODS];


  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogy (ci);
  }

  t_u = (timeperiod == 0) ? 0 : C[ci]->tvals[timeperiod - 1];
  t_d =
    (timeperiod ==
     (lastperiodnumber - 1)) ? TIMEMAX : C[ci]->tvals[timeperiod + 1];
  t_d_prior = DMIN (T[timeperiod].pr.max, t_d);
  t_u_prior = DMAX (T[timeperiod].pr.min, t_u);
  if (modeloptions[SPLITTINGRATEPARAMETER])
  {
    t_u = t_u_prior;
    t_d = t_d_prior;
  }
  oldt = C[ci]->tvals[timeperiod];
  newt = getnewt (timeperiod, t_u_prior, t_d_prior, oldt, 1);
  /*
     if (timeperiod == lastperiodnumber - 1 && modeloptions[SPLITTINGRATEPARAMETER])
     twin = (t_d_prior - t_u_prior) / (log(nloci+2)*2* npops);  //  smaller window for last time period
     else
     twin = (t_d_prior - t_u_prior) / (log((nloci+2)/2.0)*(npops - timeperiod ));  //  smaller window for more recent times
   */

  t_u_hterm = (newt - t_u) / (oldt - t_u);
  if (timeperiod == lastperiodnumber - 1)
  {
    t_d_hterm = 1;
  }
  else
  {
    t_d_hterm = (t_d - newt) / (t_d - oldt);
  }

  copy_treeinfo (&holdallgweight_t_RY, &C[ci]->allgweight);
  copy_probcalc (&holdallpcalc_t_RY, &C[ci]->allpcalc);
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];


  pdgoldsum = C[ci]->allpcalc.pdg;
  setzero_genealogy_weights (&C[ci]->allgweight);
  ecd = ecu = emd = emu = 0;
  pdgnewsum = 0;
  probgnewsum = 0;
  storetreestats_all_loci (ci, 0);
  C[ci]->tvals[timeperiod] = newt;
  for (i = 0; i < nurates; i++)
    pdgnew[i] = 0;
  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;
    copy_treeinfo (&holdgweight_t_RY[li], &G->gweight);
    setedgeflag (ci, li, G, gtree, t_d, t_u, timeperiod, skipflag[li]);
    for (i = 0; i < L[li].numlines; i++)
    {
      if (gtree[i].down != -1
          && skipflag[li][i] == 0 /* use with setedgeflag */ )
      {
        if (gtree[i].time <= oldt && gtree[i].time > t_u)

        {
          assert (skipflag[li][i] == 0);
          gtree[i].time =
            beforesplit (timeperiod, oldt, newt, t_d, t_u, gtree[i].time);
          assert (gtree[i].time != newt);
          ecu++;
        }
        else
        {
          if (gtree[i].time > oldt && gtree[i].time < t_d)
          {
            assert (skipflag[li][i] == 0);
            gtree[i].time =
              aftersplit (timeperiod, oldt, newt, t_d, t_u, gtree[i].time);
            assert (gtree[i].time != newt);
            ecd++;
          }
          //else  do not change the time
        }
        j = 0;
        while (gtree[i].mig[j].mt > -0.5)
        {
          assert (gtree[i].mig[j].mt < C[0]->tvals[lastperiodnumber]);
          if (gtree[i].mig[j].mt <= oldt && gtree[i].mig[j].mt > t_u)
          {
            gtree[i].mig[j].mt =
              beforesplit (timeperiod, oldt, newt, t_d, t_u,
                           gtree[i].mig[j].mt);
            emu++;
          }
          else
          {
            assert (oldt < C[0]->tvals[lastperiodnumber]);
            if (gtree[i].mig[j].mt > oldt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplit (timeperiod, oldt, newt, t_d, t_u,
                            gtree[i].mig[j].mt);
              emd++;
            }
            // else no need to change the time
          }
          j++;
        }
      }
    }
    if (G->roottime <= oldt && G->roottime > t_u
        && skipflag[li][G->root] == 0)
      G->roottime =
        beforesplit (timeperiod, oldt, newt, t_d, t_u, G->roottime);
    else if (G->roottime > oldt && G->roottime < t_d
             && skipflag[li][G->root] == 0)
      G->roottime =
        aftersplit (timeperiod, oldt, newt, t_d, t_u, G->roottime);
    setzero_genealogy_weights (&G->gweight);
        
    treeweight (ci, li);

    sum_treeinfo (&C[ci]->allgweight, &G->gweight);
    ai = 0;
    ui = L[li].uii[ai];

    switch (L[li].model)
    {
      assert (pdgnew[ui] == 0);
    case HKY:
      if (assignmentoptions[JCMODEL] == 1)
      {
        temppdg = pdgnew[ui] =
          likelihoodJC (ci, li, G->uvals[0]);
      }
      else
      {
        temppdg = pdgnew[ui] =
          likelihoodHKY (ci, li, G->uvals[0], G->kappaval, -1, -1, -1, -1);
      }
      break;
    case INFINITESITES:
      temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      break;
    case STEPWISE:
      temppdg = 0;
      for (; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    case JOINT_IS_SW:
      temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      for (ai = 1; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    }
    pdgnewsum += temppdg;
  }
  tpw = pdgnewsum - pdgoldsum;
  assert (!ODD (ecd));
  assert (!ODD (ecu));
  ecd /= 2;
  ecu /= 2;
  integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight_t_RY,
                       &C[ci]->allpcalc, &holdallpcalc_t_RY, &holdt[0]);
  //printf ("pdg: %lf %lf (%lf), ", pdgnewsum, pdgoldsum, tpw);
  tpw += C[ci]->allpcalc.probg - holdallpcalc_t_RY.probg;
  //printf ("pg: %lf %lf (Total: %lf)", C[ci]->allpcalc.probg, holdallpcalc_t_RY.probg, tpw);

  metropolishastingsterm = beta[ci] * tpw + (ecd + emd) * log (t_d_hterm) +
    (ecu + emu) * log (t_u_hterm);
  //printf (" - term: %lf\n", metropolishastingsterm);

  //assert(metropolishastingsterm >= -1e200 && metropolishastingsterm < 1e200);
  U = log (uniform ());
  if (metropolishastingsterm >= 0.0 || metropolishastingsterm > U)
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].pdg = 0;
      for (ai = 0; ai < L[li].nlinked; ai++)
      {
        C[ci]->G[li].pdg_a[ai] = pdgnew[L[li].uii[ai]];
        C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
      }
      if (L[li].model == HKY)
      {
        storescalefactors (ci, li);
        copyfraclike (ci, li);
      }
    }
    C[ci]->allpcalc.pdg = pdgnewsum;
    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;

    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
    return 1;
  }
  else
  {
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_RY);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_RY);
    assert (pdgoldsum == C[ci]->allpcalc.pdg);
    C[ci]->tvals[timeperiod] = oldt;
    for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      storetreestats_all_loci (ci, 1);
      copy_treeinfo (&G->gweight, &holdgweight_t_RY[li]);
      for (i = 0; i < L[li].numlines; i++)
      {
        if (gtree[i].down != -1
            && skipflag[li][i] == 0 /* use with setedgeflag */ )
        {
          if (gtree[i].time <= newt && gtree[i].time > t_u)
          {
            assert (skipflag[li][i] == 0);
            gtree[i].time =
              beforesplit (timeperiod, newt, oldt, t_d, t_u, gtree[i].time);
            //cecu++;
          }

          else
          {
            if (gtree[i].time > newt && gtree[i].time < t_d)
            {
              assert (skipflag[li][i] == 0);
              gtree[i].time =
                aftersplit (timeperiod, newt, oldt, t_d, t_u, gtree[i].time);
              //cecl++;
            }
          }
          j = 0;
          while (gtree[i].mig[j].mt > -0.5)
          {
            if (gtree[i].mig[j].mt <= newt && gtree[i].mig[j].mt > t_u)
            {
              gtree[i].mig[j].mt =
                beforesplit (timeperiod, newt, oldt, t_d,
                             t_u, gtree[i].mig[j].mt);
              //cemu++;
            }
            else if (gtree[i].mig[j].mt > newt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplit (timeperiod, newt, oldt, t_d, t_u,
                            gtree[i].mig[j].mt);
              //ceml++;
            }
            j++;
          }
        }
      }
//        assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
    }
    /*    assert(ecu==cecu/2);
       assert(ecd==cecl/2);
       assert(emu==cemu);
       assert(emd==ceml); */
    for (li = 0; li < nloci; li++)
    {
      if (L[li].model == HKY)
        restorescalefactors (ci, li);
      /* have to reset the dlikeA values in the trees for stepwise model */
      if (L[li].model == STEPWISE)
        for (ai = 0; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      if (L[li].model == JOINT_IS_SW)
        for (ai = 1; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      // assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
    }
    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
    return 0;
  }
}                               /* changet_RY1 */
