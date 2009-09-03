/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS

#include "imamp.h"
#include "update_gtree_common.h"
#include "update_gtree.h"
#include "updateassignment.h"

/* this file holds stuff for changet_NW()
This is a t update modeled on the way Nielsen and Wakeley  (2001)did it as implemented in the midv program,  by moving a split time up or down
and adding or erasing migration events. 
check out mathematica file: newt_t_updating_8_21_08.nb 
*/

/* declarded in update_gtree.h */
/* extern int holddownA[MAXLINKED]; not used here */
/* extern int medgedrop; not used here */
/* extern int mrootdrop; not used here */
/* extern int rootmove; not used here */
/* extern double lmedgedrop; not used here */
/* extern double lmrootdrop; not used here */
/* extern double holdsisdlikeA[MAXLINKED]; not used here */
/* extern struct genealogy_weights holdgweight_updategenealogy;
extern struct genealogy_weights holdallgweight_updategenealogy;
extern struct probcalc holdallpcalc_updategenealogy; 
extern struct genealogy holdgtree; not used here */


/*********** local to this file  ***************/

/* struct migrationinfo_tNW  
holds  information about an edge that is required for the update 
a static array, large enough to hold the largest of genealogies in the data set, 
is set up the first time through 
*/

struct migrationinfo_tNW
{
  int upk;                      //- upedge state: -1 if top is in time split interval and state is unknown, otherwise the population #
  int dk;                       //- downedge state:  -1 if bottom is in time split interval and state is unknown, otherwise the population #
  int upa;                      // if upk > -1  this is the population state at the top of (but just below) the split interval, after the update
  int upb;                      // if upk > -1  this is the population state at the top of (but just below) the split interval, before the update
  int da;                       // if dk > -1 this is the population state just before the bottom of the split interval, after the update
  int db;                       // if dk > -1 this is the population state just before the bottom of the split interval, before the update
  double mtime;                 // time of edge that is in the interval between the new and old split times
  int mcount;                   // # of migrations in that time interval before the update
  int mnew;                     // # of migrations in that interval after the uptdate
  double mrate;                 // migration rate over that interval before udpate 
  double mrate_r;               // migration rate over that interval after udpate 
  int simbottomflag;            // 1 if bottom population state must be simulated, 0 otherwise
  int done;                     // 1 if edge has been simulated or does not need to be simulated
  int skip;                     // edge spans relevant time but not in relevant populations and can be skipped,  
  // skip has the same value as 'done' to begin with but unlike 'done' it does not ever get reset
  int linkto;                   // sister edge # if they should be simulated together,  -1 if not
} *minfo;


static struct genealogy *holdallgtree_t_NW;
static int largestsamp;
static int ci;
static struct genealogy *G;
static double oldt, newt, tu, td;

static struct genealogy_weights holdallgweight_t_NW;
static struct genealogy_weights holdgweight_t_NW[MAXLOCI];
static struct probcalc holdallpcalc_t_NW;

/******* local prototypes ********/
static void copy_all_gtree (int mode);
static double pathcondition2 (int issame, int moves, int popsless1);
static double getmprob_tNW (double mrate, int period, double mtimeavail,
                            int mnew, int curpop, int endpop,
                            int simbottomflag, int twotargetflag);
static int simmpath_tNW (int edge, struct edge *gtree, int numm, int lastm,
                         double timein, double upt, int period, int pop,
                         int constrainpop);
static int simmpath_tNW_2 (int edge, struct edge *gtree, int numm, int lastm,
                           double timein, double upt, int period, int pop,
                           int constrainpop1, int constrainpop2);
static int addmigration_tNW_2 (int li, struct edge *gtree, int edge, int sis,
                               double mrate, int period);
static void addmigration_tNW (struct edge *gtree, int i, double mrate,
                              int period, double uptime);
static double update_mig_tNW (int li, struct edge *gtree, int period);

// copy_all_gtree  makes a copy of all genealogies in a given chain 
// just like copy_gtree() but for all loci in a chain
// could set it up so it calls copy_gtree() but that function 
// uses a local holdgtree variable - could change this to cut down on code
// if mode==1  copy *G into holdgtree  if mode==0 copy holdgtree into *G
void
copy_all_gtree (int mode)
{
  struct edge *togtree, *fromgtree;
  struct genealogy *fromG, *toG;
  int li, i, j, ai;

  for (li = 0; li < nloci; li++)
  {
    if (mode)
    {
      toG = &holdallgtree_t_NW[li];
      fromG = &C[ci]->G[li];
    }
    else
    {
      toG = &C[ci]->G[li];
      fromG = &holdallgtree_t_NW[li];
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
        checkmig (j, &(togtree[i].mig), &(togtree[i].cmm));
        togtree[i].mig[j] = fromgtree[i].mig[j];
      }
      while (fromgtree[i].mig[j].mt > -0.5);
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
  }
}                               /* copy_all_gtree */

/* pathcondition2() similar to pathcondition() but used for the case when there are two target populations 
see new_t_updating.nb     uses recursion */
double
pathcondition2 (int issame, int moves, int popsless1)
{
  if ((popsless1 == 2 && moves >= 12) || (popsless1 > 2 && moves >= 6))
    return 2.0 / (popsless1 + 1);
  if (issame)
  {
    if (moves > 0)
      return (1.0 - 1.0 / popsless1) * pathcondition2 (0, moves - 1, popsless1) + (1.0 / popsless1) * pathcondition2 (1, moves - 1, popsless1);
    else
      return 1;
  }
  else
  {
    if (moves > 0)
      return (2.0 / popsless1) * pathcondition2 (1, moves - 2, popsless1) + (1.0 - 2.0 / popsless1) * pathcondition2 (0, moves - 1, popsless1);
    else
      return 0;
  }
}                               /*pathcondition2 */


/* calculate the probability of the migration event on an edge update for migration in the interval between new and old splittimes*/
/* closely follows mathematica notebook  'new_t_updating_8_10_08.nb' simbottomflag is 1 if the bottom population was simulated */
double
getmprob_tNW (double mrate, int period, double mtimeavail, int mnew,
              int curpop, int endpop, int simbottomflag, int twotargetflag)
{
  double tempp, n, d;
  double r, logmrate;
  int popc;

  assert (mtimeavail > 0);
  logmrate = log (mrate);
  popc = npops - period - 1;
  if (simbottomflag)
  {
    r = mtimeavail * mrate * popc;
    tempp = mnew * logmrate - r;
  }
  else
  {
    if (popc == 1)              // only 2 pops in period
    {
      assert (mtimeavail > 0);
      if (ODD (mnew))
        tempp = mnew * logmrate - mylogsinh (mrate * mtimeavail);
      else
        tempp = mnew * logmrate - mylogcosh (mrate * mtimeavail);
    }
    else                        // 3 or more pops in the period 
    {
      r = mtimeavail * mrate * popc;
      if (endpop == C[ci]->addpop[period + 1] && twotargetflag) // migration was simulated with two target populations
      {
        if (curpop != C[ci]->droppops[period + 1][0]
            && curpop != C[ci]->droppops[period + 1][1])
        {
          d = log ((1 - exp (-r)) * pathcondition2 (0, mnew, popc));
          n = mnew * logmrate - r;
        }
        else
        {
          n = mnew * logmrate - r;
          d = log (pathcondition2 (1, mnew, popc));
        }
      }
      else
      {
        if (curpop == endpop || endpop == C[ci]->poptree[curpop].down)
        {
          assert (mnew != 1);
          if (mnew == 0)
          {
            n = -r;
            d = log (1 - r * exp (-r));
          }
          else                  // >=2 
          {
            d = log ((1 - r * exp (-r)) * pathcondition (1, mnew, popc));
            n = mnew * logmrate - r;
          }
        }
        else
        {
          assert (mnew > 0);
          if (mnew == 1)
          {
            d = log ((1 - exp (-r)) / popc);
            n = logmrate - r;
          }
          else                  // >= 2
          {
            d = log ((1 - exp (-r)) * pathcondition (0, mnew, popc));
            n = mnew * logmrate - r;
          }
        }
      }
      tempp = n - d;
    }
  }
  assert (tempp > -1e200 && tempp < 1e200);
  return tempp;
}                               /* getmprob_tNW */

  /* simulate the migration path  -  very similar to simmpath() */
int
simmpath_tNW (int edge, struct edge *gtree, int numm, int lastm,
              double timein, double upt, int period, int pop,
              int constrainpop)
{
  int i, lastpop, startm;
  assert (numm > 0);
  startm = lastm + 1;
  lastm = lastm + numm;

  for (i = startm; i <= lastm; i++)
  {
    checkmig (i, &(gtree[edge].mig), &(gtree[edge].cmm));
    gtree[edge].mig[i].mt = upt + uniform () * timein;
  }
  gtree[edge].mig[i].mt = -1;
  if (numm > 1)
    hpsortmig (&gtree[edge].mig[startm] - 1, numm);
  lastpop = pop;
  if (constrainpop < 0)
  {
    for (i = startm; i <= lastm; i++)
    {
      gtree[edge].mig[i].mp =
        picktopop (lastpop, C[ci]->plist[period], npops - period);
      lastpop = gtree[edge].mig[i].mp;
    }
  }
  else
  {
    if (numm >= 2)
    {
      i = startm;
      while (i < lastm - 1)
      {
        gtree[edge].mig[i].mp =
          picktopop (lastpop, C[ci]->plist[period], npops - period);
        lastpop = gtree[edge].mig[i].mp;
        i++;
      }
      gtree[edge].mig[lastm - 1].mp =
        picktopop2 (lastpop, C[ci]->plist[period], npops - period,
                    constrainpop);
      gtree[edge].mig[lastm].mp = constrainpop;
    }
    else
    {
      if (numm == 1)
      {
        gtree[edge].mig[lastm].mp = constrainpop;
        lastpop = constrainpop;
      }
    }
  }
  return lastpop;
}                               /* simmpath_tNW */


/* very similar to simmpath_tNW()  used when there are a two populations that the migration path can end in */
int
simmpath_tNW_2 (int edge, struct edge *gtree, int numm, int lastm,
                double timein, double upt, int period, int pop,
                int constrainpop1, int constrainpop2)
{
  int i, lastpop, startm;

  assert (numm > 0);
  startm = lastm + 1;
  lastm = lastm + numm;

  for (i = startm; i <= lastm; i++)
  {
    checkmig (i, &(gtree[edge].mig), &(gtree[edge].cmm));
    gtree[edge].mig[i].mt = upt + uniform () * timein;
  }
  gtree[edge].mig[i].mt = -1;
  if (numm > 1)
    hpsortmig (&gtree[edge].mig[startm] - 1, numm);
  lastpop = pop;
  if (numm >= 1)
  {
    i = startm;
    while (i < lastm)
    {
      gtree[edge].mig[i].mp =
        picktopop (lastpop, C[ci]->plist[period], npops - period);
      lastpop = gtree[edge].mig[i].mp;
      i++;
    }
  }
  if (lastpop == constrainpop1)
    gtree[edge].mig[lastm].mp = constrainpop2;
  else
  {
    if (lastpop == constrainpop2)
    {
      gtree[edge].mig[lastm].mp = constrainpop1;
    }
    else
    {
      if (bitran ())
        gtree[edge].mig[lastm].mp = constrainpop2;
      else
        gtree[edge].mig[lastm].mp = constrainpop1;
    }
  }
  lastpop = gtree[edge].mig[lastm].mp;
  return lastpop;
}                               /* simmpath_tNW_2 */


/* addmigration_tNW_2 picks number of migrations on edge and its sister edge and simulates the path, returns population state where the two join*/
/* used for when an edge and its sister both need migration simulations */
/* returns the population state of the node at junction of edge and sis */
int
addmigration_tNW_2 (int li, struct edge *gtree, int edge, int sis,
                    double mrate, int period)
{
  int i, j, k, mnew, mdone;
  double uptime_edge, uptime_sis;
  int m_up_sis, m_up_edge, lastpopedge, lastpopsis;
  double r, bothtimes, t;
  double migtimes[ABSMIGMAX];   // very wasteful,  could use dynamic memory but would complicate and slow
  int dstate;

  // count migrations on edge
  uptime_edge = (edge < L[li].numgenes) ? 0 : gtree[gtree[edge].up[0]].time;
  m_up_edge = 0;                // becomes the position in migration array to add new migration events
  if (uptime_edge < tu)
    while (gtree[edge].mig[m_up_edge].mt > -0.5
           && gtree[edge].mig[m_up_edge].mt < tu)
    {
      m_up_edge++;
    }
  // count migrations on sis
  uptime_sis = (sis < L[li].numgenes) ? 0 : gtree[gtree[sis].up[0]].time;
  m_up_sis = 0;                 // becomes the position in migration array to add new migration events
  if (uptime_sis < tu)
    while (gtree[sis].mig[m_up_sis].mt > -0.5
           && gtree[sis].mig[m_up_sis].mt < tu)
    {
      m_up_sis++;
    }
  /* there are npops-period populations,  one less than this are migration targets */
  if (newt < td && period == lastperiodnumber)  // do not add migration
  {
    mnew = 0;
    gtree[edge].mig[m_up_edge].mt = gtree[sis].mig[m_up_sis].mt = -1;
    // this was a bug, fixed 2/16/09:  dstate = C[ci]->addpop[period + 1];   
    dstate = C[ci]->addpop[period];
  }
  else
  {
    bothtimes = minfo[edge].mtime + minfo[sis].mtime;
    r = bothtimes * mrate * (npops - period - 1);
    if (npops - period == 2)
    {
      mnew = poisson (r, minfo[edge].upa != minfo[sis].upa);    /* odd if different pops,  even if same */
    }
    else
    {
      if (minfo[edge].upa == minfo[sis].upa)    //cannot be just one migration event
        mnew = poisson (r, 3);
      else                      // cannot be zero migration events 
        mnew = poisson (r, 2);
    }
    if (mnew > 0)
    {
      for (i = 0; i < mnew; i++)
      {
        migtimes[i] = uniform () * bothtimes;
      }
      if (mnew > 1)
        hpsortreg (migtimes - 1, mnew);
      uptime_edge = DMAX (tu, uptime_edge);
      j = m_up_edge;
      i = 0;
      t = uptime_edge + migtimes[i];
      lastpopedge = minfo[edge].upa;
      while (t < gtree[edge].time && i < mnew)  //JH  I think this is fully initiatlized, Sang Chul agreed.
      {
        checkmig (j, &(gtree[edge].mig), &(gtree[edge].cmm));
        gtree[edge].mig[j].mt = t;
        if (npops - period == 2)
        {
          gtree[edge].mig[j].mp = lastpopedge = picktopop2 (lastpopedge, C[ci]->plist[period], 2, lastpopedge);
        }
        else                    // more than 2 populations in the period
        {
          if (i == mnew - 2)
            lastpopedge =
              picktopop2 (lastpopedge, C[ci]->plist[period],
                          npops - period, minfo[sis].upa);
          else
          {
            if (i == mnew - 1)
              lastpopedge = minfo[sis].upa;
            else
              lastpopedge = picktopop (lastpopedge, C[ci]->plist[period], npops - period);
          }
          gtree[edge].mig[j].mp = lastpopedge;
        }
        i++;
        t = uptime_edge + migtimes[i];
        j++;
      }
      gtree[edge].mig[j].mt = -1;

      // now do the sister branch
      lastpopsis = minfo[sis].upa;
      j = m_up_sis;
      if (i < mnew)
      {
        mdone = i;
        uptime_sis = DMAX (tu, uptime_sis);
        // can do this either way,  from the last migration event in migtimes and move backwards,  or forward from the next one in migtimes
        for (i = mnew - 1, k = 0; i >= mdone; k++, i--, j++)    // i starts at position of last migtimes,  there are mnew - mdone migrations to do 
          // for (k=0; i < mnew; k++,j++, i++)   // starts at next position in migtimes
        {
          t = uptime_sis + (bothtimes - migtimes[i]);   // if starting at last position
          // t = migtimes[i] + uptime_sis -  (gtree[sis].time - uptime_edge);  // if starting at next position
          assert (t <= gtree[sis].time);        // should not need to make this <=  
          assert (t > uptime_sis);

          checkmig (j, &(gtree[sis].mig), &(gtree[sis].cmm));
          gtree[sis].mig[j].mt = t;
          if (npops - period == 2)
          {
            gtree[sis].mig[j].mp = lastpopsis = picktopop (lastpopsis, C[ci]->plist[period], 2);
          }
          else
          {
            if (k == mnew - mdone - 2)
              lastpopsis = picktopop2 (lastpopsis, C[ci]->plist[period], npops - period, lastpopedge);
            else
            {
              if (k == mnew - mdone - 1)
                lastpopsis = lastpopedge;
              else
                lastpopsis = picktopop (lastpopsis, C[ci]->plist[period], npops - period);
            }
            gtree[sis].mig[j].mp = lastpopsis;
          }
        }
      }
      gtree[sis].mig[j].mt = -1;
      dstate = lastpopedge;
      assert (lastpopsis == lastpopedge);
    }
    else
    {
      gtree[edge].mig[m_up_edge].mt = gtree[sis].mig[m_up_sis].mt = -1;
      dstate = minfo[edge].upa;
    }
  }
  minfo[edge].mnew = mnew;
  assert (dstate >= 0 && dstate < numtreepops);
  return dstate;
}                               /* addmigration_tNW_2 */

/*addmigration_tNW picks number of migrations and calls simmpath_tNW */
/* read notes in new_t_updating_8_21_08.nb */
void
addmigration_tNW (struct edge *gtree, int i, double mrate, int period,
                  double uptime)
{
  int k, j, jj, mnew, lastabovem, lastbelowm, numskip, numheld;
  double r;
  struct migstruct holdmig[ABSMIGMAX];  // very wasteful of space  - could use dynamic memory and checkmig() ? 
  int simmpathtype;

  /* determine position in mig array of the first migration event in the relevant period */
  j = 0;
  if (uptime < tu)
    while (gtree[i].mig[j].mt > -0.5 && gtree[i].mig[j].mt < tu)
    {
      j++;
    }
  lastabovem = j - 1;           // position to start simulating

  /* count how many migrations must be erased for the update */
  jj = j;
  while (gtree[i].mig[jj].mt > -0.5 && gtree[i].mig[jj].mt < td)
  {
    jj++;
  }
  numskip = jj - j;

  /* count and save the migrations after the relevant period  */
  if (gtree[i].mig[jj].mt < 0)
  {
    lastbelowm = -1;            // no migration events to save after interval
    numheld = 0;
  }
  else
  {
    lastbelowm = jj;
    k = 0;
    while (gtree[i].mig[jj].mt > -0.5)
    {
      holdmig[k] = gtree[i].mig[jj];
      jj++;
      k++;
    }
    holdmig[k].mt = -1;
    numheld = k;
  }

  /* there are npops-period populations,  one less than this are migration targets */
  if (newt < td && period == lastperiodnumber)  // do not add migration
  {
    minfo[i].mnew = 0;
    gtree[i].mig[lastabovem + 1].mt = -1;
  }
  else
  {
    r = minfo[i].mtime * mrate * (npops - period - 1);
    simmpathtype = 1;           // sets which simulation path function gets called 
    if (minfo[i].simbottomflag) // can do migration to any available pop
    {
      mnew = poisson (r, -1);
    }
    else                        // bottom pop is known,  either as a coalescent before td  or as a population state after td
    {
      if (newt > tu && (minfo[i].dk == C[ci]->addpop[period + 1]))      // target is one of the populations that join
      {
        simmpathtype = 0;
        if (minfo[i].upa == C[ci]->droppops[period + 1][0]
            || minfo[i].upa == C[ci]->droppops[period + 1][1])
          mnew = poisson (r, -1);
        else
          mnew = poisson (r, 2);
      }
      else
      {
        if (npops - period == 2)
        {
          mnew = poisson (r, minfo[i].upa != minfo[i].dk);      /* odd if different pops,  even if same */
        }
        else
        {
          if (minfo[i].upa == minfo[i].dk || (minfo[i].dk == C[ci]->addpop[period + 1]))        //cannot be just one migration event
            mnew = poisson (r, 3);
          else                  // cannot be zero migration events 
            mnew = poisson (r, 2);
        }
      }
    }
    if (mnew > 0)
    {
      if (simmpathtype)
      {
        if (minfo[i].simbottomflag)
          minfo[i].da =
            simmpath_tNW (i, gtree, mnew, lastabovem, minfo[i].mtime,
                          DMAX (tu, uptime), period, minfo[i].upa, -1);
        else
          minfo[i].da =
            simmpath_tNW (i, gtree, mnew, lastabovem, minfo[i].mtime,
                          DMAX (tu, uptime), period, minfo[i].upa,
                          minfo[i].dk);
      }
      else
        minfo[i].da =
          simmpath_tNW_2 (i, gtree, mnew, lastabovem, minfo[i].mtime,
                          DMAX (tu, uptime), period, minfo[i].upa,
                          C[ci]->droppops[period + 1][0],
                          C[ci]->droppops[period + 1][1]);

      if (C[ci]->poptree[minfo[i].dk].e <= period
          && C[ci]->poptree[minfo[i].dk].e != -1)
        minfo[i].dk = C[ci]->poptree[minfo[i].dk].down;

      /* put back stored migration events */
      if (numheld > 0)
      {
        j = lastabovem + mnew + 1;
        for (k = 0; k < numheld; k++, j++)
        {
          checkmig (j, &(gtree[i].mig), &(gtree[i].cmm));
          gtree[i].mig[j] = holdmig[k];
        }
        gtree[i].mig[j].mt = -1;
      }
    }
    else
    {
      assert (minfo[i].dk >= 0);
      //if (minfo[i].dk < 0)
      if (minfo[i].dk == -1)
      {
        minfo[i].dk = minfo[i].upa;
      }
      if (numskip > 0)
      {
        /* put back stored migration events */
        if (numheld > 0)
        {
          for (j = lastabovem + 1, k = 0; k < numheld; k++, j++)
            gtree[i].mig[j] = holdmig[k];
          gtree[i].mig[j].mt = -1;
        }
        else
        {
          gtree[i].mig[lastabovem + 1].mt = -1;
        }
      }
    }
    minfo[i].mnew = mnew;
  }
}                               /* addmigration_tNW */

/* updating migration events during a changet_NW update */
/* this is a bug function with multiple parts */

/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre tree - e.g. help to avoid promoting
correlations and improve mixing  

added mrate and mrate_r  to struct migrationinfo_tNW   
*/
double
update_mig_tNW (int li, struct edge *gtree, int period)
{
  double uptime;                //, tu, td;
  int i, j, k, kk, sis, mpart;
  double tlengthpart;
  int n;
  double mproposedenom, mproposenum;
  int topopf, topopr;
  int period_a, period_b;       // the period the oldt/newt time interval falls in after, and before, the update

  /* tlengthpart and mpart will change as migration is update, 
     they are needed for calculating the migration rate for the reverse update */
  tlengthpart = G->tlength;
  mpart = G->mignum;
  //mrate = calcmrate (G->mignum, G->tlength);
  if (newt < oldt)
  {
    tu = newt;
    td = oldt;
    period_a = period + 1;
    period_b = period;
  }
  else
  {
    tu = oldt;
    td = newt;
    period_a = period;
    period_b = period + 1;
  }

    /*********** Setup minfo - indicate which branches need attention  ***********/
  // read about struct migrationinfo_tNW at the top of this file
  for (i = 0; i < L[li].numlines; i++)
  {
    uptime = (i < L[li].numgenes) ? 0 : gtree[gtree[i].up[0]].time;
    if ((gtree[i].time > tu && uptime <= td) && i != G->root)
    {
      minfo[i].skip = minfo[i].done = 0;
      minfo[i].linkto = -1;
      if (uptime < tu)
      {
        minfo[i].upk = nowedgepop (ci, &gtree[i], tu);
        if (minfo[i].upk == C[ci]->droppops[period + 1][0]
            || minfo[i].upk == C[ci]->droppops[period + 1][1])
        {
          if (newt < oldt)
          {
            minfo[i].upb = minfo[i].upk;
            minfo[i].upa = C[ci]->addpop[period + 1];
          }
          else
          {
            minfo[i].upb = C[ci]->addpop[period + 1];
            minfo[i].upa = minfo[i].upk;
          }
        }
        else
        {
          minfo[i].upb = minfo[i].upa = minfo[i].upk;
        }
      }
      else
      {
        minfo[i].upb = nowedgepop (ci, &gtree[i], uptime);
        minfo[i].upa = minfo[i].upk = -1;
      }
      if (gtree[i].time > td)
      {
        /* nowedgepop() returns the population an edge is in at a particular time
           the use of 1+DBL_EPSILON is to get the population immediately after a split time event
           without it nowedgepop() would return the population right at the split time event */
        minfo[i].dk = nowedgepop (ci, &gtree[i], td * (1 + DBL_EPSILON));       // some risk that this won't work for all values of td
        if (minfo[i].dk == C[ci]->addpop[period + 1])
        {
          if (newt < oldt)
          {
            minfo[i].db =
              nowedgepop (ci, &gtree[i], oldt / (1 + DBL_EPSILON));
            minfo[i].da = minfo[i].dk;
          }
          else
          {
            minfo[i].db = minfo[i].dk;
            minfo[i].da = -1;
          }
        }
        else
        {
          minfo[i].da = minfo[i].db = minfo[i].dk;
        }

      }
      else
      {
        minfo[i].da = minfo[i].db = minfo[i].dk = -1;
      }
      minfo[i].mtime = DMIN (td, gtree[i].time) - DMAX (tu, uptime);
      if (period == lastperiodnumber - 1)       // set amount of length of tree in which migration can occur has been altered. 
      {
        if (newt < td)
          tlengthpart -= minfo[i].mtime;
        else
          tlengthpart += minfo[i].mtime;
      }
      j = 0;
      k = 0;
      while (gtree[i].mig[j].mt > -0.5 && gtree[i].mig[j].mt < td)
      {
        if (gtree[i].mig[j].mt > tu)
          k++;
        j++;
      }
      minfo[i].mcount = k;
      /* skip edges that never enter the splitpops or addpop */
      if ((minfo[i].upk >= 0) && (minfo[i].dk >= 0) &&
          (minfo[i].upk != C[ci]->droppops[period + 1][0])
          && (minfo[i].upk != C[ci]->droppops[period + 1][1])
          && (minfo[i].dk != C[ci]->addpop[period + 1]))
      {
        k = 0;
        while (gtree[i].mig[k].mt > 0.0 && gtree[i].mig[k].mt < tu)
          k++;
        kk = k;
        while (kk < minfo[i].mcount && gtree[i].mig[kk].mt > 0.0 &&
               gtree[i].mig[kk].mp != C[ci]->droppops[period + 1][0] &&
               gtree[i].mig[kk].mp != C[ci]->droppops[period + 1][1] &&
               gtree[i].mig[kk].mp != C[ci]->addpop[period + 1])
          kk++;
        if ((kk - k) == minfo[i].mcount)
          minfo[i].skip = minfo[i].done = 1;
      }
    }
    else
    {
      minfo[i].skip = minfo[i].done = 1;
    }
  }

    /********* Simulate Migrations ***********/
  /* now that minfo is all set,  go through and add migrations to branches that need it
     if the population state of the top of a branch is not known,  it is skipped until after that state
     is determined 
     Note the use of minfo[i].link and that calls for simulating two branches and and state of the bottom 
     of the two branches */

  do                            // loop thru repeatedly to add migrations as needed
  {
    n = 0;
    for (i = 0; i < L[li].numlines; i++)
    {
      if (minfo[i].done == 0)
      {
        if (minfo[i].upk >= 0 && minfo[i].dk >= 0)      // both up and down pop states known
        {
          uptime = (i < L[li].numgenes) ? 0 : gtree[gtree[i].up[0]].time;
          minfo[i].mrate = calcmrate (minfo[i].mcount, minfo[i].mtime);
          if (period == (lastperiodnumber - 1) && minfo[i].dk == C[ci]->addpop[period + 1])     // is this last condition redundant?
          {
            assert (minfo[i].dk == C[ci]->addpop[period + 1]);
            minfo[i].simbottomflag = 1;
            addmigration_tNW (gtree, i, minfo[i].mrate, period_a, uptime);
          }
          else
          {
            minfo[i].simbottomflag = 0;
            addmigration_tNW (gtree, i, minfo[i].mrate, period_a, uptime);
          }
          minfo[i].mrate_r = calcmrate (minfo[i].mnew, minfo[i].mtime);
          mpart += minfo[i].mnew - minfo[i].mcount;
          minfo[i].done = 1;
        }
        else
        {
          if ((sis = gtree[gtree[i].down].up[0]) == i)
            sis = gtree[gtree[i].down].up[1];
          if (minfo[i].upk >= 0 && minfo[sis].upk >= 0) // upper states known,  bottom state not known
          {
            assert (minfo[i].dk == -1);
            assert (!minfo[sis].done);
            n++;
            minfo[i].simbottomflag = minfo[sis].simbottomflag = 0;
            minfo[i].linkto = sis;
            minfo[sis].linkto = i;
            assert (minfo[i].da == -1 && minfo[sis].da == -1);
            minfo[i].mrate = minfo[sis].mrate = calcmrate (minfo[i].mcount + minfo[sis].mcount, minfo[i].mtime + minfo[sis].mtime);
            minfo[i].da = minfo[sis].da = addmigration_tNW_2 (li, gtree, i, sis, minfo[i].mrate, period_a);

            if (minfo[i].dk == -1)
            {
              minfo[i].dk = minfo[i].da;
              while /*if */ (C[ci]->poptree[minfo[i].dk].e <= period
                             && C[ci]->poptree[minfo[i].dk].e != -1)
                minfo[i].dk = C[ci]->poptree[minfo[i].dk].down;
            }
            if (minfo[sis].dk == -1)
            {
              minfo[sis].dk = minfo[sis].da;
              while /*if */ (C[ci]->poptree[minfo[sis].dk].e <= period
                             && C[ci]->poptree[minfo[sis].dk].e != -1)
                minfo[sis].dk = C[ci]->poptree[minfo[sis].dk].down;
            }
            minfo[i].mrate_r = minfo[sis].mrate_r = calcmrate (minfo[i].mnew, minfo[i].mtime + minfo[sis].mtime);

            minfo[sis].mnew = minfo[i].mnew;
            minfo[i].mtime = minfo[sis].mtime = minfo[i].mtime + minfo[sis].mtime;
            minfo[i].mcount = minfo[sis].mcount = minfo[i].mcount + minfo[sis].mcount;
            mpart += minfo[i].mnew - minfo[i].mcount;
            assert (mpart >= 0);
            assert (gtree[i].down == G->root
                    || minfo[gtree[i].down].upk == -1);
            minfo[gtree[i].down].upk = minfo[gtree[i].down].upa = gtree[gtree[i].down].pop = minfo[i].dk;       // this might not work for gtree[gtree[i].down].pop for > 2 sampled populations 
            minfo[i].done = minfo[sis].done = 1;
          }
        }
      }
    }
  }
  while (n > 0);

    /******* Calculate the Hastings ratio *******/

  /* determine the migration rate for the reverse update */
  G->mignum = mpart;
  G->tlength = tlengthpart;
  // mrate_r = calcmrate (G->mignum, G->tlength);
  mproposedenom = mproposenum = 0;
  for (i = 0; i < L[li].numlines; i++)  // calculate forward and reverse hastings terms 
  {
    uptime = (i < L[li].numgenes) ? 0 : gtree[gtree[i].up[0]].time;
    /* calculate the probability for edges that span the interval and are not the root
       If an edge is linked to another,  only do this calculation for the lower numbered edge */
    if (minfo[i].skip == 0 && (gtree[i].time > tu && uptime <= td)
        && i != G->root && (minfo[i].linkto < 0 || i < minfo[i].linkto))
    {
      if (i < minfo[i].linkto)
      {
        topopf = minfo[minfo[i].linkto].upb;
        topopr = minfo[minfo[i].linkto].upa;
      }
      else
      {
        topopf = topopr = minfo[i].dk;
      }

      if (newt < td)            // split time is moving up, do not add migration for period prior to last split,  in forward direction
      {
        // forward
        if (period < lastperiodnumber - 1)
          mproposedenom +=
            getmprob_tNW (minfo[i].mrate, period_a, minfo[i].mtime,
                          minfo[i].mnew, minfo[i].upa, topopr,
                          minfo[i].simbottomflag, 0);
        // reverse
        mproposenum +=
          getmprob_tNW (minfo[i].mrate_r, period_b, minfo[i].mtime,
                        minfo[i].mcount, minfo[i].upb, topopf,
                        minfo[i].simbottomflag, 1);
      }
      else                      // do not add migration for period prior to last split, in reverse direction
      {
        //reverse
        if (period < lastperiodnumber - 1)
          mproposenum +=
            getmprob_tNW (minfo[i].mrate_r, period_b, minfo[i].mtime,
                          minfo[i].mcount, minfo[i].upb, topopf,
                          minfo[i].simbottomflag, 0);
        //forward
        mproposedenom +=
          getmprob_tNW (minfo[i].mrate, period_a, minfo[i].mtime,
                        minfo[i].mnew, minfo[i].upa, topopr,
                        minfo[i].simbottomflag, 1);
      }


      assert (mproposedenom > -1e200 && mproposedenom < 1e200
              && mproposenum > -1e200 && mproposenum < 1e200);
      i = i;
    }
  }
  return mproposenum - mproposedenom;
}                               /* update_mig_tNW */

/**************************************/
/*  GLOBAL FUNCTIONS in update_t_NW.c */
/**************************************/


void
init_t_NW (void)
{
  int li, j;
  init_genealogy_weights (&holdallgweight_t_NW);
  holdallgtree_t_NW = calloc ((size_t) nloci, sizeof (struct locus));
  for (li = 0; li < nloci; li++)
  {
    init_genealogy_weights (&holdgweight_t_NW[li]);
    init_holdgtree (&holdallgtree_t_NW[li], L[li].numgenes);
  }
  init_probcalc (&holdallpcalc_t_NW);
  for (largestsamp = 0, j = 0; j < nloci; j++)
    if (largestsamp < L[j].numlines)
      largestsamp = L[j].numlines;
  minfo = malloc ((2 * largestsamp - 1) * sizeof (struct migrationinfo_tNW));
  NW_t_upinf = calloc ((size_t) npops - 1, sizeof (struct update_rate_calc));
}                               /* init_t_NW */

void
free_t_NW (void)
{
  int li;
  free_genealogy_weights (&holdallgweight_t_NW);
  for (li = 0; li < nloci; li++)
  {
    free_genealogy_weights (&holdgweight_t_NW[li]);
    free_holdgtree (&holdallgtree_t_NW[li], L[li].numgenes);
  }
  XFREE (holdallgtree_t_NW);
  free_probcalc (&holdallpcalc_t_NW);
  XFREE (minfo);
  XFREE (NW_t_upinf);
}                               /* init_t_NW */


/* changet_NW
   this is modelled on the update of t in Nielsen and Wakeley (2001) 
   does nothing whatsover to branch lengths  so no change in pdg
   
   This update applies to all genealogies in the chain 
   */

/*  The update is done to the actual genealogies
    Could improve speed by doing the update to a copy of the genealogies, 
    This is because most updates are rejected, so if the copy is being changed and it is rejected 
    then there would be no need to restore it 

    to make this change would have to change a few lines of code (below)
    also have to pass treeweight a genealogy pointer and not just the genealogy number
    there will no doubt be other stuff
    spent a couple hours trying this on 8/28/08,  but wasn't working so gave up.  maybe try later

        
        Current Setup -  work on the original:
    --------------------------------------

        copy_all_gtree(1);  -  copy each C[ci]->G  to holdallgtree_t_NW
        copy_treeinfo (&holdallgweight_t_NW, &C[ci]->allgweight);  - copy allgweight to holdallgweight_t_NW - then empty allgweight
        copy_probcalc (&holdallpcalc_t_NW, &C[ci]->allpcalc);  - copy allpcalc to holdallpcalc_t_NW  (stuff for integrating)
        setzero_genealogy_weights (&C[ci]->allgweight); - empty allgweight

        
        going thru the loop:
            copy_treeinfo (&holdgweight_t_NW[li], &G->gweight);  - copy gweight to holdgweight_t_NW for each locus
            update C[ci]->G
            setzero_genealogy_weights (&G->gweight);   -  set gweight to zero
            treeweight (ci, li);   -   calculate weights
            sum_treeinfo (&C[ci]->allgweight, &G->gweight);  -     sum allgweight
        
        integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc); - integrate and reset allpcalc 

        if accept, 
            everything is already in gweight, allgweight and C[ci]->G 
        if reject:
            copy_all_gtree(0);  -  copy holdallgtree_t_NW back to each C[ci]->G
            copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_NW);  -  copy holdallgweight_t_NW back to allgweight
            copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_NW);   -  copy holdallpcalc_t_NW back to allpcalc
            for (li = 0; li < nloci; li++)   - loop thru loci
            {
                copy_treeinfo (&C[ci]->G[li].gweight, &holdgweight_t_NW[li]);  - copy each locus's from holdgweight_t_NW  back to gweight
            }

       If we reverse it :
       ------------------------------------
        copy_all_gtree(1);  -  copy each C[ci]->G  to holdallgtree_t_NW
        setzero_genealogy_weights (&holdallgweight_t_NW); - empty holdallgweight_t_NW


        going thru the loop:
            update holdgallgtree[li]
            setzero_genealogy_weights (&holdgweight_t_NW[li]);   -  set holdgweight_t_NW to zero
            local_treeweight (ci, li, holdallgtree_t_NW[li]);   -   calculate weights
            sum_treeinfo (&holdallgweight_t_NW, &holdgweight_t_NW[li]);  -     sum allgweight

        integrate_tree_prob (ci, &holdallgweight_t_NW, &holdallpcalc_t_NW); - integrate and reset allpcalc 

        if accept, 
            copy_all_gtree(0);  -  copy holdallgtree_t_NW to each C[ci]->G
            copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_NW);  -  copy holdallgweight_t_NW into allgweight
            copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_NW);   -  copy holdallpcalc_t_NW to allpcalc
            for (li = 0; li < nloci; li++)   - loop thru loci
            {
                copy_treeinfo (&C[ci]->G[li].gweight, &holdgweight_t_NW[li]);  - copy each locus's from holdgweight_t_NW  back to gweight
            }

        if reject:
            don't do anything because nothing in C[ci]->G, allgweight or allpcalc was changed.
 */

/* let u refer to the more recent time  and d to the older time  */

/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre tree - e.g. help to avoid promoting
correlations and improve mixing  */
/* only works for nonzero migration priors */
int
changet_NW (int chain, int timeperiod)
{

  double metropolishastingsterm, tpw;   //newt, oldt, tpw;
  int li, i;
  double U;
  struct edge *gtree;
  double t_d, t_u, t_u_prior, t_d_prior;
  double migweight;
  double holdt[MAXPERIODS];

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogy (chain);
  }

#ifdef IMPROFILE
  startclock (&clock_changet_NW);
#endif

  ci = chain;
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];

  /* select a new time */
  t_u = (timeperiod == 0) ? 0 : C[ci]->tvals[timeperiod - 1];
  t_d = (timeperiod == (lastperiodnumber - 1)) ? TIMEMAX : C[ci]->tvals[timeperiod + 1];
  t_d_prior = DMIN (T[timeperiod].pr.max, t_d);
  t_u_prior = DMAX (T[timeperiod].pr.min, t_u);
  oldt = C[ci]->tvals[timeperiod];
  newt = getnewt (timeperiod, t_u_prior, t_d_prior, oldt, 0);
  assert (newt < T[timeperiod].pr.max);
  /*if (timeperiod == (lastperiodnumber - 1) && modeloptions[SPLITTINGRATEPARAMETER])
     twin = (t_d_prior - t_u_prior) / (20*log(nloci+1) * npops);  
     else
     twin = (t_d_prior - t_u_prior) / (3*log(nloci+1) * (npops - timeperiod) );   */
  /* store stuff and prepare for adding to storing allgweight */
  copy_all_gtree (1);
  copy_treeinfo (&holdallgweight_t_NW, &C[ci]->allgweight);
  copy_probcalc (&holdallpcalc_t_NW, &C[ci]->allpcalc);
/*
    copy_treeinfo (&C[ci]->savedallgtinfo, &C[ci]->allgtinfo);
    copy_probcalc (&C[ci]->savedallpcalc, &C[ci]->allpcalc);
*/
  setzero_genealogy_weights (&C[ci]->allgweight);

  /* initialize migration hastings term and loop thru the loci */
  migweight = 0;
  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;
    copy_treeinfo (&holdgweight_t_NW[li], &G->gweight);
/*
        copy_treeinfo (&G->savedgtinfo, &G->gtinfo);
*/

    /* if the root of the genealogy is younger than oldt and newt no updating of this genealogy is needed */
    if ((newt > oldt && G->roottime > oldt)
        || (newt < oldt && G->roottime > newt))
    {
      gtree = G->gtree;
      migweight += update_mig_tNW (li, gtree, timeperiod);
      C[ci]->tvals[timeperiod] = newt;  // reset for treeweight() calculations
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
        C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;
      setzero_genealogy_weights (&G->gweight);
//if (step== 5)         gtreeprint(ci,0, step);
      treeweight (ci, li);
      C[ci]->tvals[timeperiod] = oldt;  // put old value back back for now because update_mig_tNW() depends on old value 
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time = oldt;
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = oldt;
    }
    sum_treeinfo (&C[ci]->allgweight, &G->gweight);
//      assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
  }

  /* calculate the prior for the new genealogy - needed even if genealogy not changed */
  C[ci]->tvals[timeperiod] = newt;      // reset for integrate calculations
  C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time = newt;
  C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;
  //integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
  integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight_t_NW,
                       &C[ci]->allpcalc, &holdallpcalc_t_NW, &holdt[0]);

  /* calculate the MH term - depends ratio of prior probabilities and hastings term for simulated migration events */
  /* does not depend on P(D|G)  because branch lengths are not changed */
  tpw = C[ci]->allpcalc.probg - holdallpcalc_t_NW.probg;
  metropolishastingsterm = exp (beta[ci] * tpw + migweight);
  U = uniform ();
  if (metropolishastingsterm >= 1.0 || metropolishastingsterm > U)
  {
    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (chain);
    }
    return 1;
  }
  else
  {
    /* copy stuff back from where it was held */
    C[ci]->tvals[timeperiod] = oldt;
    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = oldt;
    copy_all_gtree (0);
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_NW);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_NW);
    for (li = 0; li < nloci; li++)
    {
      copy_treeinfo (&C[ci]->G[li].gweight, &holdgweight_t_NW[li]);
      //      assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
    }
    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (chain);
    }
    return 0;
  }
}                               /* changet_NW  - after Nielsen and Wakeley (2001) */
