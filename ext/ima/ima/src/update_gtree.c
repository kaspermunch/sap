/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"
#include "update_gtree_common.h"
#include "update_gtree.h"
#include "updateassignment.h"

/*********** LOCAL STUFF **********/

struct edgemiginfo oldedgemig;
struct edgemiginfo oldsismig;
struct edgemiginfo newedgemig;
struct edgemiginfo newsismig;

//struct edge *copyedge;

/* declarded in update_gtree.h */
/* not used here
extern int holddownA[MAXLINKED];
extern int medgedrop;
extern double lmedgedrop;
extern double holdsisdlikeA[MAXLINKED];
extern struct genealogy holdgtree;
*/
static int mrootdrop;
static double lmrootdrop;
static struct genealogy_weights holdgweight_updategenealogy;
static struct genealogy_weights holdallgweight_updategenealogy;
static struct probcalc holdallpcalc_updategenealogy;
int rootmove;                   /* used in update_gtree.c and update_gtree_covar.c */


/* find the time when two populations join */
double
findjointime (int ci, int slidepop, int sispop, double edgeuptime,
              double sisuptime)
{
  int edgeperiod, sisperiod;
  double jointime;
  edgeperiod = findperiod (ci, edgeuptime);
  sisperiod = findperiod (ci, sisuptime);
  while (edgeperiod < sisperiod)

  {
    edgeperiod++;
    if (slidepop == C[ci]->droppops[edgeperiod][0]
        || slidepop == C[ci]->droppops[edgeperiod][1])
      slidepop = C[ci]->poptree[slidepop].down;
  }

  while (sisperiod < edgeperiod)
  {
    sisperiod++;
    if (sispop == C[ci]->droppops[sisperiod][0]
        || sispop == C[ci]->droppops[sisperiod][1])
      sispop = C[ci]->poptree[sispop].down;
  }

  // at this point edgeperiod == sisperiod 
  while (slidepop != sispop)
  {
    edgeperiod++;
    if (slidepop == C[ci]->droppops[edgeperiod][0]
        || slidepop == C[ci]->droppops[edgeperiod][1])
      slidepop = C[ci]->poptree[slidepop].down;
    if (sispop == C[ci]->droppops[edgeperiod][0]
        || sispop == C[ci]->droppops[edgeperiod][1])
      sispop = C[ci]->poptree[sispop].down;
  }

  if (edgeperiod == 0)
    jointime = 0;
  else
    jointime = C[ci]->tvals[edgeperiod - 1];
  return jointime;
}                               //findjointime

void
slider_nomigration (int ci, int li, int slidingedge, int *sis,
                    double *timepoint, double *slidedist)
/* this is very similar to slider, but with a an extra if/else on the slides up for the case of zero migration */
/* timepoint points at C[ci]->G[li].gtree[*slidingedge].time and is the current position of the sliding point, slidedist is the distance it must move 
with multiple populations,  and no migration, the sliding edge can only be in the population it started in, or an ancestral population
so the same goes for the point on the edge on which the slide is currently at. 
At the beginning of a slide the point is necessarily valid  - All down slides are valid for any distance
the upper limit at any point is the 
maximum of the top of the sliding edge and the time at which the sliding edge and the sister edge are in different populations  
if distance is negative  move up
	determine the upper limit ( MAX(top of edge, split time of sis and edge,  time of node of sis) )
	if upper limit is at a node,  pick left or right  at random
		call slider and continue up
	else reflect  (switch sign on remaining distance)
		call slider and and move down
else  move down
	if reach a node,  pick down or up at random
		if down,  call slider and continue down
		if up,  switch sign on remaining distance
			call slider and move up
kinds of upper limits
 top of sliding edge
 top of sister edge (a node) 
 beginning of period when the population of the edge and the sister branch come together into the same population 
*/
{
  double edgeuptime, sisuptime, popjointime;
  struct edge *gtree = C[ci]->G[li].gtree;
  int slidepop, sispop;
  if (*slidedist < 0)
  {
    /* go up */
    *slidedist = -*slidedist;
    if (slidingedge < L[li].numgenes)
      edgeuptime = 0;
    else
      edgeuptime = gtree[gtree[slidingedge].up[0]].time;

    if (*sis < L[li].numgenes)
      sisuptime = 0;
    else
      sisuptime = gtree[gtree[*sis].up[0]].time;

    assert (*timepoint >= edgeuptime);
    slidepop = gtree[slidingedge].pop;
    sispop = gtree[*sis].pop;
    if (slidepop != sispop)
      popjointime = findjointime (ci, slidepop, sispop, edgeuptime, sisuptime);
    else
      popjointime = 0;

    if (popjointime > edgeuptime && popjointime > sisuptime)
    {
      if (*slidedist < *timepoint - popjointime)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > popjointime);
        return;
      }
      else
      {

        /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
        *slidedist -= *timepoint - popjointime;
        *timepoint = popjointime;
        assert (*slidedist > 0);
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
    else
    {
      if (sisuptime == 0 || edgeuptime >= sisuptime)
      {
        if (*slidedist < *timepoint - edgeuptime)

          /* slide up and stop,  sis remains the same */
        {
          *timepoint -= *slidedist;
          *slidedist = 0;
          assert (*timepoint > edgeuptime);
          return;
        }
        else
        {
          /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
          *slidedist -= *timepoint - edgeuptime;
          *timepoint = edgeuptime;
          assert (*slidedist > 0);
          slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
          return;
        }
      }
      else
      {
        /* edgeuptime is less than sis up time, and thus slidingedge can reach a node */
        if (*slidedist < *timepoint - sisuptime)

          /* slide up and stop,  sis remains the same */
        {
          *timepoint -= *slidedist;
          *slidedist = 0;
          assert (*timepoint > sisuptime);
          return;
        }
        else
        {

          /* slide up and reach a node, pick one side at random and recurse */
          *slidedist -= *timepoint - sisuptime;
          *timepoint = sisuptime;
          if (bitran () /*uniform() < 0.5 */ )
          {
            *sis = gtree[*sis].up[0];
          }
          else
          {
            *sis = gtree[*sis].up[1];
          }

          /* reset slidedist to negative, so slidingedge continues up the gtree in next call to slider */
          *slidedist = -*slidedist;
          assert (*slidedist < 0);
          slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
          return;
        }
      }
    }
  }
  else
  {
    /* go down */
    if (gtree[*sis].down == -1 || *timepoint + *slidedist < gtree[*sis].time)
    {

      /* if sis is the root, or distance is less than to next down node, just slide down */
      *timepoint += *slidedist;
      if (*timepoint >= TIMEMAX)
        *timepoint = TIMEMAX;
      *slidedist = 0;
      return;
    }
    else
    {

      /* a down node is reached */
      *slidedist -= gtree[*sis].time - *timepoint;
      *timepoint = gtree[*sis].time;
      if (bitran ())
      {
        /* begin to slide down the down node */
        *sis = gtree[*sis].down;
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
      else
      {
        /* begin to slide up the sis  */
        if (gtree[gtree[*sis].down].up[0] == *sis)
        {
          *sis = gtree[gtree[*sis].down].up[1];
        }
        else
        {
          *sis = gtree[gtree[*sis].down].up[0];
        }
        *slidedist = -*slidedist;
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
}                               /* slider_nomigration */

void
slider (int ci, int li, int slidingedge, int *sis, double *timepoint,
        double *slidedist)
/* this is not ready for case of multiple populations and zero migration */
/* timepoint points at C[ci]->G[li].gtree[*slidingedge].time and is the current position of the sliding point, 
slidedist is the distance it must move 
do not restructure the gtree. just figure out when and on which branch timepoint ends up on, 
this will be the new sisterbranch
use recursion */
{
  double uplimit;
  struct edge *gtree = C[ci]->G[li].gtree;
  if (*slidedist < 0)
  {
    /* go up */
    *slidedist = -*slidedist;
    if (slidingedge < L[li].numgenes)
      uplimit = 0;
    else
      uplimit = gtree[gtree[slidingedge].up[0]].time;
    assert (*timepoint >= uplimit);

    /* if uplimit >= sis up time  - slidingedge cannot reach a node */
    if (gtree[*sis].up[0] == -1 || uplimit >= gtree[gtree[*sis].up[0]].time)
    {
      if (*slidedist < *timepoint - uplimit)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > uplimit);
        return;
      }
      else
      {
        /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
        *slidedist -= *timepoint - uplimit;
        *timepoint = uplimit;
        assert (*slidedist > 0);
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
    else
    {
      /* uplimit is less than sis up time, and thus slidingedge can reach a node */
      if (*slidedist < *timepoint - gtree[gtree[*sis].up[0]].time)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > gtree[gtree[*sis].up[0]].time);
        return;
      }
      else
      {
        /* slide up and reach a node, pick one side at random and recurse */
        *slidedist -= *timepoint - gtree[gtree[*sis].up[0]].time;
        *timepoint = gtree[gtree[*sis].up[0]].time;
        if (bitran () /*uniform() < 0.5 */ )
        {
          *sis = gtree[*sis].up[0];
        }
        else
        {
          *sis = gtree[*sis].up[1];
        }

        /* reset slidedist to negative, so slidingedge continues up the gtree in next call to slider */
        *slidedist = -*slidedist;
        assert (*slidedist < 0);
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
  else
  {
    /* go down */
    if (gtree[*sis].down == -1 || *timepoint + *slidedist < gtree[*sis].time)
    {
      /* if sis is the root, or distance is less than to next down node, just slide down */
      *timepoint += *slidedist;
      if (*timepoint >= TIMEMAX)
        *timepoint = TIMEMAX;
      *slidedist = 0;
      return;
    }
    else
    {
      /* a down node is reached */
      *slidedist -= gtree[*sis].time - *timepoint;
      *timepoint = gtree[*sis].time;
      if (bitran ())
      {
        /* begin to slide down the down node */
        *sis = gtree[*sis].down;
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
      else
      {
        /* begin to slide up the sis  */
        if (gtree[gtree[*sis].down].up[0] == *sis)
        {
          *sis = gtree[gtree[*sis].down].up[1];
        }
        else
        {
          *sis = gtree[gtree[*sis].down].up[0];
        }
        *slidedist = -*slidedist;
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
}                               /* slider */

void
joinsisdown (int ci, int li, int sis, int *tmrcachange)
{

  /* extend sis, and XFREE up the down edge */
  int i, j, ai, downdown, down;
  double uptime;
  struct edge *gtree = C[ci]->G[li].gtree;
  down = gtree[sis].down;
  i = 0;
  while (gtree[sis].mig[i].mt > -0.5)
    i++;
  j = -1;

  do
  {
    j++;
    checkmig (i + 1, &(gtree[sis].mig), &(gtree[sis].cmm));
    gtree[sis].mig[i] = gtree[down].mig[j];
    i++;
  } while (gtree[down].mig[j].mt > -0.5);
  gtree[sis].time = gtree[down].time;

  /* set the up to which sis now connects */
  gtree[sis].down = gtree[down].down;
  downdown = gtree[sis].down;
  if (downdown != -1)
  {
    rootmove = 0;
    if (gtree[downdown].up[0] == down)
      gtree[downdown].up[0] = sis;
    else
      gtree[downdown].up[1] = sis;
    mrootdrop = 0;
    lmrootdrop = 0;
  }
  else
  {
    rootmove = 1;
    *tmrcachange += 1;
/* figure out total time and number of migrants being dropped */
    i = -1;
    do
    {
      i++;
    } while (gtree[sis].mig[i].mt > -0.5);

    /* mrootdrop and lmrootdrop do not seem to do anything. We may want to
     * delete two variables? */
    mrootdrop = i; 
    if (sis < L[li].numgenes)
    {
      uptime = 0;
    }
    else
    {
      uptime = gtree[gtree[sis].up[0]].time;
    }

    if (uptime < C[ci]->tvals[lastperiodnumber - 1])
    {
      if (C[ci]->G[li].roottime < C[ci]->tvals[lastperiodnumber - 1])
        lmrootdrop = C[ci]->G[li].roottime - uptime;
      else
        lmrootdrop = C[ci]->tvals[lastperiodnumber - 1] - uptime;
    }
    else
    {
      lmrootdrop = 0;
    }

    C[ci]->G[li].root = sis;
    gtree[sis].down = -1;
    gtree[sis].time = TIMEMAX;
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
        gtree[sis].dlikeA[ai] = 0;
    gtree[sis].mig[0].mt = -1;
    C[ci]->G[li].roottime = uptime;
  }
}                               /* joinsisdown */

void
splitsisdown (int ci, int li, int slidingedge, int down, int newsis)
{

  /* split newsis into two parts, and make a new down edge out of the lower part */
  int i, j, downdown, nowpop;
  double curt;
  struct edge *gtree = C[ci]->G[li].gtree;
  curt = gtree[slidingedge].time;
  gtree[down].time = gtree[newsis].time;
  gtree[newsis].time = curt;

  /* set the up  of the edge to which down now connects, depends on whether newsis is the root */
  downdown = gtree[newsis].down;
  if (downdown != -1)
  {
    if (gtree[downdown].up[0] == newsis)
      gtree[downdown].up[0] = down;
    else
      gtree[downdown].up[1] = down;
  }
  else
  {
    /* newsis is the current root so the root must move down */
    C[ci]->G[li].root = down;
    C[ci]->G[li].roottime = curt;
    rootmove = 1;
    if (C[ci]->G[li].roottime > TIMEMAX)
      IM_err(IMERR_ROOTTIMEMAXFAIL, "roottime greater than TIMEMAX, chain: %d,locus: %d, roottime %lf, TIMEMAX %lf",ci,li,C[ci]->G[li].roottime,TIMEMAX);
    gtree[down].mig[0].mt = -1;
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      for (i = (L[li].model == JOINT_IS_SW); i < L[li].nlinked; i++)
        gtree[down].dlikeA[i] = 0;

    /*if (L[li].model== STEPWISE)
       for (i=0;i< L[li].nlinked;i++)
       gtree[down].dlikeA[i] = 0;
       if (L[li].model== JOINT_IS_SW)
       for (i=1;i< L[li].nlinked;i++)
       gtree[down].dlikeA[i] = 0; */
  }
  gtree[down].down = downdown;

  /* divide the migration events along newsis into upper part for newsis and lower part for down */
  /* this might have bugs setting the population of gtree[down] */
  i = 0;
  while (gtree[newsis].mig[i].mt > -0.5 && gtree[newsis].mig[i].mt < curt)
    i++;
  if (i > 0)
    nowpop = gtree[newsis].mig[i - 1].mp;
  else
    nowpop = gtree[newsis].pop;
  j = findperiod (ci, curt);
  while (C[ci]->poptree[nowpop].e <= j && C[ci]->poptree[nowpop].e != -1)
    nowpop = C[ci]->poptree[nowpop].down;
  gtree[down].pop = nowpop;
  j = 0;
  if (downdown != -1)
  {
    while (gtree[newsis].mig[j + i].mt > -0.5)
    {
      checkmig (j, &(gtree[down].mig), &(gtree[down].cmm));
      gtree[down].mig[j] = gtree[newsis].mig[j + i];
      assert (nowpop != gtree[down].mig[j].mp);
      nowpop = gtree[down].mig[j].mp;
      j++;
    }
  }
  gtree[newsis].mig[i].mt = -1;
  gtree[down].mig[j].mt = -1;
  gtree[newsis].down = gtree[slidingedge].down = down;
  gtree[down].up[0] = newsis;
  gtree[down].up[1] = slidingedge;
  return;
}                               /* splitsisdown */


/* called by addmigration(),  does most of the migration either by adding events or by calling mwork() 
    mwork is called for the simpler cases, 
    this function handles the case where both edges need migration and the state of the population at the bottom of the edges needs to be 
    simulated  */
void
getm (int ci, double mrate, struct edgemiginfo *edgem,
      struct edgemiginfo *sisem)
{
  int i, ii, j;
  int mp, mpall1, mpall2, pop, pop1, pop2, lastpop1;
  double tempt;
  double mlist[2 * ABSMIGMAX], uptime1, uptime2;        // mlist very wasteful of space - could use dynamic memory and checkmig() 
  if (sisem->mtall <= 0)
  {
    edgem->mpall = mwork (ci, edgem, edgem->e, mrate);
  }
  else
  {
    assert (edgem->e == sisem->e);
    if (edgem->mtall <= 0)
    {
      sisem->mpall = mwork (ci, sisem, sisem->e, mrate);
    }
    else
    {
      if (edgem->e == lastperiodnumber)
      {
        assert (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0);
        edgem->mpall = mwork (ci, edgem, edgem->e, mrate);
        sisem->mpall = mwork (ci, sisem, edgem->e, mrate);

        // set the fpop values of both edges to be the root population 
        while (C[ci]->poptree[edgem->temppop].e != -1)
          edgem->temppop = C[ci]->poptree[edgem->temppop].down;
        edgem->fpop = sisem->fpop = edgem->temppop;
      }
      else                      // for the remaining cases, there is less generality so do the work in this function rather than in mwork() and simmpath()
      {
        /* first do periods < lastperiodnumber-1 */
        if (edgem->e > 0)
          edgem->mpall = mpall1 = mwork (ci, edgem, edgem->e - 1, mrate);
        else
          edgem->mpall = mpall1 = 0;

        if (sisem->e > 0)
          sisem->mpall = mpall2 = mwork (ci, sisem, sisem->e - 1, mrate);
        else
          sisem->mpall = mpall2 = 0;

        /* now do the last period,  when the edges join */
        pop1 = edgem->temppop;
        pop2 = sisem->temppop;
        if (edgem->e > 0)
          uptime1 = DMAX (edgem->upt, C[ci]->tvals[edgem->e - 1]);
        else
          uptime1 = edgem->upt;
        if (sisem->e > 0)
          uptime2 = DMAX (sisem->upt, C[ci]->tvals[sisem->e - 1]);
        else
          uptime2 = sisem->upt;
        tempt = edgem->mtimeavail[edgem->e] + sisem->mtimeavail[sisem->e];
        /* if assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1  do not enter here */
        if (edgem->e == lastperiodnumber - 1
            && (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
                || npops == 2))
          // the last period of the branches is when there are only 2 pops
        {
          mp = poisson (mrate * tempt, pop1 != pop2);   /* odd if different pops,  even if same */

          //oddcheck += pop1 != pop2;
          //check++;
          for (i = 0; i < mp; i++)
            mlist[i] = uniform () * (tempt);
          if (mp > 1)
            sort (mlist - 1, (unsigned long) mp);

          /* now split these migration events onto the two branches, based on where they fall in the total time interval */
          i = 0;
          ii = mpall1;
          pop = pop1;
          while (mlist[i] + uptime1 < edgem->dnt && i < mp)

          {
            edgem->mig[ii].mt = mlist[i] + uptime1;
            edgem->mig[ii].mp = pop =
              picktopop2 (pop, C[ci]->plist[edgem->e], 2, pop);
            i++;
            ii++;
          }
          edgem->mp[edgem->e] = i;
          edgem->mpall += i;
          if (ii > 0)           // at least one migration 
          {
            edgem->mig[ii].mt = -1;
            edgem->fpop = edgem->mig[ii - 1].mp;
          }
          else
          {
            edgem->fpop = pop;
          }
          while (C[ci]->poptree[edgem->fpop].e <=
                 (lastperiodnumber - 1)
                 && edgem->dnt >
                 C[ci]->tvals[C[ci]->poptree[edgem->fpop].e - 1])
          {
            assert (C[ci]->poptree[edgem->fpop].e != -1);
            edgem->fpop = C[ci]->poptree[edgem->fpop].down;
          }
          sisem->fpop = edgem->fpop;
          ii = mpall2;
          pop = pop2;
          for (j = i; j < mp; j++, ii++)
          {
            sisem->mig[ii].mt = mlist[j] + uptime2 - (edgem->dnt - uptime1);
            sisem->mig[ii].mp = pop =
              picktopop2 (pop, C[ci]->plist[edgem->e], 2, pop);
          }
          sisem->mp[sisem->e] = j - i;
          sisem->mpall += j - i;
          if (j - i > 0)
            assert (edgem->mp[edgem->e] + sisem->mp[sisem->e] == mp);
          if (ii > 0)
          {
            sisem->mig[ii].mt = -1;
            sisem->fpop = sisem->mig[ii - 1].mp;
          }
          else
          {
            sisem->fpop = pop;
          }
          while (C[ci]->poptree[sisem->fpop].e <=
                 (lastperiodnumber - 1)
                 && sisem->dnt >
                 C[ci]->tvals[C[ci]->poptree[sisem->fpop].e - 1])
          {
            assert (C[ci]->poptree[sisem->fpop].e != -1);
            sisem->fpop = C[ci]->poptree[sisem->fpop].down;
          }
          assert (sisem->fpop == edgem->fpop);
        }
        else                    // the last period of the branches is when there are more than two pops, (there are npops - edgem->e populations) 
          // code is mostly the same as before this else
        {
          if (pop1 != pop2)
            mp = poisson (mrate * (npops - (edgem->e + 1)) * tempt, 2); /* must pick a number > 0 */
          else
            mp = poisson (mrate * (npops - (edgem->e + 1)) * tempt, 3); /* must pick a number != 1 */
          for (i = 0; i < mp; i++)
            mlist[i] = uniform () * (tempt);
          if (mp > 1)
            //hpsortreg(mlist-1, mp);
            sort (mlist - 1, (unsigned long) mp);

          /* now split these migration events onto the two branches, based on where they fall in the total time interval */
          i = 0;
          ii = mpall1;
          pop = pop1;

          /* i is the number of migrations on the edgem path there mp - i is the number on the sisem path */
          while (mlist[i] + uptime1 < edgem->dnt && i < mp)
          {
            edgem->mig[ii].mt = mlist[i] + uptime1;
            if (i == mp - 2)
            {
              pop = picktopop2 (pop, C[ci]->plist[edgem->e], npops - edgem->e, pop2);
            }
            else
            {
              if (i == mp - 1)
                pop = pop2;

              else
                pop = picktopop (pop, C[ci]->plist[edgem->e], npops - edgem->e);
            }
            edgem->mig[ii].mp = pop;
            i++;
            ii++;
          }
          lastpop1 = pop;
          edgem->mp[edgem->e] = i;
          edgem->mpall += i;
          if (ii > 0)
          {
            edgem->mig[ii].mt = -1;
            edgem->fpop = edgem->mig[ii - 1].mp;
          }
          else
          {
            edgem->fpop = pop;
          }
          
          while (C[ci]->poptree[edgem->fpop].e <=
                 (lastperiodnumber - 1)
                 && edgem->dnt >
                 C[ci]->tvals[C[ci]->poptree[edgem->fpop].e - 1])

          {
            assert (C[ci]->poptree[edgem->fpop].e != -1);
            edgem->fpop = C[ci]->poptree[edgem->fpop].down;
          }
          sisem->fpop = edgem->fpop;
          ii = mpall2;
          pop = pop2;
          for (j = i; j < mp; j++, ii++)
          {
            sisem->mig[ii].mt = mlist[j] + uptime2 - (edgem->dnt - uptime1);
            if (j == mp - 2)
            {
              pop = picktopop2 (pop, C[ci]->plist[sisem->e], npops - sisem->e, lastpop1);

            }
            else
            {
              if (j == mp - 1)
              {
                assert (pop != lastpop1);
                pop = lastpop1;
              }
              else
              {
                pop = picktopop (pop, C[ci]->plist[sisem->e], npops - sisem->e);
              }
            }
            sisem->mig[ii].mp = pop;
          }
          sisem->mp[sisem->e] = j - i;
          sisem->mpall += j - i;
          assert (edgem->mp[edgem->e] + sisem->mp[sisem->e] == mp);
          if (ii > 0)
          {
            sisem->mig[ii].mt = -1;
          }
          else
          {
            sisem->fpop = pop;
          }
          while (C[ci]->poptree[sisem->fpop].e <=
                 (lastperiodnumber - 1)
                 && sisem->dnt >
                 C[ci]->tvals[C[ci]->poptree[sisem->fpop].e - 1])

          {
            assert (C[ci]->poptree[sisem->fpop].e != -1);
            sisem->fpop = C[ci]->poptree[sisem->fpop].down;
          }
          assert (sisem->fpop == edgem->fpop);
        }
        i = findperiod (ci, edgem->dnt);
        while (C[ci]->poptree[edgem->fpop].e <= i)
          edgem->fpop = C[ci]->poptree[edgem->fpop].down;
      }
    }
  }
}                               /* getm */


/* add migration to edge, and to its sister if edge connects to the root, return the log of the hastings ratio of 
    update probabilities */
  /* add migration events to edge that just slid.  Also if it slid down the root node, and moved the root, then migration 
     events may need to be added to the sister branch as well */
  /* oldtlength is the length of the gtree more recent than the basal population split */
  /* oldmigcount is just the total number of migration events on the gtree */

/* jh changed 7/7/08 
added newmigcount and newtlength
*/

/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre tree - e.g. help to avoid promoting
correlations and improve mixing  */

double
addmigration (int ci, int li, int oldmigcount, double oldtlength,
              int *newmigcount, double *newtlength)
{
  int newsis, edge;
  double weight;
  double mproposenum, mproposedenom, temp;
  double mparamf, mparamb;
  struct edge *gtree = C[ci]->G[li].gtree;
  double mtime;
  int mcount;

  assert (C[ci]->G[li].mignum >= 0 && C[ci]->G[li].tlength > 0);

  /* determine current migration rate to use for update */
  //mparamf = calcmrate (C[ci]->G[li].mignum, C[ci]->G[li].tlength);

  mparamf = 0;
  if (oldsismig.edgeid >= 0)
  {
    if (oldsismig.mtall + oldedgemig.mtall > 0)
      mparamf = calcmrate (oldsismig.mpall + oldedgemig.mpall, oldsismig.mtall + oldedgemig.mtall);
  }
  else
  {
    if (oldedgemig.mtall > 0)
      mparamf = calcmrate (oldedgemig.mpall, oldedgemig.mtall);

  }

  //mparamf = calcmrate (C[ci]->G[li].mignum, C[ci]->G[li].tlength);


//printf(" mparamf %lf \n",mparamf);


  /* store information on edge, and sister edge if needed */
/* REMOVED!
  memset (&newedgemig, 0, sizeof (struct edgemiginfo));
*/
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
/* REMOVED!
      memset (&newsismig, 0, sizeof (struct edgemiginfo));
*/
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
  assert (mparamf > 0 || oldedgemig.mtall + oldsismig.mtall == 0);
  getm (ci, mparamf, &newedgemig, &newsismig);
  if (gtree[newedgemig.edgeid].down == C[ci]->G[li].root)
    gtree[C[ci]->G[li].root].pop = newedgemig.fpop;
  if (newedgemig.mtall == 0 && newsismig.mtall == 0)    // both new edges occur after the last split time
    temp = 0;
  else
    temp = getmprob (ci, mparamf, &newedgemig, &newsismig);
  mproposedenom = temp;
  assert (temp > -1e200 && temp < 1e200);

  /* calculate probability of reverse update    */
  /* 7/11/08  - fixed a nasty bug here, 
     was causing wrong values for mtime and mcount, which in turn was causing wrong values
     for reverse update in calculation of hastings term.
     not sure just how it might have shaped the what it was doing to results */
  mtime = oldtlength - oldedgemig.mtall - oldsismig.mtall + newedgemig.mtall;
  mcount =
    oldmigcount - oldedgemig.mpall - oldsismig.mpall + newedgemig.mpall;
  if (newsismig.edgeid >= 0)
  {
    mtime += newsismig.mtall;
    mcount += newsismig.mpall;
  }
  /* old version 
     mtime = oldtlength - oldedgemig.mtall + newedgemig.mtall;
     mcount = oldmigcount - oldedgemig.mpall + newedgemig.mpall;
     if (newsismig.edgeid >= 0)

     {
     mtime += newsismig.mtall - oldsismig.mtall;
     mcount += newsismig.mpall - oldsismig.mpall;
     }
   */
/* find the migation rate for the backward update */
  //mparamb = calcmrate (mcount, mtime); // is it correct to use the inverse multiplier for the reverse update ????
  

  mparamb = 0;
  if (newsismig.edgeid >= 0)
  {
    if (newsismig.mtall + newedgemig.mtall > 0)
      mparamb = calcmrate (newsismig.mpall + newedgemig.mpall, newsismig.mtall + newedgemig.mtall);
  }
  else
  {
    if (newedgemig.mtall > 0)
      mparamb = calcmrate (newedgemig.mpall, newedgemig.mtall);
  }

  //mparamb = calcmrate (mcount, mtime);

  assert (mparamb > 0 || newedgemig.mtall + newsismig.mtall == 0);
  if (oldedgemig.mtall == 0 && oldsismig.mtall == 0)    // both new edges occur after the last split time
    temp = 0;
  else
    temp = getmprob (ci, mparamb, &oldedgemig, &oldsismig);
  mproposenum = temp;
  assert (temp > -1e200 && temp < 1e200);
  weight = mproposenum - mproposedenom;
  *newmigcount = mcount;
  *newtlength = mtime;
  return weight;
}                               /* addmigration */

/********GLOBAL FUNCTIONS *******/

void
init_updategenealogy (void)
{
  init_genealogy_weights (&holdallgweight_updategenealogy);
  init_genealogy_weights (&holdgweight_updategenealogy);
  init_probcalc (&holdallpcalc_updategenealogy);
}                               // init_updategenealogy


void
free_updategenealogy (void)
{
  free_genealogy_weights (&holdallgweight_updategenealogy);
  free_genealogy_weights (&holdgweight_updategenealogy);
  free_probcalc (&holdallpcalc_updategenealogy);
}                               //free_updategenealogy

#define SLIDESTDVMAX 20

/* steps in picking a new genealogy 
- pick an edge, the bottom of which will slide
- save all info for that edge, sis and the down edge that will be freed up
- join the sis and down edges, freeing up an edge
	if down was the root, then sis becomes the root
-set aside the number for the down edge - this is the freed up edge and will get used again later
-do the sliding and pick a new sis and a location - but do not change the gtree. 
-split the newsis edge into sis and down edges, and divide the migration events accordingly
-connect the original edge to the new spot
-add migration events to this edge
- calculate the probability of the update, in forwrad and reverse directions
- calculate the total Metropolis Hastings terms for genealogy update,  accept or reject 
*/
/* for gtreeprint calls,  use callingsource = 0 */


/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre tree - e.g. help to avoid promoting
correlations and improve mixing  */
int
updategenealogy (int ci, int li, int *topolchange, int *tmrcachange)
{
  int ai, ui, mpart, i;
  int updateAcount;
  int edge, oldsis, newsis, freededge, accp;
  double newpdg, newplg, newpdg_a[MAXLINKED];
  double migweight, metropolishastingsterm, U;
  double tpw;
  double Aterm[MAXLINKED], Atermsum;
  /* static struct genealogy_weights holdgweight_updategenealogy, holdallgweight_updategenealogy; */
  /* static struct probcalc holdallpcalc_updategenealogy; */
  double slidedist;
  double tlengthpart;
  double slideweight, holdslidedist, slidestdv;
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  int rejectIS;
  double like;
  double holdt[MAXPERIODS];
/* change this HPDBG section 4/2/08 */
#ifdef HPDBG
// Set the debug-heap flag so that freed blocks are kept on the
// linked list, to catch any inadvertent use of freed memory
//SET_CRT_DEBUG_FIELD( _CRTDBG_DELAY_FREE_MEM_DF );
#endif /*  */

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogyloc (ci, li);
  }

// initialize and make copies structures that hold quantities for calculating prob of genealogy
  copy_treeinfo (&holdgweight_updategenealogy, &G->gweight);
  copy_treeinfo (&holdallgweight_updategenealogy, &C[ci]->allgweight);
//setzero_genealogy_weights (&G->gweight);
//treeweight (ci, li);
// store summary stats of the genealogy
  storetreestats (ci, li, 0);
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];

/* tlengthpart is the part of the tree on which migration can happen,  mpart is the number of migration event - these 
	are used when updating migration */
  tlengthpart = G->tlength;
  mpart = G->mignum;

// Atermsum only used for Stepwise mutation model
  Atermsum = 0;

// indicitor variables of types of accepted updates to genealogy
  *tmrcachange = 0;
  *topolchange = 0;

/* pick an edge, identify freedup edge (the down edge) and the sister edge */
  do
  {
    edge = randposint (L[li].numlines);
  } while (gtree[edge].down == -1);
  freededge = gtree[edge].down;
  if ((oldsis = gtree[freededge].up[0]) == edge)
    oldsis = gtree[freededge].up[1];

  /* copy information on the edge,  and if it connects to the root, then the sister edge as well */
  if (gtree[edge].down == G->root)
    fillmiginfo (ci, li, gtree, edge, oldsis);
  else
    fillmiginfo (ci, li, gtree, edge, -1);

/* store information on the genealogy before changing it */
  storeoldedges (ci, li, edge, oldsis, freededge);

// remove any migrations  from the slidingedge 
  gtree[edge].mig[0].mt = -1;

/* slide edge, pick a distance and slide it  */
  slidestdv = DMIN (SLIDESTDVMAX, G->roottime / 3);
  holdslidedist = slidedist = normdev (0.0, slidestdv);

// joint the sister and the down branches at the point where edge used to connect, this frees up the down branch 
  joinsisdown (ci, li, oldsis, tmrcachange);

// do the slide and identify the new sister branch and where new connection point for the edge is 
  newsis = oldsis;
  if (modeloptions[NOMIGRATION])
    slider_nomigration (ci, li, edge, &newsis, &(gtree[edge].time),
                        &slidedist);
  else
    slider (ci, li, edge, &newsis, &(gtree[edge].time), &slidedist);
  *topolchange += (oldsis != newsis);

// now separate the new sister branch into a shorter sis branch and a down branch 
  splitsisdown (ci, li, edge, freededge, newsis);
  if (rootmove)
  {
    slideweight = -log (normprob (0.0, slidestdv, holdslidedist));
    slidestdv = DMIN (SLIDESTDVMAX, G->roottime / 3);
    slideweight += log (normprob (0.0, slidestdv, holdslidedist));
  }
  else
  {
    slideweight = 0;
  }

// add migration events 
  if (modeloptions[NOMIGRATION])
  {
    migweight = 0;
  }
  else
  {
    migweight = addmigration (ci, li, mpart, tlengthpart, &mpart, &tlengthpart);
  }

  // copy the migration info in newedgemig and newsismig  to the genealogy
  copynewmig_to_gtree (ci, li);

  if (modeloptions[NOMIGRATION] == 0)
  {
    if (assignmentoptions[POPULATIONTREEWCH] == 1)
    {
      /* migweight += addmigrationf (ci, li); */
    }
  }

// determine all the weights needed for calculating the probability of the genealogy
  setzero_genealogy_weights (&G->gweight);
  treeweight (ci, li);
  assert (G->mignum == mpart);  //G->tlength should also be very close to tlengthpart
//    if (!modeloptions[NOMIGRATION])
  //      assert(fabs(G->tlength - tlengthpart) < 1e-8);
  sum_subtract_treeinfo (&C[ci]->allgweight, &G->gweight,
                         &holdgweight_updategenealogy);

/* calculate P(D|G)  for new genealogy */
  rejectIS = 0;                 /* in case P(D|G) for IS model is zero */
  newpdg = 0;
  newplg = 0;
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
        likelihoodHKY (ci, li, G->uvals[0], G->kappaval, edge,
                       freededge, oldsis, newsis);
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
        newpdg_a[ai] =
          G->pdg_a[ai] + finishSWupdateA (ci, li, ai, edge, freededge,
                                          oldsis, newsis,
                                          G->uvals[ai], &Aterm[ai]);

        newpdg += newpdg_a[ai];
        Atermsum += Aterm[ai];
      }

//            checklikelihoodSW(ci, li,G->u[ai].mcinf.val);  
      break;
    }
  case JOINT_IS_SW:
    newpdg = newpdg_a[0] = likelihoodIS (ci, li, G->uvals[0]);
    rejectIS = (newpdg == REJECTINFINITESITESCONSTANT);
    for (ai = 1; ai < L[li].nlinked; ai++)
    {
      newpdg_a[ai] =
        G->pdg_a[ai] + finishSWupdateA (ci, li, ai, edge, freededge,
                                        oldsis, newsis,
                                        G->uvals[ai], &Aterm[ai]);
      newpdg += newpdg_a[ai];
      Atermsum += Aterm[ai];
    }

    //checklikelihoodSW(ci, li,Q[ci]->us[li]);  
    break;
  }

  // LEVINE LIKELIHOOD: 2008-07-11
  if (calcoptions[CALCLEVINELIKELIHOOD])
  {
    if (*topolchange)
      newplg += likelihoodDG (ci, li);
    else 
      newplg += G->plg;
  }
  accp = 0;

/* final weight calculation */
/* tpw is the ratio of new and old prior probability of the genealogies.  It is actually the ratio of the total across all loci,  but
since only genealogy li is being changed at the present time,  the ratio works out to just be the ratio for genealogy li */
  copy_probcalc (&holdallpcalc_updategenealogy, &C[ci]->allpcalc);

  /* Find all internal node sequences and mutations of a full genealogy. */
  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    if (L[li].model == INFINITESITES)
    {
      accp = IMA_genealogy_findIntSeq (ci, li);
      if (accp == 0)
      {
        rejectIS = 1;
      }
      accp = 0;
    }
  }

  if (rejectIS == 0)
  {
    // the metropolis term includes p(D|G) and p(G),  
    tpw = -C[ci]->allpcalc.probg;
    integrate_tree_prob (ci, &C[ci]->allgweight,
                         &holdallgweight_updategenealogy, &C[ci]->allpcalc,
                         &holdallpcalc_updategenealogy, &holdt[0]);
    tpw += C[ci]->allpcalc.probg;
    metropolishastingsterm = tpw + gbeta * (newpdg - G->pdg + newplg - G->plg);
    U = uniform ();
    //    assert (beta[ci] * metropolishastingsterm + migweight + Atermsum > -1e200
    //             && beta[ci] * metropolishastingsterm + migweight + Atermsum < 1e200);
//printf("step %d weight %lf  slideweight %lf  migweight %lf ",step,weight, slideweight, migweight);
    metropolishastingsterm = exp (beta[ci] * metropolishastingsterm + migweight + slideweight + Atermsum);
//printf(" total weight %lf\n",weight);
    //   assert (metropolishastingsterm >= 0);
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
    // put the edges back 
    restoreedges (ci, li, edge, oldsis, freededge, newsis);

    // copy summary stats back
    storetreestats (ci, li, 1);

    // reset HKY terms
    if (L[li].model == HKY)
      restorescalefactors (ci, li);
    // copy back all the weights and results associated with calculating the probability of the genealogy 
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_updategenealogy);
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_updategenealogy);
    copy_treeinfo (&G->gweight, &holdgweight_updategenealogy);
    *topolchange = 0;
    *tmrcachange = 0;
  }

/* do updates at nodes for stepwise loci, regardless of whether slide update was accepted.  This could go somewhere else  */
  if (!calcoptions[DONTCALCLIKELIHOODMUTATION])
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    {
      C[ci]->allpcalc.pdg -= G->pdg;
      for (ui = (L[li].model == JOINT_IS_SW); ui < L[li].nlinked; ui++) // ui starts at 1 if JOINT_IS_SW otherwise 0
      {
        updateAcount = 0;
        G->pdg -= G->pdg_a[ui];
        G->pdg_a[ui] = updateA (ci, li, ui, G->uvals[ui], &updateAcount);
        G->pdg += G->pdg_a[ui];
      }
      C[ci]->allpcalc.pdg += G->pdg;
    }
//assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogyloc (ci, li);
  }
  return accp;
}                               /* update_gtree */

#undef SLIDESTDVMAX
