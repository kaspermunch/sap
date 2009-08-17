/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"

/****** LOCAL STUFF **********/

#define MAX_FOR_MULTI_T_ARRAY  3

static int **t2array;
static int ***t3array;
static double tarraybinvals[MAX_FOR_MULTI_T_ARRAY][NUMTARRAYBINS];
static int tvalarraypos[MAX_FOR_MULTI_T_ARRAY] = { 0 };
static int highpos[MAX_FOR_MULTI_T_ARRAY];
static int numt;
static int init = 0;

/* set up multidimensional arrays for t values */
/* mode 0  means record the array of t values */
/* mode 1 means return the estimated highest t value */

/******* GLOBAL FUNCTIONS ***********/

void
setup_multi_t_arrays ()
{

  int i, j;
  double tscaleadjust;

  if (modeloptions[SPLITTINGRATEPARAMETER])
    tscaleadjust = TIMERECORDPRIORFRAC;
  else
    tscaleadjust = 1;

  if (init == 0)
  {
    numt = npops - 1;
    for (i = 0; i < numt; i++)
      for (j = 0; j < NUMTARRAYBINS; j++)
        tarraybinvals[i][j] = T[i].pr.min + ((j + 0.5) * (T[i].pr.max * tscaleadjust - T[i].pr.min)) / NUMTARRAYBINS;

    if (numt == 2)
    {
      if ((t2array = malloc (NUMTARRAYBINS * sizeof (*t2array))) == NULL)
        IM_err (IMERR_MEM, "  t2array problem in setup_multi_t_arrays()");
      for (i = 0; i < NUMTARRAYBINS; i++)
        if ((t2array[i] = calloc ((size_t) (NUMTARRAYBINS - i), sizeof (*t2array[i]))) == NULL)
          IM_err (IMERR_MEM, "  t2array problem in setup_multi_t_arrays()");
    }

    if (numt == 3)
    {
      if ((t3array = malloc (NUMTARRAYBINS * sizeof (*t3array))) == NULL)
        IM_err (IMERR_MEM, "  t3array problem in setup_multi_t_arrays()");
      for (i = 0; i < NUMTARRAYBINS; i++)
      {
        if ((t3array[i] = malloc ((NUMTARRAYBINS - i) * sizeof (*t3array[i]))) == NULL)
          IM_err (IMERR_MEM, "  t2array problem in setup_multi_t_arrays()");
        for (j = 0; j < (NUMTARRAYBINS - i); j++)
          if ((t3array[i][j] =
               calloc ((size_t) (NUMTARRAYBINS - j), sizeof (*t3array[i][j]))) == NULL)
            IM_err (IMERR_MEM, "  t2array problem in setup_multi_t_arrays()");
      }
    }
    for (i = 0; i < numt; i++)
      highpos[i] = 0;
  }
  init++;
  for (i = 0; i < numt; i++)
  {
    if (C[0]->tvals[i] > T[i].pr.max * tscaleadjust)
      return;
    if (i == 0)
      tvalarraypos[i] = (int) (NUMTARRAYBINS * ((C[0]->tvals[i] - T[i].pr.min) / (T[i].pr.max * tscaleadjust - T[i].pr.min)));
    else
      tvalarraypos[i] = (int) (NUMTARRAYBINS * ((C[0]->tvals[i] - T[i].pr.min) / (T[i].pr.max * tscaleadjust - T[i].pr.min))) - tvalarraypos[i - 1];
  }
  if (numt == 2)
  {
    t2array[tvalarraypos[0]][tvalarraypos[1]]++;
    if (t2array[tvalarraypos[0]][tvalarraypos[1]] >
        t2array[highpos[0]][highpos[1]])
    {
      highpos[0] = tvalarraypos[0];
      highpos[1] = tvalarraypos[1];
    }
  }
  if (numt == 3)
  {
    t3array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]]++;
    if (t3array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]] >
        t3array[highpos[0]][highpos[1]][highpos[2]])
    {
      highpos[0] = tvalarraypos[0];
      highpos[1] = tvalarraypos[1];
      highpos[2] = tvalarraypos[2];
    }
  }
}                               //setup_multi_t_arrays()

void
free_multi_t_arrays ()
{
  int i, j;

  if (numt == 2)
  {
    free2D ((void **) t2array, NUMTARRAYBINS);
  }

  if (numt == 3)
  {
    for (i = 0; i < NUMTARRAYBINS; i++)
    {
      for (j = 0; j < (NUMTARRAYBINS - i); j++)
        XFREE (t3array[i][j]);
      XFREE (t3array[i]);
    }
    XFREE (t3array);
  }
  if (numt < 2 || numt > 3)
  {
    // should not get here
    return;
  }
}                               //free_multi_t_arrays()


void
return_joint_t (double tvals[])
{
  int i;
  for (i = 0; i < numt; i++)
  {
    if (i == 0)
      tvals[i] = tarraybinvals[i][highpos[i]];
    else
      tvals[i] = tarraybinvals[i][highpos[i] + highpos[i - 1]];
  }
}

// estimate joint posterior probability of a t value 
//  not sure what this would be used for 
double
joint_t_prob (double *tvals)
{
  int i;
  double tscaleadjust;

  if (modeloptions[SPLITTINGRATEPARAMETER])
    tscaleadjust = TIMERECORDPRIORFRAC;
  else
    tscaleadjust = 1;
  for (i = 0; i < numt; i++)
  {
    if (i == 0)
      tvalarraypos[i] = (int) (NUMTARRAYBINS * ((tvals[i] - T[i].pr.min) / (T[i].pr.max * tscaleadjust - T[i].pr.min)));
    else
      tvalarraypos[i] = (int) (NUMTARRAYBINS * ((tvals[i] - T[i].pr.min) / (T[i].pr.max * tscaleadjust - T[i].pr.min))) - tvalarraypos[i - 1];
  }
  if (numt == 2)
  {
    return (double) t2array[tvalarraypos[0]][tvalarraypos[1]] / (double) init;
  }
  if (numt == 3)
  {
    return (double) t3array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]] / (double) init;
  }
  return -1;                    // should not get here

}                               /* joint_t_prob */
