/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"


/* calculate the probability density of the product 2Nm  */
double
calc_popmig (int thetai, int mi, double x, int prob_or_like)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp, temp1, temp2, fc, fm, hc, qintg, mintg, mmax, qmax;
  double a,b;
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmax = imig[mi].pr.max;
  qmax = itheta[thetai].pr.max;

  for (sum = 0, ei = 0; ei < treessaved; ei++)

  {
    cc = (int) gsampinf[ei][ccp];
    fc = (double) gsampinf[ei][fcp];
    hc = (double) gsampinf[ei][hcp];
    mc = (int) gsampinf[ei][mcp];
    fm = (double) gsampinf[ei][fmp];
    qintg = (double) gsampinf[ei][qip];
    mintg = (double) gsampinf[ei][mip];
    if (fc == 0 && cc == 0 && fm > 0)
    {
      temp1 = LOG2 - (mc * log (fm)) - hc - qintg - mintg;
      temp2 = log (exp (uppergamma (mc, 2 * fm * x / qmax)) - exp (uppergamma (mc, mmax * fm)));
    }
    else
    {
      if (fm == 0 && mc == 0 && fc > 0)
      {
        temp1 = LOG2 - (cc * log (fc)) - hc - qintg - mintg;
        temp2 = log (exp (uppergamma (cc, 2 * fc / qmax)) - exp (uppergamma (cc, fc * mmax / x)));
      }
      else
      {
        if (fc == 0 && cc == 0 && mc == 0 && fm == 0)
        {
          temp1 = log (2 * log (mmax * qmax / (2 * x))) - hc - qintg - mintg;;
          temp2 = 0;
        }
        else// fc>0, cc>0, mc> 0, fm> 0
        {
          temp1 = LOG2 + (mc * log (x)) - ((cc + mc) * log (fc + fm * x)) - hc - qintg - mintg;
          //temp2 = log ( exp(uppergamma (cc + mc, 2 * (fc + fm * x) / qmax)) - exp(uppergamma (cc + mc, mmax * (fm + fc / x))) );
          a =uppergamma (cc + mc, 2 * (fc + fm * x) / qmax);
          b = uppergamma (cc + mc, mmax * (fm + fc / x));
		  if (a>b)
          {
			LogDiff(temp2,a,b);
          }
		  else
          {
            temp1 = temp2 = 0.0;
          }
        }
      }
    }
    if ((temp1 + temp2 < 700 ) && (temp1 + temp2 > -700 )) // skip things that cannot be exped
      sum += exp (temp1 + temp2);
  }
  sum /= treessaved;
  if (prob_or_like)
  {
    temp = 2 * (log (qmax) + log (mmax) - log (2 * x)) / (qmax * mmax);
    sum /= temp;
  }
  return sum;
}                               //calc_popmig

/* calculate the probability density of the product 2Nm when migration has an exponential prior */
#define OCUTOFF  10
double
calc_pop_expomig (int thetai, int mi, double x, int prob_or_like, struct extendnum *tempesum)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp, temp1, temp2, temp3,fc, fm, hc, qintg, mintg, mmean, qmax;
  int zadj, maxz = -1000000000;
  double acumm = 0;
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmean = imig[mi].pr.mean;
  qmax = itheta[thetai].pr.max;

  for (sum = 0, ei = 0; ei < treessaved; ei++)

  {
    cc = (int) gsampinf[ei][ccp];
    fc = (double) gsampinf[ei][fcp];
    hc = (double) gsampinf[ei][hcp];
    mc = (int) gsampinf[ei][mcp];
    fm = (double) gsampinf[ei][fmp];
    qintg = (double) gsampinf[ei][qip];
    mintg = (double) gsampinf[ei][mip];
    temp1 = x + fc *mmean + fm * mmean * x;
    temp2 = LOG2 - hc - log(mmean) - qintg - mintg;
    temp3 = 2 * temp1/(mmean * qmax);

    if (mc == 0 && cc == 0 )
    {
      temp3 = uppergamma(0,temp3);
      temp2 += temp3;
    }
    else
    {
      if (cc == 0) // mc > 0 
      {
        temp3 = uppergamma(mc,temp3);
        temp2 += temp3 - mc* log(fm + 1/mmean);
      }
      else
      {
        if (mc == 0) //cc > 0
        {
          temp3 = uppergamma(cc,temp3);
          temp2 += temp3 + -cc *log(fc + x*(fm + 1/mmean));; 
        }
        else
        {
          temp3 = uppergamma(mc+cc,temp3);
          temp2 += temp3 - cc*log(x) + (cc+mc)*log(mmean*x/temp1); 
        }
      }
    }
    eexp (temp2, &tempesum[ei].m, &tempesum[ei].z);
    if (tempesum[ei].z > maxz)
      maxz = tempesum[ei].z;
  }
  for (ei = 0; ei < treessaved; ei++)
  {
    zadj = tempesum[ei].z - (maxz - OCUTOFF);
    tempesum[ei].m *= pow (10.0, (double) zadj);
    acumm += tempesum[ei].m;
  }
  sum = log (acumm) + (maxz - OCUTOFF) * LOG10;
  sum -= log((double) treessaved);
  sum = exp(sum);
  if (prob_or_like)
  {
    temp = 2 * exp(uppergamma(0,2*x/(mmean*qmax)))/ (qmax * mmean);
    sum /= temp;
  }
  XFREE(tempesum);
  return sum;
}                               //calc_pop_expomig
#undef OCUTOFF  

#define OFFSCALEVAL 1

/*  do calculations for finding the peak of the marginal density for a 2Nm term
    this is similar to calc_popmig().  It is called by marginalopt_popmig()*/

double
marginpopmig (int mi, int firsttree, int lasttree, double x, int thetai)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp1, temp2, fc, fm, hc, qintg, mintg, mmax, qmax;
  double a,b,c,d;
  double max, min;

  max = itheta[thetai].pr.max*imig[mi].pr.max;
  min = 0;
  if (x < min || x > max)
    return OFFSCALEVAL;
  
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmax = imig[mi].pr.max;
  qmax = itheta[thetai].pr.max;

  for (sum = 0,ei = firsttree; ei < lasttree; ei++)

  {
    cc = (int) gsampinf[ei][ccp];
    fc = (double) gsampinf[ei][fcp];
    hc = (double) gsampinf[ei][hcp];
    mc = (int) gsampinf[ei][mcp];
    fm = (double) gsampinf[ei][fmp];
    qintg = (double) gsampinf[ei][qip];
    mintg = (double) gsampinf[ei][mip];
    if (fc == 0 && cc == 0 && fm > 0)
    {
      temp1 = LOG2 - (mc * log (fm)) - hc - qintg - mintg;
      temp2 = log (exp (uppergamma (mc, 2 * fm * x / qmax)) - exp (uppergamma (mc, mmax * fm)));
    }
    else
    {
      if (fm == 0 && mc == 0 && fc > 0)
      {
        temp1 = LOG2 - (cc * log (fc)) - hc - qintg - mintg;
        temp2 = log (exp (uppergamma (cc, 2 * fc / qmax)) - exp (uppergamma (cc, fc * mmax / x)));
      }
      else
      {
        if (fc == 0 && cc == 0 && mc == 0 && fm == 0)
        {
          temp1 = log (2 * log (mmax * qmax / (2 * x))) - hc - qintg - mintg;;
          temp2 = 0;
        }
        else
        {
          temp1 = LOG2 + (mc * log (x)) - ((cc + mc) * log (fc + fm * x)) - hc - qintg - mintg;
          //temp2 = log ( exp(uppergamma (cc + mc, 2 * (fc + fm * x) / qmax)) - exp(uppergamma (cc + mc, mmax * (fm + fc / x))) );
          a =uppergamma (cc + mc, 2 * (fc + fm * x) / qmax);
          b = uppergamma (cc + mc, mmax * (fm + fc / x));
          /* shortcut to avoid floating point problems that arise when taking the log of a difference  */
          if (b > a) // swap and set C=-1 so that at the end we multiply by negative 1
          {
            c = b;
            b = a;
            a = c;
            c = -1;
          }
          else
            c = 1;
          d = a-b;
          if (a - b < LOG_DBL_MAX)
          {
            temp2 = b + log (exp (d) - 1.0);
          }
          else
          {
            temp2 = a;
          }
          temp2 *= c;


          //temp2 = log (exp (uppergamma (cc + mc, 2 * (fc + fm * x) / qmax)) - exp (uppergamma (cc + mc, mmax * (fm + fc / x))));
        }
      }
    }
    if ((temp1 + temp2 < 700 ) && (temp1 + temp2 > -700 )) // skip things that cannot be exped
      sum += exp (temp1 + temp2);
  }
  sum /= (lasttree - firsttree + (firsttree == 0)); 
  return -sum; 
}                               /* marginpopmig */

#define SEARCHSTARTFRAC  4      // fraction of position in parameter range to start at
/* find the peak of the marginal density for a 2Nm term.  This is similar to the marginalopt() function
  this calls marginpopmig() */ 

void
marginalopt_popmig (int firsttree, int lasttree, double *mlval, double *peakloc, int *mpop, int *mterm)
{
  int i;
  double ftol, ax, bx, cx, fa, fb, fc, xmax, ml;
  double axt, bxt, cxt, prior;
  double max0, min0, max1, min1;
  double (*func) (int, int, int, double, int);
  ftol = 1e-7;
  func = marginpopmig;
  for (i = 0; i < nummigrateparams; i++)
  {
    prior = itheta[mpop[i]].pr.max * imig[mterm[i]].pr.max;
    ax = prior;
    bx = prior / 2;
    mnbrakmod (mterm[i], firsttree, lasttree, &ax, &bx, &cx, &fa, &fb, &fc, func, mpop[i]);
    axt = ax;
    bxt = bx;
    cxt = cx;
    bx = prior / 2;
    ax = MINPARAMVAL;
    mnbrakmod (mterm[i], firsttree, lasttree, &ax, &bx, &cx, &fa, &fb, &fc, func, mpop[i]);
    if (axt < bxt && axt < cxt)
      min0 = axt;
    if (bxt < axt && bxt < cxt)
      min0 = bxt;
    if (cxt < axt && cxt < bxt)
      min0 = cxt;
    if (axt > bxt && axt > cxt)
    {
      if (axt > prior)
        axt = prior;
      max0 = axt;
    }
    if (bxt > axt && bxt > cxt)
    {
      if (bxt > prior)
        bxt = prior;
      max0 = bxt;
    }
    if (cxt > axt && cxt > bxt)
    {
      if (cxt > prior)
        cxt = prior;
      max0 = cxt;
    }
    if (ax < bx && ax < cx)
      min1 = ax;
    if (bx < ax && bx < cx)
      min1 = bx;
    if (cx < ax && cx < bx)
      min1 = cx;
    if (ax > bx && ax > cx)
      max1 = ax;
    {
      if (ax > prior)
        ax = prior;
      max1 = ax;
    }
    if (bx > ax && bx > cx)
    {
      if (bx > prior)
        bx = prior;
      max1 = bx;
    }
    if (cx > ax && cx > bx)
    {
      if (cx > prior)
        cx = prior;
      max1 = cx;
    }
    if (max0 <= min1 || max1 <= min0)
    {
      peakloc[i] = -1;
    }
    else
    {
      ml = -goldenmod (mterm[i], firsttree, lasttree, ax, bx, cx, ftol, &xmax, func,mpop[i]);
      mlval[i] = ml;
      peakloc[i] = xmax;
    }
  }
}                               /* marginalopt_popmig */
