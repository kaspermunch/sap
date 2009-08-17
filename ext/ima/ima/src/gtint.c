/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

/* integrations for assessing if one parameter is greater than another 
print matrices of results */ 
#undef GLOBVARS
#include "imamp.h"

static int cci, ccj, wi, wj;
static int treeinc, hitreenum, numtreesused;
static double fci, fcj, hval, denom, qmax, fmi,fmj, mmax;

static double pgt_fcj_gt_0 (double qi);
static double mgt_wj_gt_0 (double mi);
static double  trapzd(double  (*func)(double ), double  a, double  b, int n);
static double  qtrap(double  (*func)(double ), double  a, double  b);
static double gtmig(int mi, int mj);
static double gtpops(int pj, int pi);
static double logdif (double a,double  b);

#define USETREESMAX  20000

/* shortcut to avoid floating point problems that arise when taking the log of a difference  */
double logdif (double a,double  b)
{
  double c, d, temp;
  if (a==b)
      return 0.0;
assert(a>b);
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
    temp = b + log (exp (d) - 1.0);
  }
  else
  {
    temp = a;
  }
  temp *= c;
  return temp;
}  // logdif


double mgt_wj_gt_0 (double mi)
{
  double temp1, temp2, a,b;
  if (mi < MINPARAMVAL)
    return 0.0;
  a = logfact[wj];
  b = uppergamma(wj+1,fmj*mi);
  temp1 = logdif(a,b);
  //temp1 = log(exp(a)-exp(b));
  temp2 = wi*log(mi) - fmi*mi -(wj+1)*log(fmj) + temp1;
  temp2 -= denom;
  return exp(temp2);
}

double pgt_fcj_gt_0 (double qi)
{
  double temp,temp1, temp2, temp3, temp4, a,b;

  if (qi < MINPARAMVAL)
    return 0.0;
  if (ccj == 0)
  {
    a = uppergamma(0,2*fcj/qi);
    b = uppergamma(0,2*fcj/qmax);
    if (a>b)
    {
      LogDiff(temp,a,b);
      temp = exp(temp);
      temp *= 2*fcj;
    }
    else 
      temp = 2*fcj*(exp(a)-exp(b));
    
    temp1 =  - 2*fci/qi + cci*log(2/qi) ; 
    temp3 = qmax*exp(-2*fcj/qmax) - qi*exp(-2*fcj/qi) + temp;
    temp4  = temp1 + log(temp3) - hval - denom; 
    return exp(temp4);
  }
  else
  {
    a = uppergamma(ccj-1,2*fcj/qmax);
    b = uppergamma(ccj-1,2*fcj/qi);
    if (a>b)
    {
      LogDiff(temp1,a,b);
    }
    else
      return 0.0;
      //temp1 = logdif(a,b);
//    temp1 =log(exp(-uppergamma(ccj-1,2*fcj/qi)) + exp(uppergamma(ccj-1,2*fcj/qmax)) );
    temp2 = LOG2 + cci*log(2/qi) + (1-ccj)*log(fcj) -2*fci/qi; 
    temp3 = temp2 + temp1 - hval - denom;
    return exp(temp3);
  }
}

#define FUNC(x) ((*func)(x))

double  trapzd(double  (*func)(double ), double  a, double  b, int n)
{
   double  x,tnm,sum,del;
   static double  s;
   int it,j;

   if (n == 1) {
      return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
   } else {
      for (it=1,j=1;j<n-1;j++) 
        it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) 
        sum += FUNC(x);
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */


#include <math.h>
#define EPS 1.0e-2
#define JMAX 10

double  qtrap(double  (*func)(double ), double  a, double  b)
{
   /* void nrerror(char error_text[]); */
   int j;
   double  s,olds;

   olds = -1.0e30;
   for (j=1;j<=JMAX;j++) {
      s=trapzd(func,a,b,j);
      /* less stringent with it looks like there is very little area under the curve */ 
      if (s< 0.0001/treessaved)
      {
        if (fabs(s-olds) < fabs(olds)) 
          return s;
      }
      else
      {
        if (fabs(s-olds) < EPS*fabs(olds)) 
          return s;
      }
      olds=s;
   }
   return s;
   nrerror("Too many steps in routine qtrap");
   return 0.0;
}
#undef EPS
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */



double gtmig(int mi, int mj)
{
  int ei;
  double sum, temp, temp1, temp2, temp3, temp4, a,b;
  double  (*func)(double );

  func = mgt_wj_gt_0;

  mmax = imig[mi].pr.max;

  sum = 0;
  //for (ei = 0; ei < treessaved; ei++)
  for (ei = 0; ei < hitreenum; ei+= treeinc)
  {
    wi = (int) gsampinf[ei][gsamp_mcp + mi];
    fmi = gsampinf[ei][gsamp_fmp + mi];
    wj = (int) gsampinf[ei][gsamp_mcp + mj];
    fmj = gsampinf[ei][gsamp_fmp + mj];
    denom = gsampinf[ei][gsamp_mip + mi]+gsampinf[ei][gsamp_mip + mj];

    if (wj==0)
    {
      if (fmj > 0)
      {
        a = logfact[wi];
        b = uppergamma(wi+1,(fmi + fmj)*mmax);
        assert(a>b);
        //temp1 = logdif(a,b);
        LogDiff(temp1,a,b);
        //temp1 = log(exp(a)-exp(b));
        temp1 += -(wi+1)*log(fmi + fmj);
        b = uppergamma(wi+1,fmi*mmax);
        assert(a>b);
        LogDiff(temp2,a,b);
        //temp2 = logdif(a,b);
        //temp2 = log(exp(a)-exp(b));
        temp2 += -(wi+1)*log(fmi);
        assert(temp2>temp1);
        LogDiff(temp3,temp2,temp1);
        temp3 = logdif(temp2,temp1);
        //temp3 = log(exp(temp1)-exp(temp2));
        temp4 = temp3 - log(fmj) - denom;
        temp = exp(temp4);
      }
      else
      {
        a = logfact[wi+1];
        b = uppergamma(wi+2,fmi*mmax);
        //temp1 = log(exp(a)-exp(b));
        assert(a>b);
        LogDiff(temp1,a,b);
        //temp1 = logdif(a,b);
        temp2 = -(2+wi)*log(fmi);
        temp3 = temp1 + temp2 - denom;
        temp = exp(temp3);
      }
    }
    else
    {
      temp = qtrap(func, MINPARAMVAL, mmax);
    }
    sum += temp;
  }
  sum /= numtreesused;
  return sum;

} // gtmig

/* return the probability that pram pj < param pi */
double gtpops(int pj, int pi)
{
  int ei;
  double sum, temp, temp1, temp2, temp3;
  double  (*func)(double );

  func = pgt_fcj_gt_0;

  qmax = itheta[pi].pr.max;

  sum = 0;
  //for (ei = 0; ei < treessaved; ei++)
  for (ei = 0; ei < hitreenum; ei+= treeinc)
  {
    cci = (int) gsampinf[ei][gsamp_ccp + pi];
    ccj = (int) gsampinf[ei][gsamp_ccp + pj];
    fci = gsampinf[ei][gsamp_fcp + pi];
    fcj = gsampinf[ei][gsamp_fcp + pj];
    hval = gsampinf[ei][gsamp_hccp + pi] + gsampinf[ei][gsamp_hccp + pj];
    denom = gsampinf[ei][gsamp_qip + pi]+gsampinf[ei][gsamp_qip + pj];

    if (ccj==0)
    {
      if (fcj == 0)
      {
        if (fci==0)
        {
          temp = exp(2.0 * log(qmax)- LOG2 - hval - denom);
        }
        else
        {
          if (cci>= 2)
          {
            temp1 = LOG2 + (1-cci)* log(fci);
            temp2 = log(qmax* exp(uppergamma(cci-1,2*fci/qmax)) - 2*fci*exp(uppergamma(cci-2,2*fci/qmax)));
            temp3 = temp1 + temp2 - hval-denom;
            temp = exp(temp3);
          }
          else
          {
            if (cci == 1)
            {
              temp1 = (2*fci + qmax)*exp(uppergamma(0,2*fci/qmax));
              temp2 = temp1 -qmax*exp(-2*fci/qmax);
              temp = LOG2 + log(temp2) - hval - denom;
            }
            else  // cci==0
            {
              temp1 =   2*LOG2 + log(fci) + log(qmax+fci) + uppergamma(0,2*fci/qmax);  // use for ei() of a negative value 
              temp2 = log(exp(-2*fci/qmax) * qmax * (2*fci + qmax));
              temp3 = temp1 + temp2- hval -denom;
              temp = temp3;
            }
          }
        }
      }
      else //ccj=0 fcj > 0 
      {
        temp = qtrap(func, MINPARAMVAL, qmax);

      }
    }
    else
    {
      temp = qtrap(func, MINPARAMVAL, qmax);
    }
    sum += temp;
  }
sum /= numtreesused;
  return sum;
} // gtpops 

void print_greater_than_tests (FILE * outfile)
{
  double **gt_popsize, **gt_mig;
  int i, j;
  if (treessaved > USETREESMAX)
  {
    treeinc = (int) treessaved/(int) USETREESMAX;
    numtreesused = USETREESMAX;
    hitreenum = treeinc * USETREESMAX;
  }
  else
  {
    treeinc = 1;
    hitreenum = numtreesused = treessaved;
  }
  gt_popsize= alloc2Ddouble (numpopsizeparams, numpopsizeparams);
  gt_mig = alloc2Ddouble (nummigrateparams, nummigrateparams);

  for (i=0;i < numpopsizeparams - 1; i++)
    for (j = i+1; j < numpopsizeparams; j++)
      gt_popsize[i][j] = gtpops(i,j); 
  for (i=0;i < nummigrateparams - 1; i++)
    for (j = i+1; j < nummigrateparams; j++)
      gt_mig[i][j] = gtmig(i,j);

  FP "\nPARAMETER COMPARISONS, PROBABILITY THAT ROW PARAMETER IS GREATER THAN COLUMN PARAMETER\n");
  FP "========================================================================================\n");
  FP "Population Sizes\n");
  for (i=0;i<numpopsizeparams;i++)
    FP "\t%s", itheta[i].str);
  FP "\n");
  for (i=0;i<numpopsizeparams;i++)
  {
    FP "%s", itheta[i].str);
    for (j=0;j<numpopsizeparams;j++)
    {
      if (i == j)
        FP "\t  - ");
      else
      {
        if (j > i)
        {
          FP "\t%.3lf", gt_popsize[i][j]);
        }
        else
        {
          FP "\t%.3lf", 1-gt_popsize[j][i]);
        }
      }
    }
    FP "\n");
  }
  FP "\n");
  
  FP "Migration Rates\n");
  for (i=0;i<nummigrateparams;i++)
  {
    FP "\t%s", imig[i].str);
  }
  FP "\n");
  for (i=0;i<nummigrateparams;i++)
  {
    FP "%s", imig[i].str);
    for (j=0;j<nummigrateparams;j++)
    {
      if (i == j)
        FP "\t  - ");
      else
      {
        if (j > i)
        {
          FP "\t%.3lf", gt_mig[i][j]);
        }
        else
        {
          FP "\t%.3lf", 1-gt_mig[j][i]);
        }
      }
    }
    FP "\n");
  } 
  FP "\n\n");
  free2D ((void **) gt_popsize, numpopsizeparams);
  free2D ((void **) gt_mig, nummigrateparams);
  return;
}                               // print_greater_than_tests
#undef USETREESMAX 
