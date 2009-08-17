/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#ifndef _UTILITIES_H_
#define _UTILITIES_H_
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

#include <stdio.h>

__BEGIN_DECLS

enum
{
  IMERR_SUCCESS = 0,
  IMERR_READFILEOPENFAIL = 1,
  IMERR_MEM = 2,
  IMERR_TIFILE = 3,
  IMERR_CREATEFILEFAIL = 4,
  IMERR_APPENDFILEFAIL = 5,
  IMERR_OUTPUTFILECHECK = 6,

  IMERR_COMMANDLINEINCOMPAT = 8,
  IMERR_MISSINGCOMMANDINFO = 9,
  IMERR_COMMANDLINEFORMAT = 10,

  IMERR_MISSINGPOPSTRING = 12,
  IMERR_POPTREESTRINGFAIL = 13,

  IMERR_INFILEFAIL_NLOCI = 15,

  IMERR_MUTSCALARPRIORRANGEFAIL = 18,
  IMERR_MUTSCALEPRODUCTFAIL = 19,

  IMERR_TOOMANYMIG = 21,
  IMERR_MIGARRAYTOOBIG = 22,
  IMERR_HPD95 = 23,

  IMERR_DATAREADOVERRUN = 25,
  IMERR_DATAERROR = 26,

  IMERR_MCFREADFAIL = 28,
  IMERR_MCFSPLITTIMEPROB = 29,

  IMERR_ROOTTIMEMAXFAIL = 33,

  IMERR_INPUTFILEINVALID = 35,
  IMERR_INFINITESITESFAIL = 36,
  
  IMERR_SWCHECK = 38,

  IMERR_LOWERGAMMA = 41,
  IMERR_UPPERGAMMA = 42,

  IMERR_MULTITPRIOR = 45,

  IMERR_NUMERICALRECIPES = 51,
  IMERR_GSL = 52,
  IMERR_ASSERT = 55,
  IMERR_GENENAME = 56,
  IMERR_ASN = 57
};

void IM_err (int i, const char *fmt, ...);
/* void skip_a_line (FILE * f); */
int comma_exists (char *s);
double logsum (int n, ...);
double logsuma (int n, double *d); 
double logdiff (double a, double b);
double cumnorm(double ox, double mean, double sd);

__END_DECLS
#endif /* _UTILITIES_H_ */
