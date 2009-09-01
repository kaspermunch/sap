#ifndef _MSC_VER
#include "config.h"
#endif /* _MSC_VER */
#include "imagsl.h"
/* err/test_results.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
/* imagsl gsl sys h: SANGCHUL - this is added in imagsl.h */
/* imagsl gsl machine h: SANGCHUL - this is the first header file added. */

#if HAVE_VPRINTF
#ifdef STDC_HEADERS
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#endif

/* imagsl gsl test h: SANGCHUL - this is added in imagsl.h */

static unsigned int tests = 0;
static unsigned int passed = 0;
static unsigned int failed = 0;

static unsigned int verbose = 0;

static void
initialise (void)
{
  const char * p = getenv("GSL_TEST_VERBOSE");

  /* 0 = show failures only (we always want to see these) */
  /* 1 = show passes and failures */

  if (p == 0)  /* environment variable is not set */
    return ;

  if (*p == '\0') /* environment variable is empty */
    return ;

  verbose = strtoul (p, 0, 0);  

  return;
}

static void 
update (int s)
{
  tests++;

  if (s == 0) 
    {
      passed++;
    }
  else
    {
      failed++;
    }
}

void
gsl_test (int status, const char *test_description,...)
{
  if (!tests) initialise();

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif

      if (status && !verbose)
        printf(" [%u]", tests);

      printf("\n");
      fflush (stdout);
    }
}


void
gsl_test_rel (double result, double expected, double relative_error,
              const char *test_description,...)
{
  int status ;

  if (!tests) initialise();

  /* Check for NaN vs inf vs number */

  if (gsl_isnan(result) || gsl_isnan(expected)) 
    {
      status = gsl_isnan(result) != gsl_isnan(expected); 
    }
  else if (gsl_isinf(result) || gsl_isinf(expected)) 
    {
      status = gsl_isinf(result) != gsl_isinf(expected); 
    }
  else if ((expected > 0 && expected < GSL_DBL_MIN)
           || (expected < 0 && expected > -(GSL_DBL_MIN)))
    {
      status = -1;
    }
  else if (expected != 0 ) 
    {
      status = (fabs(result-expected)/fabs(expected) > relative_error) ;
    }
  else
    {
      status = (fabs(result) > relative_error) ;
    }

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif

      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n") ;
      fflush (stdout);
    }
}

void
gsl_test_abs (double result, double expected, double absolute_error,
              const char *test_description,...)
{
  int status ;

  if (!tests) initialise();

  /* Check for NaN vs inf vs number */

  if (gsl_isnan(result) || gsl_isnan(expected)) 
    {
      status = gsl_isnan(result) != gsl_isnan(expected); 
    }
  else if (gsl_isinf(result) || gsl_isinf(expected)) 
    {
      status = gsl_isinf(result) != gsl_isinf(expected); 
    }
  else if ((expected > 0 && expected < GSL_DBL_MIN)
           || (expected < 0 && expected > -(GSL_DBL_MIN)))
    {
      status = -1;
    }
  else 
    {
      status = fabs(result-expected) > absolute_error ;
    }

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif

      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n") ;
      fflush (stdout);
    }
}


void
gsl_test_factor (double result, double expected, double factor,
                 const char *test_description,...)
{
  int status;

  if (!tests) initialise();
  
  if ((expected > 0 && expected < GSL_DBL_MIN)
      || (expected < 0 && expected > -(GSL_DBL_MIN)))
    {
      status = -1;
    }
  else if (result == expected) 
    {
      status = 0;
    }
  else if (expected == 0.0) 
    {
      status = (result > expected || result < expected);
    }
  else
    {
      double u = result / expected; 
      status = (u > factor || u < 1.0 / factor) ;
    }

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status == 0)
        {
          if (strlen(test_description) < 45)
            {
              printf(" (%g observed vs %g expected)", result, expected) ;
            }
          else
            {
              printf(" (%g obs vs %g exp)", result, expected) ;
            }
        }
      else 
        {
          printf(" (%.18g observed vs %.18g expected)", result, expected) ;
        }

      if (status == -1)
        {
          printf(" [test uses subnormal value]") ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n") ;
      fflush (stdout);
    }
}

void
gsl_test_int (int result, int expected, const char *test_description,...)
{
  int status = (result != expected) ;

  if (!tests) initialise();

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status == 0)
        {
          printf(" (%d observed vs %d expected)", result, expected) ;
        }
      else 
        {
          printf(" (%d observed vs %d expected)", result, expected) ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n");
      fflush (stdout);
    }
}

void
gsl_test_str (const char * result, const char * expected, 
              const char *test_description,...)
{
  int status = strcmp(result,expected) ;

  if (!tests) initialise();

  update (status);

  if (status || verbose)
    {
      printf (status ? "FAIL: " : "PASS: ");

#if HAVE_VPRINTF
      {
        va_list ap;
        
#ifdef STDC_HEADERS
        va_start (ap, test_description);
#else
        va_start (ap);
#endif
        vprintf (test_description, ap);
        va_end (ap);
      }
#endif
      if (status)
        {
          printf(" (%s observed vs %s expected)", result, expected) ;
        }

      if (status && !verbose)
        printf(" [%u]", tests);

      printf ("\n");
      fflush (stdout);
    }
}

void
gsl_test_verbose (int v)
{
  verbose = v;
}

int
gsl_test_summary (void)
{
  if (verbose && 0)             /* FIXME: turned it off, this annoys me */
    printf ("%d tests, passed %d, failed %d.\n", tests, passed, failed);

  if (failed != 0)
    {
      return EXIT_FAILURE;
    }

  if (tests != passed + failed)
    {
      if (verbose)
        printf ("TEST RESULTS DO NOT ADD UP %d != %d + %d\n",
                tests, passed, failed);
      return EXIT_FAILURE;
    }

  if (passed == tests)
    {
      if (!verbose)         /* display a summary of passed tests */
        printf ("Completed [%d/%d]\n", passed, tests);

      return EXIT_SUCCESS;
    }

  return EXIT_FAILURE;
}
/* sys/fdiv.c
 * 
 * Copyright (C) 2001, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>

double 
gsl_fdiv (const double x, const double y)
{
  return x / y;
}
/* sys/infnan.c
 * 
 * Copyright (C) 2001, 2004, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>

#if HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

double gsl_nan (void)
{
  return gsl_fdiv (0.0, 0.0);
}

double gsl_posinf (void)
{
  return gsl_fdiv (+1.0, 0.0);
}

double gsl_neginf (void)
{
  return gsl_fdiv (-1.0, 0.0);
}


int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

#if defined(_MSC_VER) /* Microsoft Visual C++ */
#include <float.h>
int
gsl_isnan (const double x)
{
  return _isnan(x);
}

int
gsl_isinf (const double x)
{
  int fpc = _fpclass(x);

  if (fpc == _FPCLASS_PINF)
    return +1;
  else if (fpc == _FPCLASS_NINF)
    return -1;
  else 
    return 0;
}

int
gsl_finite (const double x)
{
  return _finite(x);
}
#else

# if HAVE_DECL_ISFINITE
int
gsl_finite (const double x)
{
  return isfinite(x);
}
# elif HAVE_DECL_FINITE
int
gsl_finite (const double x)
{
  return finite(x);
}
# elif HAVE_IEEE_COMPARISONS
int
gsl_finite (const double x)
{
  const double y = x - x;
  int status = (y == y);
  return status;
}
# endif

# if HAVE_DECL_ISNAN
int
gsl_isnan (const double x)
{
  return isnan(x);
}
#elif HAVE_IEEE_COMPARISONS
int
gsl_isnan (const double x)
{
  int status = (x != x);
  return status;
}
# endif

# if HAVE_DECL_ISINF
int
gsl_isinf (const double x)
{
  /* isinf(3): In glibc 2.01 and earlier, isinf() returns a
     non-zero value (actually: 1) if x is an infinity (positive or
     negative).  (This is all that C99 requires.) */

  if (isinf(x)) 
    {
      return (x > 0) ? 1 : -1;
    } 
  else 
    {
      return 0;
    }
}
# else

int
gsl_isinf (const double x)
{
  if (! gsl_finite(x) && ! gsl_isnan(x)) 
    {
      return (x > 0 ? +1 : -1); 
    } 
  else 
    {
      return 0;
    }
}

# endif
#endif



/* err/error.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

gsl_error_handler_t * gsl_error_handler = NULL;

static void no_error_handler (const char *reason, const char *file, int line, int gsl_errno);

void
gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_error_handler) 
    {
      (*gsl_error_handler) (reason, file, line, gsl_errno);
      return ;
    }

  gsl_stream_printf ("ERROR", file, line, reason);

  fflush (stdout);
  fprintf (stderr, "Default GSL error handler invoked.\n");
  fflush (stderr);

  abort ();
}

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler)
{
  gsl_error_handler_t * previous_handler = gsl_error_handler;
  gsl_error_handler = new_handler;
  return previous_handler;
}


gsl_error_handler_t *
gsl_set_error_handler_off (void)
{
  gsl_error_handler_t * previous_handler = gsl_error_handler;
  gsl_error_handler = no_error_handler;
  return previous_handler;
}

static void
no_error_handler (const char *reason, const char *file, int line, int gsl_errno)
{
  /* do nothing */
  reason = 0;
  file = 0;
  line = 0;
  gsl_errno = 0;
  return;
}


/* err/message.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

unsigned int gsl_message_mask = GSL_MESSAGE_MASK;

void
gsl_message (const char * reason, const char * file, int line, 
             unsigned int mask)
{
  if (mask & gsl_message_mask)
    {
      gsl_stream_printf ("MESSAGE", file, line, reason);
    }
}
/* err/stream.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

FILE * gsl_stream = NULL ;
gsl_stream_handler_t * gsl_stream_handler = NULL;

void
gsl_stream_printf (const char *label, const char *file, int line, 
                   const char *reason)
{
  if (gsl_stream == NULL)
    {
      gsl_stream = stderr;
    }
  if (gsl_stream_handler)
    {
      (*gsl_stream_handler) (label, file, line, reason);
      return;
    }
  fprintf (gsl_stream, "gsl: %s:%d: %s: %s\n", file, line, label, reason);

}

gsl_stream_handler_t *
gsl_set_stream_handler (gsl_stream_handler_t * new_handler)
{
  gsl_stream_handler_t * previous_handler = gsl_stream_handler;
  gsl_stream_handler = new_handler;
  return previous_handler;
}

FILE *
gsl_set_stream (FILE * new_stream)
{
  FILE * previous_stream;
  if (gsl_stream == NULL) {
    gsl_stream = stderr;
  }
  previous_stream = gsl_stream;
  gsl_stream = new_stream;
  return previous_stream;
}
/* err/strerror.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */

const char *
gsl_strerror (const int gsl_errno)
{
  switch (gsl_errno)
    {
    case GSL_SUCCESS:
      return "success" ;
    case GSL_FAILURE:
      return "failure" ;
    case GSL_CONTINUE:
      return "the iteration has not converged yet";
    case GSL_EDOM:
      return "input domain error" ;
    case GSL_ERANGE:
      return "output range error" ;
    case GSL_EFAULT:
      return "invalid pointer" ;
    case GSL_EINVAL:
      return "invalid argument supplied by user" ;
    case GSL_EFAILED:
      return "generic failure" ;
    case GSL_EFACTOR:
      return "factorization failed" ;
    case GSL_ESANITY:
      return "sanity check failed - shouldn't happen" ;
    case GSL_ENOMEM:
      return "malloc failed" ;
    case GSL_EBADFUNC:
      return "problem with user-supplied function";
    case GSL_ERUNAWAY:
      return "iterative process is out of control";
    case GSL_EMAXITER:
      return "exceeded max number of iterations" ;
    case GSL_EZERODIV:
      return "tried to divide by zero" ;
    case GSL_EBADTOL:
      return "specified tolerance is invalid or theoretically unattainable" ;
    case GSL_ETOL:
      return "failed to reach the specified tolerance" ;
    case GSL_EUNDRFLW:
      return "underflow" ;
    case GSL_EOVRFLW:
      return "overflow" ;
    case GSL_ELOSS:
      return "loss of accuracy" ;
    case GSL_EROUND:
      return "roundoff error" ;
    case GSL_EBADLEN:
      return "matrix/vector sizes are not conformant" ;
    case GSL_ENOTSQR:
      return "matrix not square" ;
    case GSL_ESING:
      return "singularity or extremely bad function behavior detected" ;
    case GSL_EDIVERGE:
      return "integral or series is divergent" ;
    case GSL_EUNSUP:
      return "the required feature is not supported by this hardware platform";
    case GSL_EUNIMPL:
      return "the requested feature is not (yet) implemented";
    case GSL_ECACHE:
      return "cache limit exceeded";
    case GSL_ETABLE:
      return "table limit exceeded";
    case GSL_ENOPROG:
      return "iteration is not making progress towards solution";
    case GSL_ENOPROGJ:
      return "jacobian evaluations are not improving the solution";
    case GSL_ETOLF:
      return "cannot reach the specified tolerance in F";
    case GSL_ETOLX:
      return "cannot reach the specified tolerance in X";
    case GSL_ETOLG:
      return "cannot reach the specified tolerance in gradient";
    case GSL_EOF:
      return "end of file";
    default:
      return "unknown error code" ;
    }
}
/* ieee-utils/env.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <stdlib.h>

void
gsl_ieee_env_setup (void)
{
  const char * p = getenv("GSL_IEEE_MODE") ;

  int precision = 0, rounding = 0, exception_mask = 0 ;

  int comma = 0 ;

  if (p == 0)  /* GSL_IEEE_MODE environment variable is not set */
    return ;

  if (*p == '\0') /* GSL_IEEE_MODE environment variable is empty */
    return ;

  gsl_ieee_read_mode_string (p, &precision, &rounding, &exception_mask) ;

  gsl_ieee_set_mode (precision, rounding, exception_mask) ;
  
  fprintf(stderr, "GSL_IEEE_MODE=\"") ;

  /* Print string with a preceeding comma if the list has already begun */

#define PRINTC(x) do {if(comma) fprintf(stderr,","); fprintf(stderr,x); comma++ ;} while(0)
  
  switch (precision) 
    {
    case GSL_IEEE_SINGLE_PRECISION:
      PRINTC("single-precision") ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      PRINTC("double-precision") ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      PRINTC("extended-precision") ;
      break ;
    }

  switch (rounding) 
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      PRINTC("round-to-nearest") ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      PRINTC("round-down") ;
      break ;
    case GSL_IEEE_ROUND_UP:
      PRINTC("round-up") ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      PRINTC("round-to-zero") ;
      break ;
    }

  if ((exception_mask & GSL_IEEE_MASK_ALL) == GSL_IEEE_MASK_ALL)
    {
      PRINTC("mask-all") ;
    }
  else if ((exception_mask & GSL_IEEE_MASK_ALL) == 0)
    {
      PRINTC("trap-common") ;
    }
  else 
    {
      if (exception_mask & GSL_IEEE_MASK_INVALID)
        PRINTC("mask-invalid") ;
      
      if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
        PRINTC("mask-denormalized") ;
      
      if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
        PRINTC("mask-division-by-zero") ;
      
      if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
        PRINTC("mask-overflow") ;
      
      if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
        PRINTC("mask-underflow") ;
    }

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    PRINTC("trap-inexact") ;
  
  fprintf(stderr,"\"\n") ;
}





/* #include <config.h> */

#if HAVE_GNUSPARC_IEEE_INTERFACE
/* ieee-utils/fp-gnusparc.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <fpu_control.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  unsigned short mode = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      mode |= _FPU_SINGLE ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      mode |= _FPU_DOUBLE ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      mode |= _FPU_EXTENDED ;
      break ;
    default:
      mode |= _FPU_EXTENDED ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case GSL_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_IM ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("sparc does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.", GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_ZM ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OM ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UM ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode &= ~ _FPU_MASK_PM ;
    }
  else
    {
      mode |= _FPU_MASK_PM ;
    }

  _FPU_SETCW(mode) ;

  return GSL_SUCCESS ;
}
#elif HAVE_GNUM68K_IEEE_INTERFACE
/* ieee-utils/fp-gnum68k.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <fpu_control.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  unsigned short mode = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      mode |= _FPU_SINGLE ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      mode |= _FPU_DOUBLE ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      mode |= _FPU_EXTENDED ;
      break ;
    default:
      mode |= _FPU_EXTENDED ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case GSL_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  /* FIXME: I don't have documentation for the M68K so I'm not sure
     about the mapping of the exceptions below. Maybe someone who does
     know could correct this. */

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_OPERR ;
  
  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    } 
  else
    {
      GSL_ERROR ("the denormalized operand exception has not been implemented for m68k yet. Use 'mask-denormalized' to work around this.", GSL_EUNSUP) ;
      /*mode |= _FPU_MASK_DM ; ???? */ 
    }
  
  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_DZ ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OVFL ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UNFL ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode &= ~ (_FPU_MASK_INEX1 | _FPU_MASK_INEX2) ;
    }
  else
    {
      mode |= (_FPU_MASK_INEX1 | _FPU_MASK_INEX2) ;
    }

  _FPU_SETCW(mode) ;

  return GSL_SUCCESS ;
}
#elif HAVE_GNUPPC_IEEE_INTERFACE
/* ieee-utils/fp-gnuppc.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough, John Fisher
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <fpu_control.h>


/*
 * Identical to fp-gnux86.c, except with references to
 * _FPU_SINGLE, _FPU_DOUBLE, _FPU_EXTENDED, _FPU_MASK_DM
 * and _FPU_MASK_PM converted to errors.
 */

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  unsigned short mode = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case GSL_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_IM ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("powerpc does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.", GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_ZM ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OM ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UM ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
     GSL_ERROR ("powerpc does not support traps for inexact operations", GSL_EUNSUP) ;
    }

  _FPU_SETCW(mode) ;

  return GSL_SUCCESS ;
}
#elif HAVE_GNUX86_IEEE_INTERFACE
/* ieee-utils/fp-gnux86.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <fpu_control.h>

  /* Handle libc5, where _FPU_SETCW is not available, suggested by
     OKUJI Yoshinori <okuji@gnu.org> and Evgeny Stambulchik
     <fnevgeny@plasma-gate.weizmann.ac.il> */

#ifndef _FPU_SETCW
#include <i386/fpu_control.h>
#define _FPU_SETCW(cw) __setfpucw(cw)
#endif

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  unsigned short mode = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      mode |= _FPU_SINGLE ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      mode |= _FPU_DOUBLE ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      mode |= _FPU_EXTENDED ;
      break ;
    default:
      mode |= _FPU_EXTENDED ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case GSL_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_IM ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode |= _FPU_MASK_DM ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_ZM ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OM ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UM ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode &= ~ _FPU_MASK_PM ;
    }
  else
    {
      mode |= _FPU_MASK_PM ;
    }

  _FPU_SETCW(mode) ;

#if HAVE_FPU_X86_SSE
#define _FPU_SETMXCSR(cw_sse) asm volatile ("ldmxcsr %0" : : "m" (*&cw_sse))
  {
    unsigned int mode_sse = 0;

    mode_sse |= (mode & 0x3f)<<7;  /* exception masks */
    mode_sse |= (mode & 0xc00)<<3;    /* rounding control */
    
    _FPU_SETMXCSR(mode_sse);
  }
#endif

  return GSL_SUCCESS ;
}
#elif HAVE_HPUX11_IEEE_INTERFACE
/* ieee-utils/fp-hpux11.c
 * 
 * Copyright (C) 2001, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <stdio.h>
#include <fenv.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  int mode;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("HPUX PA-RISC only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("HPUX PA-RISC only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("HPUX PA-RISC only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    }


  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      fesetround (FE_TONEAREST) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      fesetround (FE_DOWNWARD) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      fesetround (FE_UPWARD) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      fesetround (FE_TOWARDZERO) ;
      break ;
    default:
      fesetround (FE_TONEAREST) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ FE_INVALID ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("HP-UX does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.",
                 GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ FE_DIVBYZERO ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ FE_OVERFLOW ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ FE_UNDERFLOW ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= FE_INEXACT ;
    }
  else
    {
      mode &= ~ FE_INEXACT ;
    }

  fesettrapenable (mode) ;

  return GSL_SUCCESS ;
}
#elif HAVE_HPUX_IEEE_INTERFACE
/* ieee-utils/fp-hpux.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <stdio.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fp_except mode = 0 ;
  fp_rnd    rnd  = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("HPUX PA-RISC only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("HPUX PA-RISC only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("HPUX PA-RISC only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    }


  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RN ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RM ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RP ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RZ ;
      fpsetround (rnd) ;
      break ;
    default:
      rnd = FP_RN ;
      fpsetround (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = FP_X_INV | FP_X_DZ | FP_X_OFL | FP_X_UFL ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ FP_X_INV ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("HP-UX does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.",
                 GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ FP_X_DZ ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ FP_X_OFL ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ FP_X_UFL ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= FP_X_IMP ;
    }
  else
    {
      mode &= ~ FP_X_IMP ;
    }

  fpsetmask (mode) ;

  return GSL_SUCCESS ;
}
#elif HAVE_SUNOS4_IEEE_INTERFACE
/* ieee-utils/fp-sunos4.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <sys/ieeefp.h>
#include <floatingpoint.h>
#include <signal.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  char * out ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      ieee_flags ("set", "precision", "single", out) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      ieee_flags ("set", "precision", "double", out) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      ieee_flags ("set", "precision", "extended", out) ;
      break ;
    default:
      ieee_flags ("set", "precision", "extended", out) ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      ieee_flags ("set", "direction", "nearest", out) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      ieee_flags ("set", "direction", "negative", out) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      ieee_flags ("set", "direction", "positive", out) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      ieee_flags ("set", "direction", "tozero", out) ;
      break ;
    default:
      ieee_flags ("set", "direction", "nearest", out) ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    {
      ieee_handler ("set", "invalid", SIGFPE_IGNORE) ;
    }
  else 
    {
      ieee_handler ("set", "invalid", SIGFPE_ABORT) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      ieee_handler ("set", "denormalized", SIGFPE_IGNORE) ;
    }
  else
    {
      GSL_ERROR ("sunos4 does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.",
                 GSL_EUNSUP) ;
    }


  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    {
      ieee_handler ("set", "division", SIGFPE_IGNORE) ;
    } 
  else
    {
      ieee_handler ("set", "division", SIGFPE_ABORT) ;
    }
  
  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    {
      ieee_handler ("set", "overflow", SIGFPE_IGNORE) ;
    }
  else 
    {
      ieee_handler ("set", "overflow", SIGFPE_ABORT) ;
    }

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    {
      ieee_handler ("set", "underflow", SIGFPE_IGNORE) ;
    }
  else
    {
      ieee_handler ("set", "underflow", SIGFPE_ABORT) ;
    }

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      ieee_handler ("set", "inexact", SIGFPE_ABORT) ;
    }
  else
    {
      ieee_handler ("set", "inexact", SIGFPE_IGNORE) ;
    }

  return GSL_SUCCESS ;
}
#elif HAVE_SOLARIS_IEEE_INTERFACE
/* ieee-utils/fp-solaris.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <ieeefp.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fp_except mode = 0 ;
  fp_rnd    rnd  = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("solaris only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("solaris only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("solaris only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RN ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RM ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RP ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RZ ;
      fpsetround (rnd) ;
      break ;
    default:
      rnd = FP_RN ;
      fpsetround (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = FP_X_INV | FP_X_DZ | FP_X_OFL | FP_X_UFL ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ FP_X_INV ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("solaris does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.",
                 GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ FP_X_DZ ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ FP_X_OFL ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ FP_X_UFL ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= FP_X_IMP ;
    }
  else
    {
      mode &= ~ FP_X_IMP ;
    }

  fpsetmask (mode) ;

  return GSL_SUCCESS ;

}
#elif HAVE_IRIX_IEEE_INTERFACE
/* ieee-utils/fp-irix.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Tim Mooney
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <ieeefp.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fp_except mode = 0 ;
  fp_rnd    rnd  = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("IRIX only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("IRIX only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("IRIX only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RN ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RM ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RP ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RZ ;
      fpsetround (rnd) ;
      break ;
    default:
      rnd = FP_RN ;
      fpsetround (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = FP_X_INV | FP_X_DZ | FP_X_OFL | FP_X_UFL ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ FP_X_INV ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("IRIX does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.",
                 GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ FP_X_DZ ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ FP_X_OFL ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ FP_X_UFL ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= FP_X_IMP ;
    }
  else
    {
      mode &= ~ FP_X_IMP ;
    }

  fpsetmask (mode) ;

  return GSL_SUCCESS ;

}
#elif HAVE_AIX_IEEE_INTERFACE
/* ieee-utils/fp-aix.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Tim Mooney
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <fptrap.h>
#include <float.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fptrap_t   mode = 0 ;
  fprnd_t    rnd  = 0 ;

  switch (precision)
    {

    /* I'm not positive about AIX only supporting default precision rounding,
     * but this is the best assumption until it's proven otherwise. */

    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("AIX only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("AIX only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("AIX only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RND_RN ;
      fp_swap_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RND_RM ;
      fp_swap_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RND_RP ;
      fp_swap_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RND_RZ ;
      fp_swap_rnd (rnd) ;
      break ;
    default:
      rnd = FP_RND_RN ;
      fp_swap_rnd (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = TRP_INVALID | TRP_DIV_BY_ZERO | TRP_OVERFLOW | TRP_UNDERFLOW ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ TRP_INVALID ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else 
    {
      GSL_ERROR ("AIX does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.",
                 GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ TRP_DIV_BY_ZERO ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ TRP_OVERFLOW ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ TRP_UNDERFLOW ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= TRP_INEXACT ;
    }
  else
    {
      mode &= ~ TRP_INEXACT ;
    }

  /* AIX appears to require two steps -- first enable floating point traps
   * in general... */
  fp_trap(FP_TRAP_SYNC);

  /* next, enable the traps we're interested in */
  fp_enable(mode);

  return GSL_SUCCESS ;

}
#elif HAVE_TRU64_IEEE_INTERFACE
/* ieee-utils/fp-tru64.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Tim Mooney
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


/*
 * Under Compaq's Unix with the silly name, read the man pages for read_rnd,
 * write_rnd, and ieee(3) for more information on the functions used here.
 *
 * Note that enabling control of dynamic rounding mode (via write_rnd) requires
 * that you pass a special flag to your C compiler.  For Compaq's C compiler
 * the flag is `-fprm d', for gcc it's `-mfp-rounding-mode=d'.
 *
 * Enabling the trap control (via ieee_set_fp_control) also requires a
 * flag be passed to the C compiler.  The flag for Compaq's C compiler
 * is `-ieee' and for gcc it's `-mieee'.

 * We have not implemented the `inexact' case, since it is rarely used
 * and requires the library being built with an additional compiler
 * flag that can degrade performance for everything else. If you need
 * to add support for `inexact' the relevant flag for Compaq's
 * compiler is `-ieee_with_inexact', and the flag for gcc is
 * `-mieee-with-inexact'.
 *
 * Problem have been reported with the "fixed" float.h installed with
 * gcc-2.95 lacking some of the definitions in the system float.h (the
 * symptoms are errors like: `FP_RND_RN' undeclared). To work around
 * this we can include the system float.h before the gcc version, e.g. 
 *
 *  #include "/usr/include/float.h"
 *  #include <float.h>
 */

#include <float.h>

#ifndef FP_RND_RN
#  undef _FLOAT_H_
#  include "/usr/include/float.h"
#  undef _FLOAT_H_
#  include <float.h>
#endif

#include <machine/fpu.h>
#include <stdio.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  unsigned long int mode = 0 ;
  unsigned int    rnd  = 0 ;

/* I'm actually not completely sure that the alpha only supports default
 * precisions rounding, but I couldn't find any information regarding this, so
 * it seems safe to assume this for now until it's proven otherwise.
 */

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("Tru64 Unix on the alpha only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("Tru64 Unix on the alpha only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("Tru64 Unix on the alpha only supports default precision rounding",
                 GSL_EUNSUP) ;
      break ;
    }


  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RND_RN ;
      write_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RND_RM ;
      write_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RND_RP ;
      write_rnd (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RND_RZ ;
      write_rnd (rnd) ;
      break ;
    default:
      rnd = FP_RND_RN ;
      write_rnd (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  /* from the ieee(3) man page:
   * IEEE_TRAP_ENABLE_INV       ->      Invalid operation
   * IEEE_TRAP_ENABLE_DZE       ->      Divide by 0
   * IEEE_TRAP_ENABLE_OVF       ->      Overflow
   * IEEE_TRAP_ENABLE_UNF       ->      Underflow
   * IEEE_TRAP_ENABLE_INE       ->      Inexact (requires special option to C compiler)
   * IEEE_TRAP_ENABLE_DNO       ->      denormal operand
   * Note: IEEE_TRAP_ENABLE_DNO is not supported on OSF 3.x or Digital Unix
   * 4.0 - 4.0d(?).
   * IEEE_TRAP_ENABLE_MASK      ->      mask of all the trap enables
   * IEEE_MAP_DMZ                       ->      map denormal inputs to zero
   * IEEE_MAP_UMZ                       ->      map underflow results to zero
   */

  mode = IEEE_TRAP_ENABLE_INV | IEEE_TRAP_ENABLE_DZE | IEEE_TRAP_ENABLE_OVF
                | IEEE_TRAP_ENABLE_UNF ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ IEEE_TRAP_ENABLE_INV ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
#ifdef IEEE_TRAP_ENABLE_DNO
     mode &= ~ IEEE_TRAP_ENABLE_DNO ;
#else
     GSL_ERROR ("Sorry, this version of Digital Unix does not support denormalized operands", GSL_EUNSUP) ;
#endif
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ IEEE_TRAP_ENABLE_DZE ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ IEEE_TRAP_ENABLE_OVF ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ IEEE_TRAP_ENABLE_UNF ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      /* To implement this would require a special flag to the C
       compiler which can cause degraded performance */

      GSL_ERROR ("Sorry, GSL does not implement trap-inexact for Tru64 Unix on the alpha - see fp-tru64.c for details", GSL_EUNSUP) ;

      /* In case you need to add it, the appropriate line would be 
       *  
       *  mode |= IEEE_TRAP_ENABLE_INE ; 
       *
       */

    }
  else
    {
      mode &= ~ IEEE_TRAP_ENABLE_INE ;
    }

  ieee_set_fp_control (mode) ;

  return GSL_SUCCESS ;
}
#elif HAVE_FREEBSD_IEEE_INTERFACE
/* ieee-utils/fp-freebsd.c
 * 
 * Copyright (C) 2000 Vladimir Kushnir
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <ieeefp.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fp_prec_t prec = 0 ;
  fp_except_t mode = 0 ;
  fp_rnd_t    rnd  = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      prec = FP_PS;
      fpsetprec(prec);      
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      prec = FP_PD;
      fpsetprec(prec);
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      prec = FP_PE;
      fpsetprec(prec);
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RN ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RM ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RP ;
      fpsetround (rnd) ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RZ ;
      fpsetround (rnd) ;
      break ;
    default:
      rnd = FP_RN ;
      fpsetround (rnd) ;
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = FP_X_INV | FP_X_DNML | FP_X_DZ | FP_X_OFL | FP_X_UFL ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ FP_X_INV ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode &= ~ FP_X_DNML ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ FP_X_DZ ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ FP_X_OFL ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ FP_X_UFL ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= FP_X_IMP ;
    }
  else
    {
      mode &= ~ FP_X_IMP ;
    }

  fpsetmask (mode) ;

  return GSL_SUCCESS ;

}
#elif HAVE_OS2EMX_IEEE_INTERFACE
/* ieee-utils/fp-os2.c
 * 
 * Copyright (C) 2001 Henry Sobotka <sobotka@axess.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <float.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  unsigned mode = 0;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      _control87(PC_24, MCW_PC);    
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      _control87(PC_53, MCW_PC);
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      _control87(PC_64, MCW_PC);
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      _control87(RC_NEAR, MCW_RC);
      break ;
    case GSL_IEEE_ROUND_DOWN:
      _control87(RC_DOWN, MCW_RC);
      break ;
    case GSL_IEEE_ROUND_UP:
      _control87(RC_UP, MCW_RC);
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      _control87(RC_CHOP, MCW_RC);
      break ;
    default:
      _control87(RC_NEAR, MCW_RC);
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = EM_INVALID | EM_DENORMAL | EM_ZERODIVIDE | EM_OVERFLOW | EM_UNDERFLOW;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ EM_INVALID;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode &= ~ EM_DENORMAL;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ EM_ZERODIVIDE;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ EM_OVERFLOW;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &= ~ EM_UNDERFLOW;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= EM_INEXACT;
    }
  else
    {
      mode &= ~ EM_INEXACT;
    }

  _control87(mode, MCW_EM);

  return GSL_SUCCESS ;
}
#elif HAVE_NETBSD_IEEE_INTERFACE
/* fp-netbsd.c
 * 
 * Copyright (C) 2001 Jason Beegan
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <ieeefp.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fp_except mode = 0;
  fp_rnd    rnd  = 0;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("NetBSD only supports default precision rounding",
                 GSL_EUNSUP);
      break;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("NetBSD only supports default precision rounding",
                 GSL_EUNSUP);
      break;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("NetBSD only supports default precision rounding",
                 GSL_EUNSUP);
      break;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RN;
      fpsetround (rnd);
      break;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RM;
      fpsetround (rnd);
      break;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RP;
      fpsetround (rnd);
      break;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RZ;
      fpsetround (rnd);
      break;
    default:
      rnd = FP_RN;
      fpsetround (rnd);
    }

/* Turn on all available exceptions apart from 'inexact'.
   Denormalized operand exception not available on all platforms. */

  mode = FP_X_INV | FP_X_DZ | FP_X_OFL | FP_X_UFL;
#ifdef FP_X_DNML
  mode = mode | FP_X_DNML;
#endif

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ FP_X_INV;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
#ifdef FP_X_DNML
      mode &= ~ FP_X_DNML;
#endif
    }
  else
    {
#ifndef FP_X_DNML
      GSL_ERROR ("NetBSD does not support the denormalized operand exception on this platform. "
                 "Use 'mask-denormalized' to work around this.",
                 GSL_EUNSUP);
#endif
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ FP_X_DZ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ FP_X_OFL;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ FP_X_UFL;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= FP_X_IMP;
    }
  else
    {
      mode &= ~ FP_X_IMP;
    }

  fpsetmask (mode);

  return GSL_SUCCESS;

}

#elif HAVE_OPENBSD_IEEE_INTERFACE
/* fp-openbsd.c
 * 
 * Copyright (C) 2001 Jason Beegan
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <ieeefp.h>

/* This is a copy of fp-netbsd.c, modified for openbsd by Toby White
   --- Brian Gough */

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fp_except mode = 0;
  fp_rnd    rnd  = 0;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("OpenBSD only supports default precision rounding",
                 GSL_EUNSUP);
      break;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("OpenBSD only supports default precision rounding",
                 GSL_EUNSUP);
      break;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("OpenBSD only supports default precision rounding",
                 GSL_EUNSUP);
      break;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      rnd = FP_RN;
      fpsetround (rnd);
      break;
    case GSL_IEEE_ROUND_DOWN:
      rnd = FP_RM;
      fpsetround (rnd);
      break;
    case GSL_IEEE_ROUND_UP:
      rnd = FP_RP;
      fpsetround (rnd);
      break;
    case GSL_IEEE_ROUND_TO_ZERO:
      rnd = FP_RZ;
      fpsetround (rnd);
      break;
    default:
      rnd = FP_RN;
      fpsetround (rnd);
    }

/* Turn on all available exceptions apart from 'inexact'.
   Denormalized operand exception not available on all platforms. */

  mode = FP_X_INV | FP_X_DZ | FP_X_OFL | FP_X_UFL;
#ifdef FP_X_DNML
  mode = mode | FP_X_DNML;
#endif

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode &= ~ FP_X_INV;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
#ifdef FP_X_DNML
      mode &= ~ FP_X_DNML;
#endif
    }
  else
    {
#ifndef FP_X_DNML
      GSL_ERROR ("OpenBSD does not support the denormalized operand exception on this platform. "
                 "Use 'mask-denormalized' to work around this.",
                 GSL_EUNSUP);
#endif
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode &= ~ FP_X_DZ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode &= ~ FP_X_OFL;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode &=  ~ FP_X_UFL;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode |= FP_X_IMP;
    }
  else
    {
      mode &= ~ FP_X_IMP;
    }

  fpsetmask (mode);

  return GSL_SUCCESS;

}

/* Try to handle universal binaries */
#elif HAVE_DARWIN_IEEE_INTERFACE
# if defined(__i386__)
/* ieee-utils/fp-darwin86.c
 * 
 * Copyright (C) 2006 Erik Schnetter
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */

/* Here is the dirty part. Set up your 387 through the control word
 * (cw) register.
 *
 *     15-13    12  11-10  9-8     7-6     5    4    3    2    1    0
 * | reserved | IC | RC  | PC | reserved | PM | UM | OM | ZM | DM | IM
 *
 * IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask
 *
 * Mask bit is 1 means no interrupt.
 *
 * PC: Precision control
 * 11 - round to extended precision
 * 10 - round to double precision
 * 00 - round to single precision
 *
 * RC: Rounding control
 * 00 - rounding to nearest
 * 01 - rounding down (toward - infinity)
 * 10 - rounding up (toward + infinity)
 * 11 - rounding toward zero
 *
 * IC: Infinity control
 * That is for 8087 and 80287 only.
 *
 * The hardware default is 0x037f which we use.
 */

/* masking of interrupts */
#define _FPU_MASK_IM  0x01
#define _FPU_MASK_DM  0x02
#define _FPU_MASK_ZM  0x04
#define _FPU_MASK_OM  0x08
#define _FPU_MASK_UM  0x10
#define _FPU_MASK_PM  0x20

/* precision control */
#define _FPU_EXTENDED 0x300     /* libm requires double extended precision.  */
#define _FPU_DOUBLE   0x200
#define _FPU_SINGLE   0x0

/* rounding control */
#define _FPU_RC_NEAREST 0x0    /* RECOMMENDED */
#define _FPU_RC_DOWN    0x400
#define _FPU_RC_UP      0x800
#define _FPU_RC_ZERO    0xC00

#define _FPU_RESERVED 0xF0C0  /* Reserved bits in cw */


/* The fdlibm code requires strict IEEE double precision arithmetic,
   and no interrupts for exceptions, rounding to nearest.  */

#define _FPU_DEFAULT  0x037f

/* IEEE:  same as above.  */
#define _FPU_IEEE     0x037f

/* Type of the control word.  */
typedef unsigned int fpu_control_t __attribute__ ((__mode__ (__HI__)));

/* Macros for accessing the hardware control word.

   Note that the use of these macros is no sufficient anymore with
   recent hardware.  Some floating point operations are executed in
   the SSE/SSE2 engines which have their own control and status register.  */
#define _FPU_GETCW(cw) __asm__ __volatile__ ("fnstcw %0" : "=m" (*&cw))
#define _FPU_SETCW(cw) __asm__ __volatile__ ("fldcw %0" : : "m" (*&cw))

/* Default control word set at startup.  */
extern fpu_control_t __fpu_control;



#define _FPU_GETMXCSR(cw_sse) asm volatile ("stmxcsr %0" : "=m" (cw_sse))
#define _FPU_SETMXCSR(cw_sse) asm volatile ("ldmxcsr %0" : : "m" (cw_sse))



int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fpu_control_t mode, mode_sse;

  _FPU_GETCW (mode) ;
  mode &= _FPU_RESERVED ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      mode |= _FPU_SINGLE ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      mode |= _FPU_DOUBLE ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      mode |= _FPU_EXTENDED ;
      break ;
    default:
      mode |= _FPU_EXTENDED ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case GSL_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_IM ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode |= _FPU_MASK_DM ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_ZM ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OM ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UM ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode &= ~ _FPU_MASK_PM ;
    }
  else
    {
      mode |= _FPU_MASK_PM ;
    }

  _FPU_SETCW (mode) ;

  _FPU_GETMXCSR (mode_sse) ;
  mode_sse &= 0xFFFF0000 ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode_sse |= _FPU_MASK_IM << 7 ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode_sse |= _FPU_MASK_DM << 7 ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode_sse |= _FPU_MASK_ZM << 7 ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode_sse |= _FPU_MASK_OM << 7 ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode_sse |= _FPU_MASK_UM << 7 ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode_sse &= ~ _FPU_MASK_PM << 7 ;
    }
  else
    {
      mode_sse |= _FPU_MASK_PM << 7 ;
    }

  _FPU_SETMXCSR (mode_sse) ;

  return GSL_SUCCESS ;
}
#else
/* ieee-utils/fp-darwin.c
 * 
 * Copyright (C) 2001 Rodney Sparapani <rsparapa@mcw.edu>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <architecture/ppc/fp_regs.h> 

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  ppc_fp_scr_t fp_scr = get_fp_scr() ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      fp_scr.rn = RN_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      fp_scr.rn = RN_TOWARD_MINUS ;
      break ;
    case GSL_IEEE_ROUND_UP:
      fp_scr.rn = RN_TOWARD_PLUS ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      fp_scr.rn = RN_TOWARD_ZERO ;
      break ;
    default:
      fp_scr.rn = RN_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    fp_scr.ve = 0 ;                             //ve bit:  invalid op exception enable

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("powerpc does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.", GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    fp_scr.ze = 0 ;                             //ze bit:  zero divide exception enable

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    fp_scr.oe = 0 ;                             //oe bit:  overflow exception enable

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    fp_scr.ue  = 0 ;                            //ue bit:  underflow exception enable

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      fp_scr.xe = 1 ;                           //xe bit:  inexact exception enable
    }
  else
    {
      fp_scr.xe = 01 ;                  
    }

  set_fp_scr(fp_scr);

  return GSL_SUCCESS ;

}
# endif
#elif HAVE_DARWIN86_IEEE_INTERFACE
# if defined(__ppc__)
/* ieee-utils/fp-darwin.c
 * 
 * Copyright (C) 2001 Rodney Sparapani <rsparapa@mcw.edu>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <architecture/ppc/fp_regs.h> 

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  ppc_fp_scr_t fp_scr = get_fp_scr() ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("powerpc only supports default precision rounding", GSL_EUNSUP);
      break ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      fp_scr.rn = RN_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      fp_scr.rn = RN_TOWARD_MINUS ;
      break ;
    case GSL_IEEE_ROUND_UP:
      fp_scr.rn = RN_TOWARD_PLUS ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      fp_scr.rn = RN_TOWARD_ZERO ;
      break ;
    default:
      fp_scr.rn = RN_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    fp_scr.ve = 0 ;                             //ve bit:  invalid op exception enable

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("powerpc does not support the denormalized operand exception. "
                 "Use 'mask-denormalized' to work around this.", GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    fp_scr.ze = 0 ;                             //ze bit:  zero divide exception enable

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    fp_scr.oe = 0 ;                             //oe bit:  overflow exception enable

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    fp_scr.ue  = 0 ;                            //ue bit:  underflow exception enable

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      fp_scr.xe = 1 ;                           //xe bit:  inexact exception enable
    }
  else
    {
      fp_scr.xe = 01 ;                  
    }

  set_fp_scr(fp_scr);

  return GSL_SUCCESS ;

}
# else
/* ieee-utils/fp-darwin86.c
 * 
 * Copyright (C) 2006 Erik Schnetter
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */

/* Here is the dirty part. Set up your 387 through the control word
 * (cw) register.
 *
 *     15-13    12  11-10  9-8     7-6     5    4    3    2    1    0
 * | reserved | IC | RC  | PC | reserved | PM | UM | OM | ZM | DM | IM
 *
 * IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask
 *
 * Mask bit is 1 means no interrupt.
 *
 * PC: Precision control
 * 11 - round to extended precision
 * 10 - round to double precision
 * 00 - round to single precision
 *
 * RC: Rounding control
 * 00 - rounding to nearest
 * 01 - rounding down (toward - infinity)
 * 10 - rounding up (toward + infinity)
 * 11 - rounding toward zero
 *
 * IC: Infinity control
 * That is for 8087 and 80287 only.
 *
 * The hardware default is 0x037f which we use.
 */

/* masking of interrupts */
#define _FPU_MASK_IM  0x01
#define _FPU_MASK_DM  0x02
#define _FPU_MASK_ZM  0x04
#define _FPU_MASK_OM  0x08
#define _FPU_MASK_UM  0x10
#define _FPU_MASK_PM  0x20

/* precision control */
#define _FPU_EXTENDED 0x300     /* libm requires double extended precision.  */
#define _FPU_DOUBLE   0x200
#define _FPU_SINGLE   0x0

/* rounding control */
#define _FPU_RC_NEAREST 0x0    /* RECOMMENDED */
#define _FPU_RC_DOWN    0x400
#define _FPU_RC_UP      0x800
#define _FPU_RC_ZERO    0xC00

#define _FPU_RESERVED 0xF0C0  /* Reserved bits in cw */


/* The fdlibm code requires strict IEEE double precision arithmetic,
   and no interrupts for exceptions, rounding to nearest.  */

#define _FPU_DEFAULT  0x037f

/* IEEE:  same as above.  */
#define _FPU_IEEE     0x037f

/* Type of the control word.  */
typedef unsigned int fpu_control_t __attribute__ ((__mode__ (__HI__)));

/* Macros for accessing the hardware control word.

   Note that the use of these macros is no sufficient anymore with
   recent hardware.  Some floating point operations are executed in
   the SSE/SSE2 engines which have their own control and status register.  */
#define _FPU_GETCW(cw) __asm__ __volatile__ ("fnstcw %0" : "=m" (*&cw))
#define _FPU_SETCW(cw) __asm__ __volatile__ ("fldcw %0" : : "m" (*&cw))

/* Default control word set at startup.  */
extern fpu_control_t __fpu_control;



#define _FPU_GETMXCSR(cw_sse) asm volatile ("stmxcsr %0" : "=m" (cw_sse))
#define _FPU_SETMXCSR(cw_sse) asm volatile ("ldmxcsr %0" : : "m" (cw_sse))



int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  fpu_control_t mode, mode_sse;

  _FPU_GETCW (mode) ;
  mode &= _FPU_RESERVED ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      mode |= _FPU_SINGLE ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      mode |= _FPU_DOUBLE ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      mode |= _FPU_EXTENDED ;
      break ;
    default:
      mode |= _FPU_EXTENDED ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case GSL_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_IM ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode |= _FPU_MASK_DM ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_ZM ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OM ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UM ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode &= ~ _FPU_MASK_PM ;
    }
  else
    {
      mode |= _FPU_MASK_PM ;
    }

  _FPU_SETCW (mode) ;

  _FPU_GETMXCSR (mode_sse) ;
  mode_sse &= 0xFFFF0000 ;

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode_sse |= _FPU_MASK_IM << 7 ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode_sse |= _FPU_MASK_DM << 7 ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode_sse |= _FPU_MASK_ZM << 7 ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode_sse |= _FPU_MASK_OM << 7 ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode_sse |= _FPU_MASK_UM << 7 ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode_sse &= ~ _FPU_MASK_PM << 7 ;
    }
  else
    {
      mode_sse |= _FPU_MASK_PM << 7 ;
    }

  _FPU_SETMXCSR (mode_sse) ;

  return GSL_SUCCESS ;
}
#endif
#elif HAVE_DECL_FEENABLEEXCEPT || HAVE_DECL_FESETTRAPENABLE
/* ieee-utils/fp-gnuc99.c
 * 
 * Copyright (C) 2003, 2004, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#define _GNU_SOURCE 1

#include <math.h>
#include <stdio.h>
#include <fenv.h>

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  int mode;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      GSL_ERROR ("single precision rounding is not supported by <fenv.h>",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      GSL_ERROR ("double precision rounding is not supported by <fenv.h>",
                 GSL_EUNSUP) ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      GSL_ERROR ("extended precision rounding is not supported by <fenv.h>",
                 GSL_EUNSUP) ;
      break ;
    }


  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
#ifdef FE_TONEAREST
      fesetround (FE_TONEAREST) ;
#else
      GSL_ERROR ("round-to-nearest is not supported by <fenv.h>", GSL_EUNSUP) ;
#endif
      break ;
    case GSL_IEEE_ROUND_DOWN:
#ifdef FE_DOWNWARD
      fesetround (FE_DOWNWARD) ;
#else
      GSL_ERROR ("round-down is not supported by <fenv.h>", GSL_EUNSUP) ;
#endif
      break ;
    case GSL_IEEE_ROUND_UP:
#ifdef FE_UPWARD
      fesetround (FE_UPWARD) ;
#else
      GSL_ERROR ("round-up is not supported by <fenv.h>", GSL_EUNSUP) ;
#endif
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
#ifdef FE_TOWARDZERO
      fesetround (FE_TOWARDZERO) ;
#else
      GSL_ERROR ("round-toward-zero is not supported by <fenv.h>", GSL_EUNSUP) ;
#endif
      break ;
    default:
#ifdef FE_TONEAREST
      fesetround (FE_TONEAREST) ;
#else
      GSL_ERROR ("default round-to-nearest mode is not supported by <fenv.h>", GSL_EUNSUP) ;
#endif
    }

  /* Turn on all the exceptions apart from 'inexact' */

  mode = 0;

#ifdef FE_INVALID 
  mode |= FE_INVALID;
#endif

#ifdef FE_DIVBYZERO
  mode |= FE_DIVBYZERO;
#endif
  
#ifdef FE_OVERFLOW
  mode |= FE_OVERFLOW ;
#endif

#ifdef FE_UNDERFLOW
  mode |= FE_UNDERFLOW ;
#endif

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    {
#ifdef FE_INVALID
    mode &= ~ FE_INVALID ;
#else
    GSL_ERROR ("invalid operation exception not supported by <fenv.h>", 
               GSL_EUNSUP);
#endif
    }

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    {
      /* do nothing */
    }
  else
    {
      GSL_ERROR ("denormalized operand exception not supported by <fenv.h>. "
                 "Use 'mask-denormalized' to work around this.", GSL_EUNSUP) ;
    }

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    {
#ifdef FE_DIVBYZERO
      mode &= ~ FE_DIVBYZERO ;
#else
      GSL_ERROR ("division by zero exception not supported by <fenv.h>", 
                 GSL_EUNSUP);
#endif
    }

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    {
#ifdef FE_OVERFLOW
      mode &= ~ FE_OVERFLOW ;
#else
      GSL_ERROR ("overflow exception not supported by <fenv.h>", GSL_EUNSUP);
#endif
    }

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    {
#ifdef FE_UNDERFLOW
      mode &=  ~ FE_UNDERFLOW ;
#else
      GSL_ERROR ("underflow exception not supported by <fenv.h>", GSL_EUNSUP);
#endif
    }

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
#ifdef FE_INEXACT
      mode |= FE_INEXACT ;
#else
      GSL_ERROR ("inexact exception not supported by <fenv.h>", GSL_EUNSUP);
#endif
    }
  else
    {
#ifdef FE_INEXACT
      mode &= ~ FE_INEXACT ;
#else
      /* do nothing */
#endif
    }

#if HAVE_DECL_FEENABLEEXCEPT
  feenableexcept (mode) ;
#elif HAVE_DECL_FESETTRAPENABLE
  fesettrapenable (mode);
#else
  GSL_ERROR ("unknown exception trap method", GSL_EUNSUP)
#endif

  return GSL_SUCCESS ;
}
#else
/* ieee-utils/fp-unknown.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */

int
gsl_ieee_set_mode (int precision, int rounding, int exception_mask)
{
  GSL_ERROR (
"the IEEE interface for this platform is unsupported or could not be "
"determined at configure time\n", GSL_EUNSUP) ;
}
#endif
/* ieee-utils/make_rep.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */

/* ieee-utils/endian.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <stdlib.h>

static int little_endian_p (void) ;

static int 
little_endian_p (void) {
  /* Are we little or big endian?  From Harbison & Steele.  */
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}

/* ieee-utils/standardize.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

static void make_float_bigendian (float * x);
static void make_double_bigendian (double * x);

static void
make_float_bigendian (float * x)
{
  union { 
    float f;
    unsigned char b[4];
  } u,v;

  u.f = *x ;

  v.b[0]=u.b[3] ;
  v.b[1]=u.b[2] ;
  v.b[2]=u.b[1] ;
  v.b[3]=u.b[0] ;

  *x=v.f ;
}

static void
make_double_bigendian (double * x)
{
  union { 
    double d;
    unsigned char b[8];
  } u,v;

  u.d = *x ;

  v.b[0]=u.b[7] ;
  v.b[1]=u.b[6] ;
  v.b[2]=u.b[5] ;
  v.b[3]=u.b[4] ;
  v.b[4]=u.b[3] ;
  v.b[5]=u.b[2] ;
  v.b[6]=u.b[1] ;
  v.b[7]=u.b[0] ;

  *x=v.d ;
}

static void sprint_nybble(int i, char *s) ;
static void sprint_byte(int i, char *s) ;
static int determine_ieee_type (int non_zero, int exponent, int max_exponent);


/* For the IEEE float format the bits are found from the following
   masks,
   
   sign      = 0x80000000  
   exponent  = 0x7f800000 
   mantisssa = 0x007fffff  

   For the IEEE double format the masks are,

   sign      = 0x8000000000000000  
   exponent  = 0x7ff0000000000000 
   mantissa  = 0x000fffffffffffff

   */

void 
gsl_ieee_float_to_rep (const float * x, gsl_ieee_float_rep * r)
{
  int e, non_zero;

  union { 
    float f;
    struct  { 
      unsigned char byte[4] ;
    } ieee ;
  } u;
  
  u.f = *x ; 

  if (little_endian_p())
    make_float_bigendian(&(u.f)) ;
  
  /* note that r->sign is signed, u.ieee.byte is unsigned */

  if (u.ieee.byte[3]>>7)
    {
      r->sign = 1 ;
    }
  else
    {
      r->sign = 0 ;
    }

  e = (u.ieee.byte[3] & 0x7f) << 1 | (u.ieee.byte[2] & 0x80)>>7 ; 
  
  r->exponent = e - 127 ;

  sprint_byte((u.ieee.byte[2] & 0x7f) << 1,r->mantissa) ;
  sprint_byte(u.ieee.byte[1],r->mantissa + 7) ;
  sprint_byte(u.ieee.byte[0],r->mantissa + 15) ;

  r->mantissa[23] = '\0' ;

  non_zero = u.ieee.byte[0] || u.ieee.byte[1] || (u.ieee.byte[2] & 0x7f);

  r->type = determine_ieee_type (non_zero, e, 255) ;
}

void 
gsl_ieee_double_to_rep (const double * x, gsl_ieee_double_rep * r)
{

  int e, non_zero;

  union 
  { 
    double d;
    struct  { 
      unsigned char byte[8];
    } ieee ;
  } u;

  u.d= *x ; 
  
  if (little_endian_p())
    make_double_bigendian(&(u.d)) ;
  
  /* note that r->sign is signed, u.ieee.byte is unsigned */

  if (u.ieee.byte[7]>>7)
    {
      r->sign = 1 ;
    }
  else
    {
      r->sign = 0 ;
    }


  e =(u.ieee.byte[7] & 0x7f)<<4 ^ (u.ieee.byte[6] & 0xf0)>>4 ;
  
  r->exponent = e - 1023 ;

  sprint_nybble(u.ieee.byte[6],r->mantissa) ;
  sprint_byte(u.ieee.byte[5],r->mantissa + 4) ;
  sprint_byte(u.ieee.byte[4],r->mantissa + 12) ;
  sprint_byte(u.ieee.byte[3],r->mantissa + 20) ; 
  sprint_byte(u.ieee.byte[2],r->mantissa + 28) ;
  sprint_byte(u.ieee.byte[1],r->mantissa + 36) ;
  sprint_byte(u.ieee.byte[0],r->mantissa + 44) ;

  r->mantissa[52] = '\0' ;

  non_zero = (u.ieee.byte[0] || u.ieee.byte[1] || u.ieee.byte[2]
              || u.ieee.byte[3] || u.ieee.byte[4] || u.ieee.byte[5] 
              || (u.ieee.byte[6] & 0x0f)) ;

  r->type = determine_ieee_type (non_zero, e, 2047) ;
}

/* A table of character representations of nybbles */

static char nybble[16][5]={ /* include space for the \0 */
  "0000", "0001", "0010", "0011",
  "0100", "0101", "0110", "0111",
  "1000", "1001", "1010", "1011",
  "1100", "1101", "1110", "1111"
}  ;
          
static void
sprint_nybble(int i, char *s)
{
  char *c ;
  c=nybble[i & 0x0f ];
  *s=c[0] ;  *(s+1)=c[1] ;  *(s+2)=c[2] ;  *(s+3)=c[3] ;
} 

static void
sprint_byte(int i, char *s)
{
  char *c ;
  c=nybble[(i & 0xf0)>>4];
  *s=c[0] ;  *(s+1)=c[1] ;  *(s+2)=c[2] ;  *(s+3)=c[3] ;
  c=nybble[i & 0x0f];
  *(s+4)=c[0] ;  *(s+5)=c[1] ;  *(s+6)=c[2] ;  *(s+7)=c[3] ;
} 

static int 
determine_ieee_type (int non_zero, int exponent, int max_exponent)
{
  if (exponent == max_exponent)
    {
      if (non_zero)
        {
          return GSL_IEEE_TYPE_NAN ;
        }
      else
        {
          return GSL_IEEE_TYPE_INF ;
        }
    }
  else if (exponent == 0)
    {
      if (non_zero)
        {
          return GSL_IEEE_TYPE_DENORMAL ;
        }
      else
        {
          return GSL_IEEE_TYPE_ZERO ;
        }
    }
  else
    {
      return GSL_IEEE_TYPE_NORMAL ;
    }
}
/* ieee-utils/print.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <stdio.h>
#include <math.h>

/* A table of sign characters, 0=positive, 1=negative. We print a space
   instead of a unary + sign for compatibility with bc */

static char signs[2]={' ','-'} ;  

void 
gsl_ieee_fprintf_float (FILE * stream, const float * x) {
  gsl_ieee_float_rep r ;
  gsl_ieee_float_to_rep(x, &r) ;

  switch (r.type)
    {
    case GSL_IEEE_TYPE_NAN:
      fprintf(stream, "NaN") ;
      break ;
    case GSL_IEEE_TYPE_INF:
      fprintf(stream, "%cInf", signs[r.sign]) ;
      break ;
    case GSL_IEEE_TYPE_NORMAL:
      fprintf(stream, "%c1.%s*2^%d", signs[r.sign], r.mantissa, r.exponent) ;
      break ;
    case GSL_IEEE_TYPE_DENORMAL:
      fprintf(stream, "%c0.%s*2^%d", signs[r.sign], r.mantissa, r.exponent + 1) ;
      break ;
    case GSL_IEEE_TYPE_ZERO:
      fprintf(stream, "%c0", signs[r.sign]) ;
      break ;
    default:
      fprintf(stream, "[non-standard IEEE float]") ;
    }
}

void 
gsl_ieee_printf_float (const float * x)
{
  gsl_ieee_fprintf_float (stdout,x);
}

void
gsl_ieee_fprintf_double (FILE * stream, const double * x) {
  gsl_ieee_double_rep r ;
  gsl_ieee_double_to_rep (x, &r) ;

  switch (r.type)
    {
    case GSL_IEEE_TYPE_NAN:
      fprintf(stream, "NaN") ;
      break ;
    case GSL_IEEE_TYPE_INF:
      fprintf(stream, "%cInf", signs[r.sign]) ;
      break ;
    case GSL_IEEE_TYPE_NORMAL:
      fprintf(stream, "%c1.%s*2^%d", signs[r.sign], r.mantissa, r.exponent) ;
      break ;
    case GSL_IEEE_TYPE_DENORMAL:
      fprintf(stream, "%c0.%s*2^%d", signs[r.sign], r.mantissa, r.exponent + 1) ;
      break ;
    case GSL_IEEE_TYPE_ZERO:
      fprintf(stream, "%c0", signs[r.sign]) ;
      break ;
    default:
      fprintf(stream, "[non-standard IEEE double]") ;
    }
}

void 
gsl_ieee_printf_double (const double * x)
{
  gsl_ieee_fprintf_double (stdout,x);
}





/* ieee-utils/read.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <stdlib.h>
#include <string.h>

static int 
lookup_string (const char * p, int * precision, int * rounding, 
               int * exception_mask) ;

int
gsl_ieee_read_mode_string (const char * description, 
                           int * precision, 
                           int * rounding, 
                           int * exception_mask)
{
  char * start ;
  char * end;
  char * p;

  int precision_count = 0 ;
  int rounding_count = 0 ;
  int exception_count = 0 ;

  start = (char *) malloc(strlen(description) + 1) ;

  if (start == 0) 
    {
      GSL_ERROR ("no memory to parse mode string", GSL_ENOMEM) ;
    }

  strcpy (start, description) ;

  p = start ;

  *precision = 0 ;
  *rounding = 0 ;
  *exception_mask = 0 ;

  do {
    int status ;
    int new_precision, new_rounding, new_exception ;

    end = strchr (p,',') ;

    if (end) 
      {
        *end = '\0' ;
        do 
          {
            end++ ;  /* skip over trailing whitespace */
          } 
        while (*end == ' ' || *end == ',') ;
      }
        
    new_precision = 0 ; 
    new_rounding = 0 ; 
    new_exception = 0 ;

    status = lookup_string (p, &new_precision, &new_rounding, &new_exception) ;

    if (status)
      GSL_ERROR ("unrecognized GSL_IEEE_MODE string.\nValid settings are:\n\n" 
                 "  single-precision double-precision extended-precision\n"
                 "  round-to-nearest round-down round-up round-to-zero\n"
                 "  mask-invalid mask-denormalized mask-division-by-zero\n"
                 "  mask-overflow mask-underflow mask-all\n"
                 "  trap-common trap-inexact\n"
                 "\n"
                 "separated by commas. "
                 "(e.g. GSL_IEEE_MODE=\"round-down,mask-underflow\")",
                 GSL_EINVAL) ;

    if (new_precision) 
      {
        *precision = new_precision ;
        precision_count ++ ;
        if (precision_count > 1)
          GSL_ERROR ("attempted to set IEEE precision twice", GSL_EINVAL) ;
      }

    if (new_rounding) 
      {
        *rounding = new_rounding ;
        rounding_count ++ ;
        if (rounding_count > 1)
          GSL_ERROR ("attempted to set IEEE rounding mode twice", GSL_EINVAL) ;
      }

    if (new_exception) 
      {
        *exception_mask |= new_exception ;
        exception_count ++ ;
      }

    p = end ; 

  } while (end && *p != '\0') ;

  free(start) ;

  return GSL_SUCCESS ;
}

static int 
lookup_string (const char * p, int * precision, int * rounding, 
               int * exception_mask)
{
  if (strcmp(p,"single-precision") == 0) 
    {
      *precision = GSL_IEEE_SINGLE_PRECISION ;
    }
  else if (strcmp(p,"double-precision") == 0) 
    {
      *precision = GSL_IEEE_DOUBLE_PRECISION ;
    }
  else if (strcmp(p,"extended-precision") == 0) 
    {
      *precision = GSL_IEEE_EXTENDED_PRECISION ;
    }
  else if (strcmp(p,"round-to-nearest") == 0) 
    {
      *rounding = GSL_IEEE_ROUND_TO_NEAREST ;
    }
  else if (strcmp(p,"round-down") == 0) 
    {
      *rounding = GSL_IEEE_ROUND_DOWN ;
    }
  else if (strcmp(p,"round-up") == 0) 
    {
      *rounding = GSL_IEEE_ROUND_UP ;
    }
  else if (strcmp(p,"round-to-zero") == 0) 
    {
      *rounding = GSL_IEEE_ROUND_TO_ZERO ;
    }
  else if (strcmp(p,"mask-all") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_ALL ;
    }
  else if (strcmp(p,"mask-invalid") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_INVALID ;
    }
  else if (strcmp(p,"mask-denormalized") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_DENORMALIZED ;
    }
  else if (strcmp(p,"mask-division-by-zero") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_DIVISION_BY_ZERO ;
    }
  else if (strcmp(p,"mask-overflow") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_OVERFLOW ;
    }
  else if (strcmp(p,"mask-underflow") == 0) 
    {
      *exception_mask = GSL_IEEE_MASK_UNDERFLOW ;
    }
  else if (strcmp(p,"trap-inexact") == 0) 
    {
      *exception_mask = GSL_IEEE_TRAP_INEXACT ;
    }
  else if (strcmp(p,"trap-common") == 0) 
    {
      return 0 ;
    }
  else
    {
      return 1 ;
    }

  return 0 ;
}
/* version.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
/* #include <gsl_version.h> */

/* This file needs to use the top-level <gsl_version.h> due to the
   possibility of a VPATH-style build where the original source
   tree is on read-only filesystem and so will not be picked up
   by the symlinking comands in gsl/Makefile.am */


const char * gsl_version = GSL_VERSION;

/* sys/coerce.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>

double 
gsl_coerce_double (const double x)
{
  volatile double y;
  y = x;
  return y;
}

float 
gsl_coerce_float (const float x)
{
  volatile float y;
  y = x;
  return y;
}

/* The following function is not needed, but is included for completeness */

long double 
gsl_coerce_long_double (const long double x)
{
  volatile long double y;
  y = x;
  return y;
}

/* sys/expm1.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>

double gsl_expm1 (const double x)
{
  /* FIXME: this should be improved */

  if (fabs(x) < M_LN2)
    {
      /* Compute the taylor series S = x + (1/2!) x^2 + (1/3!) x^3 + ... */

      double i = 1.0;
      double sum = x;
      double term = x / 1.0;

      do
        {
          i++ ;
          term *= x/i;
          sum += term;
        }
      while (fabs(term) > fabs(sum) * GSL_DBL_EPSILON) ;
      
      return sum ;
    }
  else
    {
      return exp(x) - 1;
    }
}
/* sys/gsl_compare.c
 * 
 * Copyright (C) 2002 Gert Van den Eynde
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 * 
 * Based on fcmp 1.2.2 Copyright (c) 1998-2000 Theodore C. Belding
 * University of Michigan Center for the Study of Complex Systems
 * Ted.Belding@umich.edu
 *
 */

/* #include <config.h> */
#include <math.h>

int
gsl_fcmp (const double x1, const double x2, const double epsilon)
{
  int exponent;
  double delta, difference;

  /* Find exponent of largest absolute value */

  {
    double max = (fabs (x1) > fabs (x2)) ? x1 : x2;

    frexp (max, &exponent);
  }

  /* Form a neighborhood of size  2 * delta */

  delta = ldexp (epsilon, exponent);

  difference = x1 - x2;

  if (difference > delta)       /* x1 > x2 */
    {
      return 1;
    }
  else if (difference < -delta) /* x1 < x2 */
    {
      return -1;
    }
  else                          /* -delta <= difference <= delta */
    {
      return 0;                 /* x1 ~=~ x2 */
    }
}

/* sys/hypot.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough, Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>

double gsl_hypot (const double x, const double y)
{
  double xabs = fabs(x) ;
  double yabs = fabs(y) ;
  double min, max;

  if (xabs < yabs) {
    min = xabs ;
    max = yabs ;
  } else {
    min = yabs ;
    max = xabs ;
  }

  if (min == 0) 
    {
      return max ;
    }

  {
    double u = min / max ;
    return max * sqrt (1 + u * u) ;
  }
}

double
gsl_hypot3(const double x, const double y, const double z)
{
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);
  double w = GSL_MAX(xabs, GSL_MAX(yabs, zabs));

  if (w == 0.0)
    {
      return (0.0);
    }
  else
    {
      double r = w * sqrt((xabs / w) * (xabs / w) +
                          (yabs / w) * (yabs / w) +
                          (zabs / w) * (zabs / w));
      return r;
    }
}
/* sys/invhyp.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>

double
gsl_acosh (const double x)
{
  if (x > 1.0 / GSL_SQRT_DBL_EPSILON)
    {
      return log (x) + M_LN2;
    }
  else if (x > 2)
    {
      return log (2 * x - 1 / (sqrt (x * x - 1) + x));
    }
  else if (x > 1)
    {
      double t = x - 1;
      return log1p (t + sqrt (2 * t + t * t));
    }
  else if (x == 1)
    {
      return 0;
    }
  else
    {
      return GSL_NAN;
    }
}

double
gsl_asinh (const double x)
{
  double a = fabs (x);
  double s = (x < 0) ? -1 : 1;

  if (a > 1 / GSL_SQRT_DBL_EPSILON)
    {
      return s * (log (a) + M_LN2);
    }
  else if (a > 2)
    {
      return s * log (2 * a + 1 / (a + sqrt (a * a + 1)));
    }
  else if (a > GSL_SQRT_DBL_EPSILON)
    {
      double a2 = a * a;
      return s * log1p (a + a2 / (1 + sqrt (1 + a2)));
    }
  else
    {
      return x;
    }
}

double
gsl_atanh (const double x)
{
  double a = fabs (x);
  double s = (x < 0) ? -1 : 1;

  if (a > 1)
    {
      return GSL_NAN;
    }
  else if (a == 1)
    {
      return (x < 0) ? GSL_NEGINF : GSL_POSINF;
    }
  else if (a >= 0.5)
    {
      return s * 0.5 * log1p (2 * a / (1 - a));
    }
  else if (a > GSL_DBL_EPSILON)
    {
      return s * 0.5 * log1p (2 * a + 2 * a * a / (1 - a));
    }
  else
    {
      return x;
    }
}
/* sys/ldfrexp.c
 * 
 * Copyright (C) 2002, Gert Van den Eynde
 * Copyright (C) 2007, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>

double
gsl_ldexp (const double x, const int e)
{
  int ex;
  
  if (x == 0.0)
    {
      return x;
    }

  {
    double y = gsl_frexp (x, &ex);
    double e2 = e + ex, p2;
    
    if (e2 >= DBL_MAX_EXP)
      {
        y *= pow (2.0, e2 - DBL_MAX_EXP + 1);
        e2 = DBL_MAX_EXP - 1;
      }
    else if (e2 <= DBL_MIN_EXP)
      {
        y *= pow (2.0, e2 - DBL_MIN_EXP - 1);
        e2 = DBL_MIN_EXP + 1;
      }
    
    p2 = pow (2.0, e2);
    return y * p2;
  }
}

double
gsl_frexp (const double x, int *e)
{
  if (x == 0.0)
    {
      *e = 0;
      return 0.0;
    }
  else if (!finite (x))
    {
      *e = 0;
      return x;
    }
  else if (fabs (x) >= 0.5 && fabs (x) < 1)     /* Handle the common case */
    {
      *e = 0;
      return x;
    }
  else
    {
      double ex = ceil (log (fabs (x)) / M_LN2);
      int ei = (int) ex;
      double f;

      /* Prevent underflow and overflow of 2**(-ei) */
      if (ei < DBL_MIN_EXP)
        ei = DBL_MIN_EXP;

      if (ei > -DBL_MIN_EXP)
        ei = -DBL_MIN_EXP;

      f = x * pow (2.0, -ei);

      if (!finite (f))
        {
          /* This should not happen */
          *e = 0;
          return f;
        }

      while (fabs (f) >= 1.0)
        {
          ei++;
          f /= 2.0;
        }

      while (fabs (f) > 0 && fabs (f) < 0.5)
        {
          ei--;
          f *= 2.0;
        }

      *e = ei;
      return f;
    }
}
/* sys/log1p.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */
#include <math.h>

double gsl_log1p (const double x)
{
  volatile double y, z;
  y = 1 + x;
  z = y - 1;
  return log(y) - (z-x)/y ;  /* cancels errors with IEEE arithmetic */
}
/* sys/minmax.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* #include <config.h> */

/* Compile all the inline functions */

#define COMPILE_INLINE_STATIC
/* Compile subsequent inline functions as static functions */

#ifdef __GSL_BUILD_H__
#error build.h must not be included multiple times
#endif

#define __GSL_BUILD_H__

#ifdef COMPILE_INLINE_STATIC
#ifndef HIDE_INLINE_STATIC  /* skip if inline functions are hidden */

#undef __GSL_INLINE_H__
#define __GSL_INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif     

/* Compile range checking code for inline functions used in the library */
#undef GSL_RANGE_CHECK
#define GSL_RANGE_CHECK 1

/* Use the global variable gsl_check_range to enable/disable range checking at
   runtime */
#undef GSL_RANGE_COND
#define GSL_RANGE_COND(x) (gsl_check_range && (x))

#endif
#else 
#error must be called with COMPILE_INLINE_STATIC
#endif

/* imagsl gsl minmax h: base directory gsl_minmax.h is added to imagsl.h */

/* Define some static functions which are always available */

double gsl_max (double a, double b)
{
  return GSL_MAX (a, b);
}

double gsl_min (double a, double b)
{
  return GSL_MIN (a, b);
}
#undef __GSL_BUILD_H__

/* sys/pow_int.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
/* #include <config.h> */
#include <math.h>

/* Compile all the inline functions */

#define COMPILE_INLINE_STATIC
/* Compile subsequent inline functions as static functions */

#ifdef __GSL_BUILD_H__
#error build.h must not be included multiple times
#endif

#define __GSL_BUILD_H__

#ifdef COMPILE_INLINE_STATIC
#ifndef HIDE_INLINE_STATIC  /* skip if inline functions are hidden */

#undef __GSL_INLINE_H__
#define __GSL_INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif     

/* Compile range checking code for inline functions used in the library */
#undef GSL_RANGE_CHECK
#define GSL_RANGE_CHECK 1

/* Use the global variable gsl_check_range to enable/disable range checking at
   runtime */
#undef GSL_RANGE_COND
#define GSL_RANGE_COND(x) (gsl_check_range && (x))

#endif
#else 
#error must be called with COMPILE_INLINE_STATIC
#endif

/* imagsl gsl pow int h */

double gsl_pow_int(double x, int n)
{
  double value = 1.0;

  if(n < 0) {
    x = 1.0/x;
    n = -n;
  }

  /* repeated squaring method 
   * returns 0.0^0 = 1.0, so continuous in x
   */
  do {
     if(n & 1) value *= x;  /* for n odd */
     n >>= 1;
     x *= x;
  } while (n);

  return value;
}
#undef __GSL_BUILD_H__

/* sys/prec.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

/* #include <config.h> */
/* imagsl gsl machine h */

/* Compile all the inline functions */

#define COMPILE_INLINE_STATIC
/* Compile subsequent inline functions as static functions */

#ifdef __GSL_BUILD_H__
#error build.h must not be included multiple times
#endif

#define __GSL_BUILD_H__

#ifdef COMPILE_INLINE_STATIC
#ifndef HIDE_INLINE_STATIC  /* skip if inline functions are hidden */

#undef __GSL_INLINE_H__
#define __GSL_INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif     

/* Compile range checking code for inline functions used in the library */
#undef GSL_RANGE_CHECK
#define GSL_RANGE_CHECK 1

/* Use the global variable gsl_check_range to enable/disable range checking at
   runtime */
#undef GSL_RANGE_COND
#define GSL_RANGE_COND(x) (gsl_check_range && (x))

#endif
#else 
#error must be called with COMPILE_INLINE_STATIC
#endif

/* imagsl gsl precision h */
/* imagsl gsl mode h */

const double gsl_prec_eps[_GSL_PREC_T_NUM] = {
  GSL_DBL_EPSILON,
  GSL_FLT_EPSILON,
  GSL_SFLT_EPSILON
};

const double gsl_prec_sqrt_eps[_GSL_PREC_T_NUM] = {
  GSL_SQRT_DBL_EPSILON,
  GSL_SQRT_FLT_EPSILON,
  GSL_SQRT_SFLT_EPSILON
};

const double gsl_prec_root3_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT3_DBL_EPSILON,
  GSL_ROOT3_FLT_EPSILON,
  GSL_ROOT3_SFLT_EPSILON
};

const double gsl_prec_root4_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT4_DBL_EPSILON,
  GSL_ROOT4_FLT_EPSILON,
  GSL_ROOT4_SFLT_EPSILON
};

const double gsl_prec_root5_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT5_DBL_EPSILON,
  GSL_ROOT5_FLT_EPSILON,
  GSL_ROOT5_SFLT_EPSILON
};

const double gsl_prec_root6_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT6_DBL_EPSILON,
  GSL_ROOT6_FLT_EPSILON,
  GSL_ROOT6_SFLT_EPSILON
};
#undef __GSL_BUILD_H__

/* #include <config.h> */
#include <stdlib.h>

#define BASE_GSL_COMPLEX_LONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_DOUBLE

#define BASE_FLOAT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_FLOAT

#define BASE_ULONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_ULONG

#define BASE_LONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_LONG

#define BASE_UINT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_UINT

#define BASE_INT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_INT

#define BASE_USHORT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_USHORT

#define BASE_SHORT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_SHORT

#define BASE_UCHAR
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_UCHAR

#define BASE_CHAR
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/init_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

TYPE (gsl_block) *
FUNCTION (gsl_block, alloc) (const size_t n)
{
  TYPE (gsl_block) * b;

  if (n == 0)
    {
      GSL_ERROR_VAL ("block length n must be positive integer",
                        GSL_EINVAL, 0);
    }

  b = (TYPE (gsl_block) *) malloc (sizeof (TYPE (gsl_block)));

  if (b == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for block struct",
                        GSL_ENOMEM, 0);
    }

  b->data = (ATOMIC *) malloc (MULTIPLICITY * n * sizeof (ATOMIC));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for block data",
                        GSL_ENOMEM, 0);
    }

  b->size = n;

  return b;
}

TYPE (gsl_block) *
FUNCTION (gsl_block, calloc) (const size_t n)
{
  size_t i;

  TYPE (gsl_block) * b = FUNCTION (gsl_block, alloc) (n);

  if (b == 0)
    return 0;

  /* initialize block to zero */

  for (i = 0; i < MULTIPLICITY * n; i++)
    {
      b->data[i] = 0;
    }

  return b;
}

void
FUNCTION (gsl_block, free) (TYPE (gsl_block) * b)
{
  free (b->data);
  free (b);
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_CHAR
/* #include <config.h> */

#define BASE_GSL_COMPLEX_LONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_DOUBLE

#define BASE_FLOAT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_FLOAT

#define BASE_ULONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_ULONG

#define BASE_LONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_LONG

#define BASE_UINT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_UINT

#define BASE_INT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_INT

#define BASE_USHORT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_USHORT

#define BASE_SHORT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_SHORT

#define BASE_UCHAR
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_UCHAR

#define BASE_CHAR
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/block_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

size_t
FUNCTION(gsl_block,size) (const TYPE(gsl_block) * b)
{
  return b->size ;
}

ATOMIC *
FUNCTION(gsl_block,data) (const TYPE(gsl_block) * b)
{
  return b->data ;
}
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_CHAR
/* #include <config.h> */
#include <stdio.h>

#define BASE_GSL_COMPLEX_LONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_DOUBLE

#define BASE_FLOAT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_FLOAT

#define BASE_ULONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_ULONG

#define BASE_LONG
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_LONG

#define BASE_UINT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_UINT

#define BASE_INT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_INT

#define BASE_USHORT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_USHORT

#define BASE_SHORT
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_SHORT

#define BASE_UCHAR
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_UCHAR

#define BASE_CHAR
/* If BASE is undefined we use function names like gsl_name()
   and assume that we are using doubles.

   If BASE is defined we used function names like gsl_BASE_name()
   and use BASE as the base datatype      */

#if   defined(BASE_GSL_COMPLEX_LONG)
#define BASE gsl_complex_long_double
#define SHORT complex_long_double
#define SHORT_REAL long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0L,0.0L}}
#define ONE {{1.0L,0.0L}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX)
#define BASE gsl_complex
#define SHORT complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0,0.0}}
#define ONE {{1.0,0.0}}
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_GSL_COMPLEX_FLOAT)
#define BASE gsl_complex_float
#define SHORT complex_float
#define SHORT_REAL float
#define ATOMIC float
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO {{0.0F,0.0F}}
#define ONE {{1.0F,0.0F}}
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_LONG_DOUBLE)
#define BASE long double
#define SHORT long_double
#define ATOMIC long double
#define USES_LONGDOUBLE 1
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%Lg"
#define OUT_FORMAT "%Lg"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0L
#define ONE 1.0L
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define ATOMIC double
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lg"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON GSL_DBL_EPSILON

#elif defined(BASE_FLOAT)
#define BASE float
#define SHORT float
#define ATOMIC float
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%g"
#define OUT_FORMAT "%g"
#define ATOMIC_IO ATOMIC
#define ZERO 0.0F
#define ONE 1.0F
#define BASE_EPSILON GSL_FLT_EPSILON

#elif defined(BASE_ULONG)
#define BASE unsigned long
#define SHORT ulong
#define ATOMIC unsigned long
#define MULTIPLICITY 1
#define IN_FORMAT "%lu"
#define OUT_FORMAT "%lu"
#define ATOMIC_IO ATOMIC
#define ZERO 0UL
#define ONE 1UL
#define UNSIGNED 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define ATOMIC long
#define MULTIPLICITY 1
#define IN_FORMAT "%ld"
#define OUT_FORMAT "%ld"
#define ATOMIC_IO ATOMIC
#define ZERO 0L
#define ONE 1L

#elif defined(BASE_UINT)
#define BASE unsigned int
#define SHORT uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define ATOMIC int
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_USHORT)
#define BASE unsigned short
#define SHORT ushort
#define ATOMIC unsigned short
#define MULTIPLICITY 1
#define IN_FORMAT "%hu"
#define OUT_FORMAT "%hu"
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_SHORT)
#define BASE short
#define SHORT short
#define ATOMIC short
#define MULTIPLICITY 1
#define IN_FORMAT "%hd"
#define OUT_FORMAT "%hd"
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1

#elif defined(BASE_UCHAR)
#define BASE unsigned char
#define SHORT uchar
#define ATOMIC unsigned char
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define ATOMIC_IO unsigned int
#define ZERO 0U
#define ONE 1U
#define UNSIGNED 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define ATOMIC char
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define ATOMIC_IO int
#define ZERO 0
#define ONE 1
#ifdef __CHAR_UNSIGNED__
#define UNSIGNED 1
#endif

#else
#error unknown BASE_ directive in source.h
#endif

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

/* FIXME: SANGCHUL */
#define CONCAT4COMPLEX(a,c,d) CONCAT4x(a,complex,c,d)
/* FIXME: SANGCHUL */

#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#ifdef USE_QUALIFIER
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,QUALIFIER,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT4COMPLEX(a,QUALIFIER,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2x(dir,complex)
#define QUALIFIED_VIEW(dir,name) CONCAT4COMPLEX(dir,QUALIFIER,name)
#else
#define FUNCTION(a,c) CONCAT4(a,SHORT,QUALIFIER,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT)
#define QUALIFIED_VIEW(dir,name) CONCAT4(dir,SHORT,QUALIFIER,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) QUALIFIER CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT4(dir,SHORT_REAL,QUALIFIER,name)
#endif
#else
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT2(dir,name)
#define TYPE(dir) dir
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT2(dir,name)
#else

/* FIXME: SANGCHUL */
#if defined(BASE_GSL_COMPLEX)
#define FUNCTION(a,c) CONCAT3x(a,complex,c)
#define TYPE(dir) CONCAT2x(dir,complex)
#define VIEW(dir,name) CONCAT3x(dir,complex,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3x(dir,complex,name)
#else
#define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#define TYPE(dir) CONCAT2(dir,SHORT)
#define VIEW(dir,name) CONCAT3(dir,SHORT,name)
#define QUALIFIED_TYPE(dir) TYPE(dir)
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,SHORT,name)
#endif
/* FIXME: SANGCHUL */

#endif
#if defined(BASE_GSL_COMPLEX)
#define REAL_TYPE(dir) dir
#define REAL_VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_REAL_TYPE(dir) dir
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT2(dir,name)
#else
#define REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#define QUALIFIED_REAL_TYPE(dir) CONCAT2(dir,SHORT_REAL)
#define QUALIFIED_REAL_VIEW(dir,name) CONCAT3(dir,SHORT_REAL,name)
#endif
#endif

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

/* block/fwrite_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;

  size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fread failed", GSL_EFAILED);
    }
      
return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
  size_t n = b->size ;

  ATOMIC * data = b->data ;
  
  size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);
  
  if (items != n)
    {
      GSL_ERROR ("fwrite failed", GSL_EFAILED);
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
                                 const size_t n, const size_t stride)
{
  if (stride == 1)
    {
      size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fread failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fread (data + MULTIPLICITY * i * stride,
                               MULTIPLICITY * sizeof (ATOMIC), 1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fread failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
                                  const size_t n, const size_t stride)
{

  if (stride == 1)
    {
      size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

      if (items != n)
        {
          GSL_ERROR ("fwrite failed", GSL_EFAILED);
        }
    }
  else
    {
      size_t i;

      for (i = 0; i < n; i++)
        {
          size_t item = fwrite (data + MULTIPLICITY * i * stride,
                                MULTIPLICITY * sizeof (ATOMIC),
                                1, stream);
          if (item != 1)
            {
              GSL_ERROR ("fwrite failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}
/* block/fprintf_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION (gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
  size_t n = b->size ;
  
  ATOMIC * data = b->data ;

  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp ;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;
          
          data [MULTIPLICITY * i + k] = tmp;


          if (status != 1)
            {
              GSL_ERROR ("fscanf failed", GSL_EFAILED);
            }
        }
    }

  return GSL_SUCCESS;
}


int
FUNCTION (gsl_block, raw_fprintf) (FILE * stream, 
                                   const ATOMIC * data,
                                   const size_t n,
                                   const size_t stride, 
                                   const char *format)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      int status;

      for (k = 0; k < MULTIPLICITY; k++)
        {
          if (k > 0)
            {
              status = putc (' ', stream);

              if (status == EOF)
                {
                  GSL_ERROR ("putc failed", GSL_EFAILED);
                }
            }
          status = fprintf (stream,
                            format,
                            data[MULTIPLICITY * i * stride + k]);
          if (status < 0)
            {
              GSL_ERROR ("fprintf failed", GSL_EFAILED);
            }
        }

      status = putc ('\n', stream);

      if (status == EOF)
        {
          GSL_ERROR ("putc failed", GSL_EFAILED);
        }
    }

  return 0;
}

int
FUNCTION (gsl_block, raw_fscanf) (FILE * stream, 
                                  ATOMIC * data,
                                  const size_t n, 
                                  const size_t stride)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      int k;
      for (k = 0; k < MULTIPLICITY; k++)
        {
          ATOMIC_IO tmp;

          int status = fscanf (stream, IN_FORMAT, &tmp) ;

          data [MULTIPLICITY * i * stride + k] = tmp;

          if (status != 1)
            GSL_ERROR ("fscanf failed", GSL_EFAILED);
        }
    }

  return GSL_SUCCESS;
}

#endif
#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef REAL_TYPE
#undef REAL_TYPE
#endif

#ifdef QUALIFIED_TYPE
#undef QUALIFIED_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_TYPE
#undef QUALIFIED_REAL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND
#undef UNSIGNED

#ifdef FP
#undef FP
#endif
#undef  BASE_CHAR
