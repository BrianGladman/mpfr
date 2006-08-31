/* Miscellaneous support for test programs.

Copyright 2001, 2002, 2003, 2004, 2005, 2006 Free Software Foundation, Inc.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#ifdef HAVE_CONFIG_H
# if HAVE_CONFIG_H
#  include "config.h"     /* for a build within gmp */
# endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#if HAVE_SETLOCALE
#include <locale.h>
#endif

#if TIME_WITH_SYS_TIME
# include <sys/time.h>  /* for struct timeval */
# include <time.h>
#elif HAVE_SYS_TIME_H
#  include <sys/time.h>
#else
#  include <time.h>
#endif

#if HAVE_SYS_FPU_H
# include <sys/fpu.h>
#endif

#include "mpfr-test.h"

static void tests_rand_start (void);
static void tests_rand_end   (void);

/* We want to always import the function mpfr_dump inside the test
   suite, so that we can use it in GDB. But it doesn't work if we build
   a Windows DLL (initializer element is not a constant) */
#if !__GMP_LIBGMP_DLL
extern void (*dummy_func) (mpfr_srcptr);
void (*dummy_func)(mpfr_srcptr) = mpfr_dump;
#endif

void
tests_start_mpfr (void)
{
  /* don't buffer, so output is not lost if a test causes a segv etc */
  setbuf (stdout, NULL);

#if HAVE_SETLOCALE
  /* Added on 2005-07-09. This allows to test MPFR under various
     locales. New bugs will probably be found, in particular with
     LC_ALL="tr_TR.ISO8859-9" because of the i/I character... */
  setlocale (LC_ALL, "");
#endif

  tests_memory_start ();
  tests_rand_start ();
}

void
tests_end_mpfr (void)
{
  mpfr_free_cache ();
  tests_rand_end ();
  tests_memory_end ();
}

static void
tests_rand_start (void)
{
  gmp_randstate_ptr  rands;
  char           *perform_seed;
  unsigned long  seed;

  if (__gmp_rands_initialized)
    {
      printf ("Please let tests_start() initialize the global __gmp_rands.\n");
      printf ("ie. ensure that function is called before the first use of RANDS.\n");
      abort ();
    }

  gmp_randinit_default (__gmp_rands);
  __gmp_rands_initialized = 1;
  rands = __gmp_rands;

  perform_seed = getenv ("GMP_CHECK_RANDOMIZE");
  if (perform_seed != NULL)
    {
      seed = atoi (perform_seed);
      if (! (seed == 0 || seed == 1))
        {
          printf ("Re-seeding with GMP_CHECK_RANDOMIZE=%lu\n", seed);
          gmp_randseed_ui (rands, seed);
        }
      else
        {
#if HAVE_GETTIMEOFDAY
          struct timeval  tv;
          gettimeofday (&tv, NULL);
          seed = tv.tv_sec + tv.tv_usec;
#else
          time_t  tv;
          time (&tv);
          seed = tv;
#endif
          gmp_randseed_ui (rands, seed);
          printf ("Seed GMP_CHECK_RANDOMIZE=%lu (include this in bug reports)\n", seed);
        }
    }
}

static void
tests_rand_end (void)
{
  RANDS_CLEAR ();
}

/* initialization function for tests using the hardware floats
   Not very useful now. */
void
mpfr_test_init ()
{
  double d;
#if HAVE_FPC_CSR
  /* to get denormalized numbers on IRIX64 */
  union fpc_csr exp;

  exp.fc_word = get_fpc_csr();
  exp.fc_struct.flush = 0;
  set_fpc_csr(exp.fc_word);
#endif
#ifdef HAVE_DENORMS
  d = DBL_MIN;
  if (2.0 * (d / 2.0) != d)
    {
      printf ("Warning: no denormalized numbers\n");
      exit (1);
    }
#endif

  /* generate DBL_EPSILON with a loop to avoid that the compiler
     optimizes the code below in non-IEEE 754 mode, deciding that
     c = d is always false. */
#if 0
  for (eps = 1.0; eps != DBL_EPSILON; eps /= 2.0);
  c = 1.0 + eps;
  d = eps * (1.0 - eps) / 2.0;
  d += c;
  if (c != d)
    {
      printf ("Warning: IEEE 754 standard not fully supported\n"
              "         (maybe extended precision not disabled)\n"
              "         Some tests may fail\n");
    }
#endif
}


/* generate a random limb */
mp_limb_t
randlimb (void)
{
  mp_limb_t limb;

  _gmp_rand (&limb, RANDS, BITS_PER_MP_LIMB);
  return limb;
}

/* returns ulp(x) for x a 'normal' double-precision number */
double
Ulp (double x)
{
   double y, eps;

   if (x < 0) x = -x;

   y = x * 2.220446049250313080847263336181640625e-16 ; /* x / 2^52 */

   /* as ulp(x) <= y = x/2^52 < 2*ulp(x),
   we have x + ulp(x) <= x + y <= x + 2*ulp(x),
   therefore o(x + y) = x + ulp(x) or x + 2*ulp(x) */

   eps =  x + y;
   eps = eps - x; /* ulp(x) or 2*ulp(x) */

   return (eps > y) ? 0.5 * eps : eps;
}

/* returns the number of ulp's between a and b,
   where a and b can be any floating-point number, except NaN
 */
int
ulp (double a, double b)
{
  double twoa;

  if (a == b) return 0; /* also deals with a=b=inf or -inf */

  twoa = a + a;
  if (twoa == a) /* a is +/-0.0 or +/-Inf */
    return ((b < a) ? INT_MAX : -INT_MAX);

  return (int) ((a - b) / Ulp (a));
}

/* return double m*2^e */
double
dbl (double m, int e)
{
  if (e >=0 )
    while (e-- > 0)
      m *= 2.0;
  else
    while (e++ < 0)
      m /= 2.0;
  return m;
}

int
Isnan (double d)
{
  return (d) != (d);
}

void
d_trace (const char *name, double d)
{
  union {
    double         d;
    unsigned char  b[sizeof(double)];
  } u;
  int  i;

  if (name != NULL && name[0] != '\0')
    printf ("%s=", name);

  u.d = d;
  printf ("[");
  for (i = 0; i < (int) sizeof (u.b); i++)
    {
      if (i != 0)
        printf (" ");
      printf ("%02X", (int) u.b[i]);
    }
  printf ("] %.20g\n", d);
}

void
ld_trace (const char *name, long double ld)
{
  union {
    long double    ld;
    unsigned char  b[sizeof(long double)];
  } u;
  int  i;

  if (name != NULL && name[0] != '\0')
    printf ("%s=", name);

  u.ld = ld;
  printf ("[");
  for (i = 0; i < (int) sizeof (u.b); i++)
    {
      if (i != 0)
        printf (" ");
      printf ("%02X", (int) u.b[i]);
    }
  printf ("] %.20Lg\n", ld);
}

/* Open a file in the src directory - can't use fopen directly */
FILE *
src_fopen (const char *filename, const char *mode)
{
  const char *srcdir = getenv ("srcdir");
  char *buffer;
  FILE *f;

  if (srcdir == NULL)
    return fopen (filename, mode);
  buffer = (char*) malloc (strlen (filename) + strlen (srcdir) + 2);
  if (buffer == NULL)
    {
      printf ("src_fopen: failed to alloc memory)\n");
      exit (1);
    }
  sprintf (buffer, "%s/%s", srcdir, filename);
  f = fopen (buffer, mode);
  free (buffer);
  return f;
}

void
set_emin (mp_exp_t exponent)
{
  if (mpfr_set_emin (exponent))
    {
      printf ("set_emin: setting emin to %ld failed\n", (long int) exponent);
      exit (1);
    }
}

void
set_emax (mp_exp_t exponent)
{
  if (mpfr_set_emax (exponent))
    {
      printf ("set_emax: setting emax to %ld failed\n", (long int) exponent);
      exit (1);
    }
}
