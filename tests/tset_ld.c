/* Test file for mpfr_set_ld and mpfr_get_ld.

Copyright 2002, 2003, 2004 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

#include "mpfr-test.h"

static void
check_gcc33_bug (void)
{
  volatile long double x;
  x = (long double) 9007199254740992.0 + 1.0;
  if (x != 0.0)
    return;  /* OK */
  printf
    ("Detected optimization bug of gcc 3.3 on Alpha concerning long double\n"
     "comparisons; set_ld tests are disabled (set_ld won't work correctly).\n"
     "See http://gcc.gnu.org/ml/gcc-bugs/2003-10/msg00853.html for more\n"
     "information on this bug.\n");
  exit (0);  /* This is not a bug in MPFR, so don't fail. */
}

static int
Isnan_ld (long double d)
{
  LONGDOUBLE_NAN_ACTION (d, goto yes);
  return 0;
 yes:
  return 1;
}

/* checks that a long double converted to a mpfr (with precision >=113),
   then converted back to a long double gives the initial value,
   or in other words mpfr_get_ld(mpfr_set_ld(d)) = d.
*/
static void
check_set_get (long double d, mpfr_t x)
{
  mp_rnd_t r;
  long double e;
  int inex;

  for (r = 0; r < GMP_RND_MAX; r++)
    {
      inex = mpfr_set_ld (x, d, r);
      if (inex != 0)
        {
          printf ("Error: mpfr_set_ld should be exact\n");
          printf ("d=%1.30Le inex=%d\n", d, inex);
          printf ("emin=%ld emax=%ld\n", mpfr_get_emin (), mpfr_get_emax ());
          mpfr_dump (x);
          exit (1);
        }
      e = mpfr_get_ld (x, r);
      if (e != d && !(Isnan_ld(e) && Isnan_ld(d)))
        {
          printf ("Error: mpfr_get_ld o mpfr_set_ld <> Id\n");
          printf ("  r=%d\n", r);
          printf ("  d=%1.30Le get_ld(set_ld(d))=%1.30Le\n", d, e);
          ld_trace ("  d", d);
          printf ("  x="); mpfr_out_str (NULL, 16, 0, x, GMP_RNDN);
          printf ("\n");
          ld_trace ("  e", e);
          exit (1);
        }
    }
}

int
main (int argc, char *argv[])
{
  long double d, e;
  mpfr_t x;
  int i;
  mp_exp_t emax;

  check_gcc33_bug ();

  tests_start_mpfr ();
  mpfr_test_init ();
  tests_machine_prec_long_double ();

  mpfr_init2 (x, 113);

  /* check +0.0 and -0.0 */
  d = 0.0;
  check_set_get (d, x);
  d = -0.0;
  check_set_get (d, x);

  /* checks that sign of -0.0 is set */
  mpfr_set_ld (x, -0.0, GMP_RNDN);
  if (MPFR_SIGN(x) > 0)
    {
      printf ("Error: sign of -0.0 is not set correctly\n");
#if _GMP_IEEE_FLOATS
      exit (1);
      /* Non IEEE doesn't support negative zero yet */
#endif
    }

  /* checks NaN, Inf and -Inf */
  mpfr_set_nan (x);
  d = mpfr_get_ld (x, GMP_RNDN);
  check_set_get (d, x);

  mpfr_set_inf (x, 1);
  d = mpfr_get_ld (x, GMP_RNDN);
  check_set_get (d, x);

  mpfr_set_inf (x, -1);
  d = mpfr_get_ld (x, GMP_RNDN);
  check_set_get (d, x);

  /* checks the largest power of two */
  d = 1.0; while (d < LDBL_MAX / 2.0) d += d;
  check_set_get (d, x);
  check_set_get (-d, x);

  /* checks largest long double */
  d = LDBL_MAX;
  check_set_get (d, x);
  check_set_get (-d, x);

  /* checks the smallest power of two */
  d = 1.0; while ((e = d / 2.0) != (long double) 0.0) d = e;
  check_set_get (d, x);
  check_set_get (-d, x);

  /* checks largest 2^(2^k) that is representable as a long double */
  d = (LDBL_MAX / 2) + (LDBL_MAX / 4 * LDBL_EPSILON);
  check_set_get (d, x);

  /* checks that 2^i, 2^i+1 and 2^i-1 are correctly converted */
  d = 1.0;
  for (i = 1; i <= 113; i++)
    {
      d = 2.0 * d; /* d = 2^i */
      check_set_get (d, x);
      check_set_get (d + 1.0, x);
      check_set_get (d - 1.0, x);
    }

  for (i = 0; i < 10000; i++)
    {
      mpfr_random (x);
      d = mpfr_get_ld (x, GMP_RNDN);
      check_set_get (d, x);
    }

  /* check with reduced emax to exercise overflow */
  emax = mpfr_get_emax ();
  mpfr_set_prec (x, 2);
  set_emax (1);
  mpfr_set_ld (x, (long double) 2.0, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) > 0);
  for (d = (long double) 2.0, i = 0; i < 13; i++, d *= d);
  /* now d = 2^8192 */
  mpfr_set_ld (x, d, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) > 0);
  set_emax (emax);

  mpfr_clear (x);

  tests_end_mpfr ();

  return 0;
}
