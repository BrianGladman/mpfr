/* Test file for mpfr_sub_ui

Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "mpfr-test.h"

/* checks that x-y gives the right results with 53 bits of precision */
static void
check3 (const char *xs, unsigned long y, mp_rnd_t rnd_mode, const char *zs)
{
  mpfr_t xx,zz;

  mpfr_inits2 (53, xx, zz, (void *) 0);
  mpfr_set_str1 (xx, xs);
  mpfr_sub_ui (zz, xx, y, rnd_mode);
  if (mpfr_cmp_str1(zz, zs))
    {
      printf ("expected sum is %s, got ", zs);
      mpfr_print_binary(zz);
      printf ("\nmpfr_sub_ui failed for x=%s y=%lu with rnd_mode=%s\n",
              xs, y, mpfr_print_rnd_mode (rnd_mode));
      exit (1);
    }
  mpfr_clears (xx, zz, (void *) 0);
}

/* FastTwoSum: if EXP(x) >= EXP(y), u = o(x+y), v = o(u-x), w = o(y-v),
               then x + y = u + w
thus if u = o(y-x), v = o(u+x), w = o(v-y), then y-x = u-w */
static void
check_two_sum (mp_prec_t p)
{
  unsigned int x;
  mpfr_t y, u, v, w;
  mp_rnd_t rnd;
  int inexact;

  mpfr_inits2 (p, y, u, v, w, (void *) 0);
  do
    {
      x = randlimb ();
    }
  while (x < 1);
  mpfr_random (y);
  rnd = GMP_RNDN;
  inexact = mpfr_sub_ui (u, y, x, rnd);
  mpfr_add_ui (v, u, x, rnd);
  mpfr_sub (w, v, y, rnd);
  /* as u - (y-x) = w, we should have inexact and w of same sign */
  if (((inexact == 0) && mpfr_cmp_ui (w, 0)) ||
      ((inexact > 0) && (mpfr_cmp_ui (w, 0) <= 0)) ||
      ((inexact < 0) && (mpfr_cmp_ui (w, 0) >= 0)))
    {
      printf ("Wrong inexact flag for prec=%u, rnd=%s\n",
              (unsigned int) p, mpfr_print_rnd_mode (rnd));
      printf ("x=%u\n", x);
      printf ("y="); mpfr_print_binary(y); puts ("");
      printf ("u="); mpfr_print_binary(u); puts ("");
      printf ("v="); mpfr_print_binary(v); puts ("");
      printf ("w="); mpfr_print_binary(w); puts ("");
      printf ("inexact = %d\n", inexact);
      exit (1);
    }
  mpfr_clears (y, u, v, w, (void *) 0);
}

static void
check_nans (void)
{
  mpfr_t  x, y;

  mpfr_init2 (x, 123L);
  mpfr_init2 (y, 123L);

  /* nan - 1 == nan */
  mpfr_set_nan (x);
  mpfr_sub_ui (y, x, 1L, GMP_RNDN);
  MPFR_ASSERTN (mpfr_nan_p (y));

  /* +inf - 1 == +inf */
  mpfr_set_inf (x, 1);
  mpfr_sub_ui (y, x, 1L, GMP_RNDN);
  MPFR_ASSERTN (mpfr_inf_p (y));
  MPFR_ASSERTN (mpfr_sgn (y) > 0);

  /* -inf - 1 == -inf */
  mpfr_set_inf (x, -1);
  mpfr_sub_ui (y, x, 1L, GMP_RNDN);
  MPFR_ASSERTN (mpfr_inf_p (y));
  MPFR_ASSERTN (mpfr_sgn (y) < 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

#define TEST_FUNCTION mpfr_sub_ui
#define INTEGER_TYPE  unsigned long
#define RAND_FUNCTION(x) mpfr_random2(x, MPFR_LIMB_SIZE (x), 1)
#include "tgeneric_ui.c"

int
main (int argc, char *argv[])
{
  mp_prec_t p;
  int k;

  MPFR_TEST_USE_RANDS ();
  tests_start_mpfr ();

  check_nans ();

  for (p=2; p<200; p++)
    for (k=0; k<200; k++)
      check_two_sum (p);

  check3 ("0.9999999999", 1, GMP_RNDN,
          "-10000000827403709990903735160827636718750e-50");

  test_generic_ui (2, 1000, 100);

  tests_end_mpfr ();
  return 0;
}
