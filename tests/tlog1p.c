/* Test file for mpfr_log1p.

Copyright 2001, 2002, 2003, 2004, 2005 Free Software Foundation.
Adapted from tsinh.c.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

#ifdef CHECK_EXTERNAL
static int
test_log1p (mpfr_ptr a, mpfr_srcptr b, mp_rnd_t rnd_mode)
{
  int res;
  int ok = rnd_mode == GMP_RNDN && mpfr_number_p (b) && mpfr_get_prec (a)>=53;
  if (ok)
    {
      mpfr_print_raw (b);
    }
  res = mpfr_log1p (a, b, rnd_mode);
  if (ok)
    {
      printf (" ");
      mpfr_print_raw (a);
      printf ("\n");
    }
  return res;
}
#else
#define test_log1p mpfr_log1p
#endif

#define TEST_FUNCTION test_log1p
#include "tgeneric.c"

static void
special (void)
{
  mpfr_t x;

  mpfr_init (x);
  
  mpfr_set_nan (x);
  test_log1p (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (x));

  mpfr_set_inf (x, -1);
  test_log1p (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (x));

  mpfr_set_inf (x, 1);
  test_log1p (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) > 0);

  mpfr_set_ui (x, 0, GMP_RNDN);
  test_log1p (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 0) == 0 && MPFR_IS_POS (x));
  mpfr_neg (x, x, GMP_RNDN);
  test_log1p (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 0) == 0 && MPFR_IS_NEG (x));

  mpfr_set_si (x, -1, GMP_RNDN);
  test_log1p (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) < 0);

  mpfr_set_si (x, -2, GMP_RNDN);
  test_log1p (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (x));

  mpfr_clear (x);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  special ();

  test_generic (2, 100, 50);

  tests_end_mpfr ();
  return 0;
}
