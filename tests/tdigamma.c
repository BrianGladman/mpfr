/* test file for digamma function

Copyright 2009-2021 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_digamma
#include "tgeneric.c"

static void
special (void)
{
  mpfr_t x, y;

  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_inf (y, -1);
  mpfr_set_inf (x, 1);
  mpfr_digamma (y, x, MPFR_RNDN);
  if (mpfr_inf_p (y) == 0 || mpfr_sgn (y) < 0)
    {
      printf ("error for Psi(+Inf)\n");
      printf ("expected +Inf\n");
      printf ("got      ");
      mpfr_dump (y);
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

/* With some GMP_CHECK_RANDOMIZE values, test_generic triggers an error
     tests_addsize(): too much memory (576460752303432776 bytes)
  Each time on prec = 200, n = 3, xprec = 140.
  The following test is a more general testcase.
*/
static void
bug20210206 (void)
{
  mpfr_t x, y[2], z;
  mpfr_exp_t emin, emax;
  int i, precx, precy[2] = { 200, 400 };

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();
  set_emin (MPFR_EMIN_MIN);
  set_emax (MPFR_EMAX_MAX);

  for (i = 0; i < 2; i++)
    mpfr_init2 (y[i], precy[i]);
  mpfr_init2 (z, precy[0]);

  for (precx = MPFR_PREC_MIN; precx < 150; precx++)
    {
      mpfr_init2 (x, precx);
      mpfr_setmax (x, __gmpfr_emax);
      for (i = 0; i < 2; i++)
        mpfr_digamma (y[i], x, MPFR_RNDA);
      mpfr_set (z, y[1], MPFR_RNDA);
      MPFR_ASSERTN (mpfr_equal_p (y[0], z));
      mpfr_clear (x);
    }

  for (i = 0; i < 2; i++)
    mpfr_clear (y[i]);
  mpfr_clear (z);

  set_emin (emin);
  set_emax (emax);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  special ();
  bug20210206 ();

  test_generic (MPFR_PREC_MIN, 200, 20);

  data_check ("data/digamma", mpfr_digamma, "mpfr_digamma");

  tests_end_mpfr ();
  return 0;
}
