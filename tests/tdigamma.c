/* test file for digamma function

Copyright 2009-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

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
#define NPREC 4
  mpfr_t x, y[NPREC], z;
  mpfr_exp_t emin, emax;
  int i, precx, precy[NPREC] = { 200, 400, 520, 1416 };

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();
  set_emin (MPFR_EMIN_MIN);
  set_emax (MPFR_EMAX_MAX);

  for (i = 0; i < NPREC; i++)
    mpfr_init2 (y[i], precy[i]);
  mpfr_init2 (z, precy[0]);

  for (precx = MPFR_PREC_MIN; precx < 150; precx++)
    {
      mpfr_init2 (x, precx);
      mpfr_setmax (x, __gmpfr_emax);
      for (i = 0; i < NPREC; i++)
        mpfr_digamma (y[i], x, MPFR_RNDA);
      mpfr_set (z, y[1], MPFR_RNDA);
      MPFR_ASSERTN (mpfr_equal_p (y[0], z));
      mpfr_clear (x);
    }

  for (i = 0; i < NPREC; i++)
    mpfr_clear (y[i]);
  mpfr_clear (z);

  set_emin (emin);
  set_emax (emax);
}

/* another test that fails with GMP_CHECK_RANDOMIZE=1612741376857003
   on revision 14398 */
static void
bug20210208 (void)
{
  mpfr_t x, y;
  int inex;

  mpfr_init2 (x, 73);
  mpfr_init2 (y, 1);
  mpfr_set_str (x, "1.4613470547060071827450", 10, MPFR_RNDN);
  mpfr_clear_flags ();
  inex = mpfr_digamma (y, x, MPFR_RNDU);
  MPFR_ASSERTN (mpfr_cmp_si_2exp (y, -1, -12) == 0);
  MPFR_ASSERTN (inex > 0);
  MPFR_ASSERTN (__gmpfr_flags == MPFR_FLAGS_INEXACT);
  mpfr_clear (x);
  mpfr_clear (y);
}

/* another test that fails with GMP_CHECK_RANDOMIZE=1613197421465830
   on revision 14429 */
static void
bug20210215 (void)
{
  mpfr_t x, y;
  int inex;

  mpfr_init2 (x, 510);
  mpfr_init2 (y, 4);
  mpfr_set_str (x, "-8.2923051438433494998166335341807999322052669984208422481227138906096000469898717007386115912802685588348601663465077353194268894939972221117314512518182580e+35", 10, MPFR_RNDN);
  mpfr_clear_flags ();
  inex = mpfr_digamma (y, x, MPFR_RNDU);
  MPFR_ASSERTN (mpfr_cmp_ui0 (y, 88) == 0);
  MPFR_ASSERTN (inex > 0);
  MPFR_ASSERTN (__gmpfr_flags == MPFR_FLAGS_INEXACT);
  mpfr_clear (x);
  mpfr_clear (y);
}

/* exercise the special case in the code for very small x */
static void
test_tiny (void)
{
  int p;

  for (p = MPFR_PREC_MIN; p <= 10; p++)
    {
      /* the special code needs EXP(x) < -2, thus |x| < 2^-3 */
      mpfr_t x, y, z, t;
      mpfr_init2 (x, p);
      mpfr_init2 (y, p);
      mpfr_init2 (z, p);
      mpfr_init2 (t, p + 20);
      mpfr_set_ui_2exp (x, 1, -5, MPFR_RNDN); /* x = 2^-5 */
      while (mpfr_cmp_ui_2exp (x, 1, -3) < 0)
        {
          int rnd;
          RND_LOOP (rnd)
            {
              mpfr_digamma (y, x, (mpfr_rnd_t) rnd);
              mpfr_digamma (t, x, MPFR_RNDN);
              if (mpfr_can_round (t, 0, MPFR_RNDN, (mpfr_rnd_t) rnd, p))
                {
                  mpfr_set (z, t, (mpfr_rnd_t) rnd);
                  if (! mpfr_equal_p (y, z))
                    {
                      printf ("Error for x=");
                      mpfr_dump (x);
                      printf ("Expected  z=");
                      mpfr_dump (z);
                      printf ("Got       y=");
                      mpfr_dump (z);
                      exit (1);
                    }
                }
              /* check also -x */
              mpfr_neg (x, x, MPFR_RNDN);
              mpfr_digamma (y, x, (mpfr_rnd_t) rnd);
              mpfr_digamma (t, x, MPFR_RNDN);
              if (mpfr_can_round (t, 0, MPFR_RNDN, (mpfr_rnd_t) rnd, p))
                {
                  mpfr_set (z, t, (mpfr_rnd_t) rnd);
                  if (! mpfr_equal_p (y, z))
                    {
                      printf ("Error for x=");
                      mpfr_dump (x);
                      printf ("Expected  z=");
                      mpfr_dump (z);
                      printf ("Got       y=");
                      mpfr_dump (z);
                      exit (1);
                    }
                }
              mpfr_neg (x, x, MPFR_RNDN); /* restore x > 0 */
            }
          mpfr_nextabove (x);
        }
      mpfr_clear (x);
      mpfr_clear (y);
      mpfr_clear (z);
      mpfr_clear (t);
    }
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  test_tiny ();
  special ();
  bug20210206 ();
  bug20210208 ();
  bug20210215 ();

  test_generic (MPFR_PREC_MIN, 200, 20);

  data_check ("data/digamma", mpfr_digamma, "mpfr_digamma");

  tests_end_mpfr ();
  return 0;
}
