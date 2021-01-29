/* Test file for mpfr_asinu.

Copyright 2021 Free Software Foundation, Inc.
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

#define TEST_FUNCTION mpfr_asinu
#define ULONG_ARG2
#include "tgeneric.c"

int
main (void)
{
  mpfr_t x, y;
  int r, inex;

  tests_start_mpfr ();

  mpfr_init (x);
  mpfr_init (y);

  MPFR_SET_NAN(x);
  mpfr_asinu (y, x, 1, MPFR_RNDN);
  if (mpfr_nan_p (y) == 0)
    {
      printf ("Error: asinu (NaN, 1) != NaN\n");
      exit (1);
    }

  mpfr_set_ui (x, 2, MPFR_RNDN);
  mpfr_asinu (y, x, 1, MPFR_RNDN);
  if (mpfr_nan_p (y) == 0)
    {
      printf ("Error: asinu (2, 1) != NaN\n");
      exit (1);
    }

  mpfr_set_si (x, -2, MPFR_RNDN);
  mpfr_asinu (y, x, 1, MPFR_RNDN);
  if (mpfr_nan_p (y) == 0)
    {
      printf ("Error: asinu (-2, 1) != NaN\n");
      exit (1);
    }

  /* asinu (+0,u) = +0 */
  mpfr_set_ui (x, 0, MPFR_RNDN);
  mpfr_asinu (y, x, 1, MPFR_RNDN);
  if (MPFR_NOTZERO (y) || MPFR_IS_NEG (y))
    {
      printf ("Error: asinu(+0,1) != +0\n");
      exit (1);
    }

  /* asinu (-0,u) = -0 */
  mpfr_set_ui (x, 0, MPFR_RNDN);
  mpfr_neg (x, x, MPFR_RNDN);
  mpfr_asinu (y, x, 1, MPFR_RNDN);
  if (MPFR_NOTZERO (y) || MPFR_IS_POS (y))
    {
      printf ("Error: asinu(-0,1) != -0\n");
      exit (1);
    }

  /* asinu (1,u) = u/4 */
  for (r = 0; r < MPFR_RND_MAX; r++)
    {
      mpfr_set_si (x, 1, MPFR_RNDN); /* exact */
      mpfr_asinu (y, x, 17, (mpfr_rnd_t) r);
      mpfr_set_ui_2exp (x, 17, -2, (mpfr_rnd_t) r);
      if (mpfr_cmp (x, y))
        {
          printf ("Error: asinu(1,17) != 17/4 for rnd=%s\n",
                  mpfr_print_rnd_mode ((mpfr_rnd_t) r));
          exit (1);
        }
    }

  /* asinu (-1,u) = -u/4 */
  for (r = 0; r < MPFR_RND_MAX; r++)
    {
      mpfr_set_si (x, -1, MPFR_RNDN); /* exact */
      mpfr_asinu (y, x, 17, (mpfr_rnd_t) r);
      mpfr_set_si_2exp (x, -17, -2, (mpfr_rnd_t) r);
      if (mpfr_cmp (x, y))
        {
          printf ("Error: asinu(-1,17) != -17/4 for rnd=%s\n",
                  mpfr_print_rnd_mode ((mpfr_rnd_t) r));
          exit (1);
        }
    }

  /* asinu (1/2,u) = u/12 */
  mpfr_set_si_2exp (x, 1, -1, MPFR_RNDN); /* exact */
  inex = mpfr_asinu (y, x, 12, MPFR_RNDN);
  if (inex != 0 || mpfr_cmp_ui (y, 1) != 0)
    {
      printf ("Error: asinu(1/2,12) != 1 (inex=%d)\n", inex);
      mpfr_dump (y);
      exit (1);
    }

  /* asinu (-1/2,u) = -u/12 */
  mpfr_set_si_2exp (x, -1, -1, MPFR_RNDN); /* exact */
  inex = mpfr_asinu (y, x, 12, MPFR_RNDN);
  if (inex != 0 || mpfr_cmp_si (y, -1) != 0)
    {
      printf ("Error: asinu(-1/2,12) != -1 (inex=%d)\n", inex);
      mpfr_dump (y);
      exit (1);
    }

  test_generic (MPFR_PREC_MIN, 100, 1000);

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();
  return 0;
}
