/* Test file for mpfr_get_q.

Copyright 2017 Free Software Foundation, Inc.
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
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "mpfr-test.h"

#ifndef MPFR_USE_MINI_GMP

static void
special (void)
{
  mpfr_t f;
  mpq_t q;

  mpfr_init2 (f, MPFR_PREC_MIN);
  mpq_init (q);

  /* check NaN */
  mpfr_set_nan (f);
  mpfr_clear_erangeflag ();
  mpfr_get_q (q, f);
  MPFR_ASSERTN(mpq_cmp_ui (q, 0, 1) == 0);
  MPFR_ASSERTN(mpfr_erangeflag_p ());

  /* check +Inf */
  mpfr_set_inf (f, 1);
  mpfr_clear_erangeflag ();
  mpfr_get_q (q, f);
  MPFR_ASSERTN(mpq_cmp_ui (q, 0, 1) == 0);
  MPFR_ASSERTN(mpfr_erangeflag_p ());

  /* check -Inf */
  mpfr_set_inf (f, -1);
  mpfr_clear_erangeflag ();
  mpfr_get_q (q, f);
  MPFR_ASSERTN(mpq_cmp_ui (q, 0, 1) == 0);
  MPFR_ASSERTN(mpfr_erangeflag_p ());

  /* check +0 */
  mpfr_set_zero (f, 1);
  mpfr_clear_erangeflag ();
  mpfr_get_q (q, f);
  MPFR_ASSERTN(mpq_cmp_ui (q, 0, 1) == 0);
  MPFR_ASSERTN(!mpfr_erangeflag_p ());

  /* check -0 */
  mpfr_set_zero (f, -1);
  mpfr_clear_erangeflag ();
  mpfr_get_q (q, f);
  MPFR_ASSERTN(mpq_cmp_ui (q, 0, 1) == 0);
  MPFR_ASSERTN(!mpfr_erangeflag_p ());

  mpq_clear (q);
  mpfr_clear (f);
}

static void
random_tests (void)
{
  mpfr_t f, g;
  mpq_t q;
  int inex;
  mpfr_rnd_t rnd;
  int i;

  mpfr_init2 (f, MPFR_PREC_MIN + (randlimb() % 100));
  mpfr_init2 (g, mpfr_get_prec (f));
  mpq_init (q);

  for (i = 0; i < 1000; i++)
    {
      mpfr_urandomb (f, RANDS);
      mpfr_get_q (q, f);
      rnd = RND_RAND ();
      inex = mpfr_set_q (g, q, rnd);
      MPFR_ASSERTN(inex == 0);
      MPFR_ASSERTN(mpfr_cmp (f, g) == 0);
    }

  mpq_clear (q);
  mpfr_clear (f);
  mpfr_clear (g);
}

int
main (void)
{
  tests_start_mpfr ();

  special ();
  random_tests ();

  tests_end_mpfr ();
  return 0;
}

#else

int
main (void)
{
  return 77;
}

#endif /* MPFR_USE_MINI_GMP */
