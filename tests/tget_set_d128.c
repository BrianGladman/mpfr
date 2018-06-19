/* Test file for mpfr_set_decimal128 (and mpfr_get_decimal128 when available).

Copyright 2018 Free Software Foundation, Inc.
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

/* Needed due to the test on MPFR_WANT_DECIMAL_FLOATS */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef MPFR_WANT_DECIMAL_FLOATS

#include "mpfr-test.h"

static void
test_set (void)
{
  long v[] = { 1, -1, 2147483647, -2147483647 };
  mpfr_t x;
  mpfr_flags_t flags;
  int i, inex;

  mpfr_init2 (x, 53);
  for (i = 0; i < numberof (v); i++)
    {
      mpfr_clear_flags ();
      inex = mpfr_set_decimal128 (x, (_Decimal128) v[i], MPFR_RNDN);
      flags = __gmpfr_flags;
      if (mpfr_cmp_si (x, v[i]) != 0 || inex != 0 || flags != 0)
        {
          printf ("Error in test_set for i=%d\n", i);
          printf ("Expected %ld\n    with inex = 0 and flags =", v[i]);
          flags_out (0);
          printf ("Got      ");
          mpfr_dump (x);
          printf ("    with inex = %d and flags =", inex);
          flags_out (flags);
          exit (1);
        }
    }
  mpfr_clear (x);
}

static void
powers_of_10 (void)
{
  mpfr_t x1, x2;
  int inex1, inex2;
  mpfr_flags_t flags1, flags2;
  _Decimal128 d[2];
  int i, neg, rnd;

  mpfr_inits2 (200, x1, x2, (mpfr_ptr) 0);
  for (i = 0, d[0] = 1, d[1] = 1; i < 150; i++, d[0] *= 10, d[1] /= 10)
    for (neg = 0; neg <= 1; neg++)
      RND_LOOP_NO_RNDF (rnd)
        {
          inex1 = mpfr_set_si (x1, neg ? -i : i, MPFR_RNDN);
          MPFR_ASSERTN (inex1 == 0);

          mpfr_clear_flags ();
          inex1 = mpfr_exp10 (x1, x1, (mpfr_rnd_t) rnd);
          flags1 = __gmpfr_flags;

          mpfr_clear_flags ();
          inex2 = mpfr_set_decimal128 (x2, d[neg], (mpfr_rnd_t) rnd);
          flags2 = __gmpfr_flags;

          if (!(mpfr_equal_p (x1, x2) &&
                SAME_SIGN (inex1, inex2) &&
                flags1 == flags2))
            {
              printf ("Error in powers_of_10 for i=%d, neg=%d, %s\n",
                      i, neg, mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
              printf ("Expected ");
              mpfr_dump (x1);
              printf ("with inex = %d and flags =", inex1);
              flags_out (flags1);
              printf ("Got      ");
              mpfr_dump (x2);
              printf ("with inex = %d and flags =", inex2);
              flags_out (flags2);
              exit (1);
            }
      }
  mpfr_clears (x1, x2, (mpfr_ptr) 0);
}

int
main (int argc, char *argv[])
{
  int verbose = argc > 1;

  tests_start_mpfr ();
  mpfr_test_init ();

  if (verbose)
    {
#ifdef DPD_FORMAT
      printf ("Using DPD format\n");
#else
      printf ("Using BID format\n");
#endif
    }

  test_set ();
  powers_of_10 ();

  tests_end_mpfr ();
  return 0;
}

#else /* MPFR_WANT_DECIMAL_FLOATS */

int
main (void)
{
  return 77;
}

#endif /* MPFR_WANT_DECIMAL_FLOATS */
