/* Test file for mpfr_fmma and mpfr_fmms.

Copyright 2016 Free Software Foundation, Inc.
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

/* TODO: add more tests, with special values and exception checking. */

/* check both mpfr_fmma and mpfr_fmms */
static void
random_test (mpfr_t a, mpfr_t b, mpfr_t c, mpfr_t d, mpfr_rnd_t rnd)
{
  mpfr_t ref, res, ab, cd;
  int inex_ref, inex_res;
  mpfr_prec_t p = MPFR_PREC(a);

  mpfr_init2 (res, p);
  mpfr_init2 (ref, p);
  mpfr_init2 (ab, mpfr_get_prec (a) +  mpfr_get_prec (b));
  mpfr_init2 (cd, mpfr_get_prec (c) +  mpfr_get_prec (d));

  /* first check fmma */
  inex_res = mpfr_fmma (res, a, b, c, d, rnd);
  mpfr_mul (ab, a, b, rnd);
  mpfr_mul (cd, c, d, rnd);
  inex_ref = mpfr_add (ref, ab, cd, rnd);
  if (inex_res != inex_ref)
    {
      printf ("mpfr_fmma failed for p=%ld rnd=%s\n",
              (long int) p, mpfr_print_rnd_mode (rnd));
      printf ("a="); mpfr_dump (a);
      printf ("b="); mpfr_dump (b);
      printf ("c="); mpfr_dump (c);
      printf ("d="); mpfr_dump (d);
      printf ("expected inex %d, got %d\n", inex_ref, inex_res);
      exit (1);
    }
  if (mpfr_nan_p (res))
    MPFR_ASSERTN (mpfr_nan_p (ref));
  else
    MPFR_ASSERTN (mpfr_cmp (res, ref) == 0);

  /* then check fmms */
  inex_res = mpfr_fmms (res, a, b, c, d, rnd);
  mpfr_mul (ab, a, b, rnd);
  mpfr_mul (cd, c, d, rnd);
  inex_ref = mpfr_sub (ref, ab, cd, rnd);
  if (inex_res != inex_ref)
    {
      printf ("mpfr_fmms failed for p=%ld rnd=%s\n",
              (long int) p, mpfr_print_rnd_mode (rnd));
      printf ("a="); mpfr_dump (a);
      printf ("b="); mpfr_dump (b);
      printf ("c="); mpfr_dump (c);
      printf ("d="); mpfr_dump (d);
      printf ("expected inex %d, got %d\n", inex_ref, inex_res);
      exit (1);
    }
  if (mpfr_nan_p (res))
    MPFR_ASSERTN (mpfr_nan_p (ref));
  else
    MPFR_ASSERTN (mpfr_cmp (res, ref) == 0);

  mpfr_clear (ab);
  mpfr_clear (cd);
  mpfr_clear (res);
  mpfr_clear (ref);
}

static void
random_tests (void)
{
  mpfr_prec_t p;
  int r;
  mpfr_t a, b, c, d;

  for (p = MPFR_PREC_MIN; p <= 4096; p++)
    {
      mpfr_inits2 (p, a, b, c, d, (mpfr_ptr) 0);
      mpfr_urandomb (a, RANDS);
      mpfr_urandomb (b, RANDS);
      mpfr_urandomb (c, RANDS);
      mpfr_urandomb (d, RANDS);
      RND_LOOP (r)
        random_test (a, b, c, d, (mpfr_rnd_t) r);
      mpfr_clears (a, b, c, d, (mpfr_ptr) 0);
    }
}

static void
zero_tests (void)
{
  mpfr_exp_t emin, emax;
  mpfr_t a, b, c, d, res;
  int i, r;

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();
  set_emin (MPFR_EMIN_MIN);
  set_emax (MPFR_EMAX_MAX);

  mpfr_inits2 (64, a, b, c, d, res, (mpfr_ptr) 0);
  for (i = 0; i <= 4; i++)
    {
      switch (i)
        {
        case 0: case 1:
          mpfr_set_ui (a, i, MPFR_RNDN);
          break;
        case 2:
          mpfr_setmax (a, MPFR_EMAX_MAX);
          break;
        case 3:
          mpfr_setmin (a, MPFR_EMIN_MIN);
          break;
        case 4:
          mpfr_setmin (a, MPFR_EMIN_MIN+1);
          break;
        }
      RND_LOOP (r)
        {
          int inex;
          mpfr_flags_t flags;

          mpfr_set (b, a, MPFR_RNDN);
          mpfr_set (c, a, MPFR_RNDN);
          mpfr_neg (d, a, MPFR_RNDN);
          mpfr_clear_flags ();
          inex = mpfr_fmma (res, a, b, c, d, (mpfr_rnd_t) r);
          flags = __gmpfr_flags;
          if (flags != 0 || inex != 0 || ! mpfr_zero_p (res) ||
              (r == MPFR_RNDD ? MPFR_IS_POS (res) : MPFR_IS_NEG (res)))
            {
              printf ("Error in zero_tests on i = %d, %s\n",
                      i, mpfr_print_rnd_mode ((mpfr_rnd_t) r));
              printf ("Expected %c0, inex = 0\n", r == MPFR_RNDD ? '-' : '+');
              printf ("Got      ");
              if (MPFR_IS_POS (res))
                printf ("+");
              mpfr_out_str (stdout, 16, 0, res, MPFR_RNDN);
              printf (", inex = %d\n", inex);
              printf ("Expected flags:");
              flags_out (0);
              printf ("Got flags:     ");
              flags_out (flags);
              exit (1);
            }
        }
    }
  mpfr_clears (a, b, c, d, res, (mpfr_ptr) 0);

  set_emin (emin);
  set_emax (emax);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  random_tests ();
  zero_tests ();

  tests_end_mpfr ();
  return 0;
}
