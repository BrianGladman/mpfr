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
  if (! SAME_SIGN (inex_res, inex_ref) ||
      mpfr_nan_p (res) || mpfr_nan_p (ref) ||
      ! mpfr_equal_p (res, ref))
    {
      printf ("mpfr_fmma failed for p=%ld rnd=%s\n",
              (long int) p, mpfr_print_rnd_mode (rnd));
      printf ("a="); mpfr_dump (a);
      printf ("b="); mpfr_dump (b);
      printf ("c="); mpfr_dump (c);
      printf ("d="); mpfr_dump (d);
      printf ("Expected\n  ");
      mpfr_dump (ref);
      printf ("  with inex = %d\n", inex_ref);
      printf ("Got\n  ");
      mpfr_dump (res);
      printf ("  with inex = %d\n", inex_res);
      exit (1);
    }

  /* then check fmms */
  inex_res = mpfr_fmms (res, a, b, c, d, rnd);
  mpfr_mul (ab, a, b, rnd);
  mpfr_mul (cd, c, d, rnd);
  inex_ref = mpfr_sub (ref, ab, cd, rnd);
  if (! SAME_SIGN (inex_res, inex_ref) ||
      mpfr_nan_p (res) || mpfr_nan_p (ref) ||
      ! mpfr_equal_p (res, ref))
    {
      printf ("mpfr_fmms failed for p=%ld rnd=%s\n",
              (long int) p, mpfr_print_rnd_mode (rnd));
      printf ("a="); mpfr_dump (a);
      printf ("b="); mpfr_dump (b);
      printf ("c="); mpfr_dump (c);
      printf ("d="); mpfr_dump (d);
      printf ("Expected\n  ");
      mpfr_dump (ref);
      printf ("  with inex = %d\n", inex_ref);
      printf ("Got\n  ");
      mpfr_dump (res);
      printf ("  with inex = %d\n", inex_res);
      exit (1);
    }

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

static void
max_tests (void)
{
  mpfr_exp_t emin, emax;
  mpfr_t x, y1, y2;
  int r;
  int inex1, inex2;
  mpfr_flags_t flags1, flags2;

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();
  set_emin (MPFR_EMIN_MIN);
  set_emax (MPFR_EMAX_MAX);

  mpfr_inits2 (64, x, y1, y2, (mpfr_ptr) 0);
  mpfr_setmax (x, MPFR_EMAX_MAX);
  flags1 = MPFR_FLAGS_OVERFLOW | MPFR_FLAGS_INEXACT;
  RND_LOOP (r)
    {
      inex1 = mpfr_mul (y1, x, x, (mpfr_rnd_t) r);
      mpfr_clear_flags ();
      inex2 = mpfr_fmma (y2, x, x, x, x, (mpfr_rnd_t) r);
      flags2 = __gmpfr_flags;
      if (! (flags1 == flags2 && SAME_SIGN (inex1, inex2) &&
             mpfr_equal_p (y1, y2)))
        {
          printf ("Error in max_tests for %s",
                  mpfr_print_rnd_mode ((mpfr_rnd_t) r));
          printf ("Expected ");
          mpfr_dump (y1);
          printf ("  with inex = %d, flags =", inex1);
          flags_out (flags1);
          printf ("Got      ");
          mpfr_dump (y2);
          printf ("  with inex = %d, flags =", inex2);
          flags_out (flags2);
          exit (1);
        }
    }
  mpfr_clears (x, y1, y2, (mpfr_ptr) 0);

  set_emin (emin);
  set_emax (emax);
}

/* a^2 - (a+k)(a-k) = k^2 where a^2 overflows but k^2 usually doesn't.
 * With MPFR r10337 and gcc -fsanitize=undefined -fno-sanitize-recover,
 * this triggers the following error on a 64-bit machine:
 *   /usr/local/gmp-6.0.0-debug/include/gmp.h:2142:3:
 *   runtime error: store to null pointer of type 'mp_limb_t'
 * (mpn_add_1 is called on a null pointer).
 */
static void
near_overflow_tests (void)
{
  mpfr_exp_t old_emax;
  int i;

  old_emax = mpfr_get_emax ();

  for (i = 0; i < 40; i++)
    {
      mpfr_exp_t emax;
      mpfr_t a, k, c, d, z1, z2;
      mpfr_rnd_t rnd;
      mpfr_prec_t prec_a;
      int add = i & 1, inex, inex1, inex2;
      mpfr_flags_t flags1, flags2;

      /* In most cases, we do the test with the maximum exponent. */
      emax = i % 8 != 0 ? MPFR_EMAX_MAX : 64 + (randlimb () % 1);
      set_emax (emax);

      rnd = RND_RAND ();
      prec_a = 64 + randlimb () % 100;

      mpfr_init2 (a, prec_a);
      mpfr_urandom (a, RANDS, MPFR_RNDN);
      mpfr_set_exp (a, emax/2 + 32);

      mpfr_init2 (k, prec_a - 32);
      mpfr_urandom (k, RANDS, MPFR_RNDN);
      mpfr_set_exp (k, emax/2);

      mpfr_init2 (c, prec_a + 1);
      inex = mpfr_add (c, a, k, MPFR_RNDN);
      MPFR_ASSERTN (inex == 0);

      mpfr_init2 (d, prec_a);
      inex = mpfr_sub (d, a, k, MPFR_RNDN);
      MPFR_ASSERTN (inex == 0);
      if (add)
        mpfr_neg (d, d, MPFR_RNDN);

      mpfr_clear_flags ();
      inex1 = mpfr_sqr (z1, k, rnd);
      flags1 = __gmpfr_flags;

      mpfr_clear_flags ();
      inex2 = add ?
        mpfr_fmma (z2, a, a, c, d, rnd) :
        mpfr_fmms (z2, a, a, c, d, rnd);
      flags2 = __gmpfr_flags;

      if (! (flags1 == flags2 && SAME_SIGN (inex1, inex2) &&
             mpfr_equal_p (z1, z2)))
        {
          printf ("Error in near_overflow_tests for %s",
                  mpfr_print_rnd_mode (rnd));
          printf ("Expected ");
          mpfr_dump (z1);
          printf ("  with inex = %d, flags =", inex1);
          flags_out (flags1);
          printf ("Got      ");
          mpfr_dump (z2);
          printf ("  with inex = %d, flags =", inex2);
          flags_out (flags2);
          exit (1);
        }
    }

  set_emax (old_emax);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  random_tests ();
  zero_tests ();
  max_tests ();
  near_overflow_tests ();
  /* TODO: near_underflow_tests (); */

  tests_end_mpfr ();
  return 0;
}
