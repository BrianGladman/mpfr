/* Test file for mpfr_random_deviate

Copyright 2011-2014 Free Software Foundation, Inc.
Contributed by Charles Karney <charles@karney.com>, SRI International.

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

#ifndef WANT_MINI_GMP
#include "random_deviate.h"

#define W 32                    /* Must match value in random_deviate.c */

/* set random deviates rop from op */
static void
mpfr_random_deviate_set (mpfr_random_deviate_t rop, mpfr_random_deviate_t op)
{
  rop->e = op->e;
  rop->h = op->h;
  mpz_set (rop->f, op->f);
}

/* set random deviate to fract * 2^-expt.  expt must be a multiple
   of W and cannot be 0.  fract must be in [0,2^W) */
static void
mpfr_random_deviate_ldexp (mpfr_random_deviate_t rop,
			   unsigned long fract, unsigned long expt)
{
  rop->h = (expt > W ? 0ul : fract);
  mpz_set_ui(rop->f, expt > W ? fract: 0ul);
  rop->e = expt;
}

/* Test mpfr_random_deviate_less.  With two initially equal deviates
   this should return true half the time.  In order to execute
   additional code paths, the two deviates are repeatedly set equal and
   the test repeated (with now a longer fraction and with the test now
   triggering the sampling of an additional chunk. */
static void
test_compare (long nbtests, int verbose)
{
  mpfr_random_deviate_t u, v;
  mpfr_random_deviate_init (u);
  mpfr_random_deviate_init (v);
  int k, i, t1, t2;
  long count;

  count = 0;
  for (k = 0; k < nbtests; ++k) {
    mpfr_random_deviate_reset (u);
    mpfr_random_deviate_reset (v);
    for (i = 0; i < 10; ++i)
      {
        t1 = mpfr_random_deviate_less (u, v, RANDS);
        t2 = mpfr_random_deviate_less (u, v, RANDS);
        if (t1 != t2)
          {
            printf ("Error: mpfr_random_deviate_less() result inconsistent.\n");
            exit (1);
          }
        if (t1) ++count;
        /* Force the test to sample an additional chunk */
        mpfr_random_deviate_set(u, v);
      }
  }
  mpfr_random_deviate_clear (v);
  mpfr_random_deviate_clear (u);
#ifdef HAVE_STDARG
  if (verbose)
    printf ("Fraction of true random_deviate_less = %.4f\n",
	    count / (double)(10 * nbtests));
#endif /* HAVE_STDARG */
}

/* Test mpfr_random_deviate_value.  Check for the leading bit in the number in
   various positions. */
static void
test_value (long nbtests, mpfr_prec_t prec, mpfr_rnd_t rnd,
            int verbose)
{
  mpfr_t x;
  mpfr_random_deviate_t u;
  int inexact, inexactc;
  mpfr_random_deviate_init (u);
  mpfr_init2 (x, prec);
  int i, k, b, neg;
  unsigned long e, f, n;
  long count, sum;
  count = 0; sum = 0;
  inexactc = 1;

  for (k = 0; k < nbtests; ++k) {
    for (i = 0; i < 32; ++i) {
      b = gmp_urandomm_ui (RANDS, 32) + 1; /* bits to sample in integer */
      n = gmp_urandomb_ui (RANDS, b);
      neg = gmp_urandomb_ui (RANDS, 1);
      inexact = mpfr_random_deviate_value (neg, n, u, x, RANDS, rnd);
      inexactc *= inexact;
      if (inexact > 0) ++count;
      ++sum;
    }
    for (i = 0; i < 32; ++i) {
      b = gmp_urandomm_ui (RANDS, W) + 1; /* bits to sample in fraction */
      f = gmp_urandomb_ui (RANDS, b);
      e = W * (gmp_urandomm_ui (RANDS, 3) + 1);
      mpfr_random_deviate_ldexp (u, f, e);
      neg = gmp_urandomb_ui (RANDS, 1);
      inexact = mpfr_random_deviate_value (neg, 0, u, x, RANDS, rnd);
      inexactc *= inexact;
      if (inexact > 0) ++count;
      ++sum;
    }
    /* this tests the else clause for "if (x->e <= (mpfr_uexp_t)(-1) >> 1)" in
       mpfr_random_deviate_value.  However this immediately leads to a gmp:
       overflow in mpz type in "mpz_mul_2exp (qd, qd, x->e)".  The probability
       that this else branch is taken is very close to zero...  So skip this
       test for now. */
    /*
      e = ((((mpfr_uexp_t)(-1) >> 1) + W) / W) * W;
      f = gmp_urandomb_ui (RANDS, W);
      mpfr_random_deviate_ldexp (u, f, e);
      neg = gmp_urandomb_ui (RANDS, 1);
      inexact = mpfr_random_deviate_value (neg, 0, u, x, RANDS, rnd);
      inexactc *= inexact;
      if (inexact > 0) ++count;
      ++sum;
    */
    if (inexactc == 0)
      {
        printf ("Error: random_deviate() returns a zero ternary value.\n");
        exit (1);
      }
    mpfr_random_deviate_reset (u);
  }
  mpfr_random_deviate_clear (u);
  mpfr_clear(x);
#ifdef HAVE_STDARG
  if (verbose)
    printf ("Fraction of inexact > 0 = %.4f\n", count / (double)(sum));
#endif /* HAVE_STDARG */
}

int
main (int argc, char *argv[])
{
  long nbtests;
  int verbose;
  tests_start_mpfr ();

  verbose = 0;
  nbtests = 10;
  if (argc > 1)
    {
      long a = atol (argv[1]);
      verbose = 1;
      if (a != 0)
        nbtests = a;
    }

  test_compare (nbtests, verbose);
  test_value (nbtests,  2, MPFR_RNDD, verbose);
  test_value (nbtests,  5, MPFR_RNDU, verbose);
  test_value (nbtests, 24, MPFR_RNDN, verbose);
  test_value (nbtests, 53, MPFR_RNDZ, verbose);
  test_value (nbtests, 64, MPFR_RNDA, verbose);

  tests_end_mpfr ();
  return 0;
}

#else

int
main (void)
{
  return 77;
}

#endif
