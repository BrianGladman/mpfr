/* Chi-squared test for mpfr_erandom

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

/* Return Phi(x) = 1 - exp(-x), the cumulative probability function for the
 * exponential distribution.  We only take differences of this function so the
 * offset doesn't matter; here Phi(0) = 0. */
static void
exponential_cumulative (mpfr_t z, mpfr_t x, mpfr_rnd_t rnd)
{
  mpfr_neg (z, x, rnd);
  mpfr_expm1 (z, z, rnd);
  mpfr_neg (z, z, rnd);
}

/* The continuous chi-squared test on with a set of bins of equal width.
 *
 * A single precision is picked for sampling and the chi-squared calculation.
 * This should picked high enough so that binning in test doesn't need to be
 * accurately aligned with possible values of the deviates.  Also we need the
 * precision big enough that chi-squared calculation itself is reliable.
 *
 * There's no particular benefit is testing with at very higher precisions;
 * because of the way terandom samples, this just adds additional barely
 * significant random bits to the deviates.  So this chi-squared test with
 * continuous equal width bins isn't a good tool for finding problems here.
 *
 * The testing of low precision exponential deviates is done by
 * test_erandom_chisq_disc. */
static void
test_erandom_chisq_cont (long num, mpfr_prec_t prec, int nu,
                         double xmin, double xmax, int verbose)
{
  mpfr_t x, a, b, dx, z, pa, pb, ps, t;
  long *counts;
  int i, inexact;
  long k;
  mpfr_rnd_t rnd, rndd;
  double chilim, xp, sqrt2dof;

  rnd = MPFR_RNDN;              /* For chi-squared calculation */
  rndd = MPFR_RNDD;             /* For sampling and figuring the bins */
  mpfr_inits2 (prec, x, a, b, dx, z, pa, pb, ps, t, (mpfr_ptr) 0);

  counts = (long *) calloc (nu + 1, sizeof (long));
  if (counts == NULL)
    {
      fprintf (stderr, "terandom_chisq: can't allocate memory\n");
      exit (1);
    }

  /* a and b are bounds of nu equally spaced bins.  Set dx = (b-a)/nu */
  mpfr_set_d (a, xmin, rnd);
  mpfr_set_d (b, xmax, rnd);

  mpfr_sub (dx, b, a, rnd);
  mpfr_div_si (dx, dx, nu, rnd);

  for (k = 0; k < num; ++k)
    {
      inexact = mpfr_erandom (x, RANDS, rndd);
      if (inexact == 0)
        {
          /* one call in the loop pretended to return an exact number! */
          printf ("Error: mpfr_erandom() returns a zero ternary value.\n");
          exit (1);
        }
      if (mpfr_signbit (x))
        {
          printf ("Error: mpfr_erandom() returns a negative deviate.\n");
          exit (1);
        }
      mpfr_sub (x, x, a, rndd);
      mpfr_div (x, x, dx, rndd);
      i = mpfr_get_si (x, rndd);
      ++counts[i >= 0 && i < nu ? i : nu];
    }

  mpfr_set (x, a, rnd);
  exponential_cumulative (pa, x, rnd);
  mpfr_add_ui (ps, pa, 1, rnd);
  mpfr_set_zero (t, 1);
  for (i = 0; i <= nu; ++i)
    {
      if (i < nu)
        {
          mpfr_add (x, x, dx, rnd);
          exponential_cumulative (pb, x, rnd);
          mpfr_sub (pa, pb, pa, rnd); /* prob for this bin */
        }
      else
        mpfr_sub (pa, ps, pa, rnd); /* prob for last bin, i = nu */

      /* Compute z = counts[i] - num * p; t += z * z / (num * p) */
      mpfr_mul_ui (pa, pa, num, rnd);
      mpfr_ui_sub (z, counts[i], pa, rnd);
      mpfr_sqr (z, z, rnd);
      mpfr_div (z, z, pa, rnd);
      mpfr_add (t, t, z, rnd);
      mpfr_swap (pa, pb);       /* i.e., pa = pb */
    }

  if (verbose)
    mpfr_printf ("num = %ld, discrete (prec = %ld) bins in [%.2f, %.2f], "
                 "nu = %d: chi2 = %.2Rf\n", num, prec, xmin, xmax, nu, t);

  /* See Knuth, TAOCP, Vol 2, 3.3.1, Table 1, approx formula for nu > 30; xp =
   * +/- 3.1 gives 0.1% and 99.9% percentile points. */
  mpfr_set_si (pa, 2*nu, rnd);
  mpfr_sqrt (pa, pa, rnd);
  sqrt2dof = mpfr_get_d (pa, rnd);
  xp = 3.1;
  chilim = nu - sqrt2dof * xp  +   2 * xp * xp / 3;
  if (mpfr_cmp_d (t, chilim) < 0)
    {
      mpfr_printf ("chi2 = %.2Rf is less than %.2f"
                   " should only happen 1 in 1000 times\n", t, chilim);
      exit (1);
    }
  chilim = nu + sqrt2dof * xp  +   2 * xp * xp / 3;
  if (mpfr_cmp_d (t, chilim) > 0)
    {
      mpfr_printf ("chi2 = %.2Rf is greater than %.2f"
                   " should only happen 1 in 1000 times\n", t, chilim);
      exit (1);
    }

  mpfr_clears (x, a, b, dx, z, pa, pb, ps, t, (mpfr_ptr) 0);
}

/* Return a sequential number for a positive low-precision x.  x is altered by
 * this fuction.  low precision means prec = 2, 3, or 4.  High values of
 * precision will result in integer overflow. */
static long
sequential (mpfr_t x)
{
  long expt, prec;

  prec = mpfr_get_prec (x);
  expt =  mpfr_get_exp (x);
  mpfr_mul_2si (x, x, prec - expt, MPFR_RNDN);

  return expt * (1 << (prec - 1)) + mpfr_get_si (x, MPFR_RNDN);
}

/* The chi-squared test on low precision exponential deviates.  wprec is the
 * working precision for the chi-squared calculation.  prec is the precision
 * for the sampling; choose this in [2,5].  The bins consist of all the
 * possible deviate values in the range [xmin, xmax] coupled with the value of
 * inexact.  Thus with prec = 2, the bins are
 *   ...
 *   (7/16, 1/2)  x = 1/2, inexact = +1
 *   (1/2 , 5/8)  x = 1/2, inexact = -1
 *   (5/8 , 3/4)  x = 3/4, inexact = +1
 *   (3/4 , 7/8)  x = 3/4, inexact = -1
 *   (7/8 , 1  )  x = 1  , inexact = +1
 *   (1   , 5/4)  x = 1  , inexact = -1
 *   (5/4 , 3/2)  x = 3/2, inexact = +1
 *   (3/2 , 7/4)  x = 3/2, inexact = -1
 *   ...
 * In addition, two bins are allocated for [0,xmin) and (xmax,inf).
 *
 * The sampling is with MPFR_RNDN.  This is the rounding mode which elicits the
 * most information.  trandom_deviate includes checks on the consistency of the
 * results extracted from a random_deviate with other rounding modes.  */
static void
test_erandom_chisq_disc (long num, mpfr_prec_t wprec, mpfr_prec_t prec,
                         double xmin, double xmax, int verbose)
{
  mpfr_t x, v, pa, pb, z, t;
  mpfr_rnd_t rnd;
  int i, inexact, nu;
  long *counts;
  long k, seqmin, seqmax, seq;
  double chilim, xp, sqrt2dof;

  rnd = MPFR_RNDN;
  mpfr_init2 (x, prec);
  mpfr_init2 (v, prec+1);
  mpfr_inits2 (wprec, pa, pb, z, t, (mpfr_ptr) 0);

  mpfr_set_d (x, xmin, rnd);
  xmin = mpfr_get_d (x, rnd);
  mpfr_set (v, x, rnd);
  seqmin = sequential (x);
  mpfr_set_d (x, xmax, rnd);
  xmax = mpfr_get_d (x, rnd);
  seqmax = sequential (x);

  /* Two bins for each sequential number (for inexact = +/- 1), plus 1 for u <
   * umin and 1 for u > umax, minus 1 for degrees of freedom */
  nu = 2 * (seqmax - seqmin + 1) + 2 - 1;
  counts = (long *) calloc (nu + 1, sizeof (long));
  if (counts == NULL)
    {
      fprintf (stderr, "terandom_chisq: can't allocate memory\n");
      exit (1);
    }

  for (k = 0; k < num; ++k)
    {
      inexact = mpfr_erandom (x, RANDS, rnd);
      if (mpfr_signbit (x))
        {
          printf ("Error: mpfr_erandom() returns a negative deviate.\n");
          exit (1);
        }
      /* Don't call sequential with small args to avoid undefined behavior with
       * zero and possibility of overflow. */
      seq = mpfr_greaterequal_p (x, v) ? sequential (x) : seqmin - 1;
      ++counts[seq < seqmin ? 0 :
               seq <= seqmax ? 2 * (seq - seqmin) + 1 + (inexact > 0 ? 0 : 1) :
               nu];
    }

  mpfr_set_zero (v, 1);
  exponential_cumulative (pa, v, rnd);
  /* Cycle through all the bin boundaries using mpfr_nextabove at precision
   * prec + 1 starting at mpfr_nextbelow (xmin) */
  mpfr_set_d (x, xmin, rnd);
  mpfr_set (v, x, rnd);
  mpfr_nextbelow (v);
  mpfr_nextbelow (v);
  mpfr_set_zero (t, 1);
  for (i = 0; i <= nu; ++i)
    {
       if (i < nu)
        mpfr_nextabove (v);
      else
        mpfr_set_inf (v, 1);
      exponential_cumulative (pb, v, rnd);
      mpfr_sub (pa, pb, pa, rnd);

      /* Compute z = counts[i] - num * p; t += z * z / (num * p). */
      mpfr_mul_ui (pa, pa, num, rnd);
      mpfr_ui_sub (z, counts[i], pa, rnd);
      mpfr_sqr (z, z, rnd);
      mpfr_div (z, z, pa, rnd);
      mpfr_add (t, t, z, rnd);
      mpfr_swap (pa, pb);       /* i.e., pa = pb */
    }
  if (verbose)
    mpfr_printf ("num = %ld, discrete (prec = %ld) bins in [%.6f, %.2f], "
                 "nu = %d: chi2 = %.2Rf\n", num, prec, xmin, xmax, nu, t);

  /* See Knuth, TAOCP, Vol 2, 3.3.1, Table 1, approx formula for nu > 30; xp =
   * +/- 3.1 gives 0.1% and 99.9% percentile points. */
  mpfr_set_si (pa, 2*nu, rnd);
  mpfr_sqrt (pa, pa, rnd);
  sqrt2dof = mpfr_get_d (pa, rnd);
  xp = 3.1;
  chilim = nu - sqrt2dof * xp  +   2 * xp * xp / 3;
  if (mpfr_cmp_d (t, chilim) < 0)
    {
      mpfr_printf ("chi2 = %.2Rf is less than %.2f"
                   " should only happen 1 in 1000 times\n", t, chilim);
      exit (1);
    }
  chilim = nu + sqrt2dof * xp  +   2 * xp * xp / 3;
  if (mpfr_cmp_d (t, chilim) > 0)
    {
      mpfr_printf ("chi2 = %.2Rf is greater than %.2f"
                   " should only happen 1 in 1000 times\n", t, chilim);
      exit (1);
    }

  free (counts);
  mpfr_clears (x, v, pa, pb, z, t, (mpfr_ptr) 0);
}

int
main (int argc, char *argv[])
{
  long nbtests;
  int verbose;
  tests_start_mpfr ();

  verbose = 0;
  nbtests = 100000;
  if (argc > 1)
    {
      long a = atol (argv[1]);
      verbose = 1;
      if (a != 0)
        nbtests = a;
    }

  test_erandom_chisq_cont (nbtests, 64, 60, 0, 7, verbose);
  test_erandom_chisq_disc (nbtests, 64, 2, 0.002, 6, verbose);
  test_erandom_chisq_disc (nbtests, 64, 3, 0.02, 7, verbose);
  test_erandom_chisq_disc (nbtests, 64, 4, 0.04, 8, verbose);

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
