/* mpfr_zeta -- compute the Riemann Zeta function

Copyright 2003, 2004, 2005 Free Software Foundation.
Contributed by Jean-Luc Re'my and the Spaces project, INRIA Lorraine.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

/* #define DEBUG */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/*
   Parameters:
   s - the input floating-point number
   n, p - parameters from the algorithm
   tc - an array of p floating-point numbers tc[1]..tc[p]
   Output:
   b is the result, i.e.
   sum(tc[i]*product((s+2j)*(s+2j-1)/n^2,j=1..i-1), i=1..p)*s*n^(-s-1)
*/
static void
mpfr_zeta_part_b (mpfr_t b, mpfr_srcptr s, int n, int p, mpfr_t *tc)
{
  int n2, l, t, precb;
  mpfr_t s1, d, u;

  if (p == 0)
    {
      mpfr_set_ui (b, 0, GMP_RNDN);
      return;
    }

  n2 = n * n;
  precb = mpfr_get_prec (b);
  mpfr_init2 (s1, precb);
  mpfr_init2 (d, precb);
  mpfr_init2 (u, precb);
  /* t equals 2p-2, 2p-3, ... ; s1 equals s+t */
  t = 2 * p - 2;
  mpfr_set (d, tc[p], GMP_RNDN);
  for (l = 1; l < p; l++)
    {
      mpfr_add_ui (s1, s, t, GMP_RNDN); /* s + (2p-2l) */
      mpfr_mul (d, d, s1, GMP_RNDN);
      t = t - 1;
      mpfr_add_ui (s1, s, t, GMP_RNDN); /* s + (2p-2l-1) */
      mpfr_mul (d, d, s1, GMP_RNDN);
      t = t - 1;
      mpfr_div_ui (d, d, n2, GMP_RNDN);
      mpfr_add (d, d, tc[p-l], GMP_RNDN);
      /* since s is positive and the tc[i] have alternate signs,
         the following is unlikely */
      if (MPFR_UNLIKELY(mpfr_cmpabs (d, tc[p-l]) > 0))
	mpfr_set (d, tc[p-l], GMP_RNDN);
    }
  mpfr_mul (d, d, s, GMP_RNDN);
  mpfr_add_ui (s1, s, 1, GMP_RNDN);
  mpfr_neg (s1, s1, GMP_RNDN);
  mpfr_ui_pow (u, n, s1, GMP_RNDN);
  mpfr_mul (b, d, u, GMP_RNDN);
  MPFR_TRACE (MPFR_DUMP (b));
  mpfr_clear (s1);
  mpfr_clear (d);
  mpfr_clear (u);
}

/* Input: p - an integer
   Output: fills tc[1..p], tc[i] = bernoulli(2i)/(2i)!
   tc[1]=1/12, tc[2]=-1/720, tc[3]=1/30240, ...
*/
static void
mpfr_zeta_c (int p, mpfr_t *tc)
{
  mpfr_t d;
  int k, l;

  if (p > 0)
    {
      mpfr_init2 (d, mpfr_get_prec (tc[1]));
      mpfr_set_ui (tc[1], 1, GMP_RNDN);
      mpfr_div_ui (tc[1], tc[1], 12, GMP_RNDN);
      for (k = 2; k <= p; k++)
	{
	  mpfr_set_ui (d, k-1, GMP_RNDN);
	  mpfr_div_ui (d, d, 12*k+6, GMP_RNDN);
	  for (l=2; l<=k-1; l++)
	    {
	      mpfr_div_ui (d, d, 4*(2*k-2*l+3)*(2*k-2*l+2), GMP_RNDN);
	      mpfr_add (d, d, tc[l], GMP_RNDN);
	    }
	  mpfr_div_ui (tc[k], d, 24, GMP_RNDN);
	  mpfr_neg (tc[k], tc[k], GMP_RNDN);
	}
      mpfr_clear(d);
    }
}

/* Input: s - a floating-point number
          n - an integer
   Output: sum - a floating-point number approximating sum(1/i^s, i=1..n-1) */
static void
mpfr_zeta_part_a (mpfr_t sum, mpfr_srcptr s, int n)
{
  int i, preca;
  mpfr_t u, s1;

  preca = MPFR_PREC (sum);
  mpfr_init2 (u, preca);
  mpfr_init2 (s1, preca);
  mpfr_neg (s1, s, GMP_RNDN);
  mpfr_ui_pow (u, n, s1, GMP_RNDN);
  mpfr_div_2exp (u, u, 1, GMP_RNDN);
  mpfr_set (sum, u, GMP_RNDN);
  for (i=n-1; i>1; i--)
    {
      mpfr_ui_pow (u, i, s1, GMP_RNDN);
      mpfr_add (sum, sum, u, GMP_RNDN);
    }
  mpfr_add_ui (sum, sum, 1, GMP_RNDN);
  MPFR_TRACE (MPFR_DUMP (sum));
  mpfr_clear (s1);
  mpfr_clear (u);
}

/* Input: s - a floating-point number >= 1/2.
          rnd_mode - a rounding mode.
          Assumes s is neither NaN nor Infinite.
   Output: z - Zeta(s) rounded to the precision of z with direction rnd_mode
*/
static int
mpfr_zeta_pos (mpfr_t z, mpfr_srcptr s, mp_rnd_t rnd_mode)
{
  int p, n, l, add;
  double beta, sd, dnep;
  mpfr_t a, b, c, z_pre, f, g, s1;
  mpfr_t *tc1;
  mp_prec_t precz, precs, d, dint;
  int inex;
  MPFR_ZIV_DECL (loop);

  MPFR_TRACE (MPFR_DUMP (s));

  precz = MPFR_PREC (z);
  precs = MPFR_PREC (s);

  d = precz + 11;
  /* we want that s1 = s-1 is exact, i.e. we should have PREC(s1) >= EXP(s) */
  mpfr_init2 (s1, MAX (precs, ((mpfr_uexp_t) MPFR_EXP (s))));

  mpfr_init2 (a, MPFR_PREC_MIN);
  mpfr_init2 (b, MPFR_PREC_MIN);
  mpfr_init2 (c, MPFR_PREC_MIN);
  mpfr_init2 (z_pre, MPFR_PREC_MIN);
  mpfr_init2 (f, MPFR_PREC_MIN);
  mpfr_init2 (g, MPFR_PREC_MIN);

  MPFR_ZIV_INIT (loop, d);
  for (;;)
    {
      /* Principal loop: we compute, in z_pre,
	 an approximation of Zeta(s), that we send to mpfr_can_round */
      mpfr_sub_ui (s1, s, 1, GMP_RNDN);
      MPFR_ASSERTN (MPFR_IS_FP (s1));

      if (MPFR_IS_ZERO (s1))
        {
          mpfr_set_inf (z, 1);
          inex = 0;
          goto clear_and_return;
        }
      else if (MPFR_GET_EXP (s1) <= -(mp_exp_t) ((mpfr_prec_t) (d-3)/2))
	/* Branch 1: when s-1 is very small, one
	  uses the approximation Zeta(s)=1/(s-1)+gamma,
	  where gamma is Euler's constant */
	{
	  dint = MAX (d + 3, precs);
	  MPFR_TRACE (printf ("branch 1\ninternal precision=%d\n", dint));
	  mpfr_set_prec (z_pre, dint);
	  mpfr_set_prec (g, dint);
	  mpfr_ui_div (z_pre, 1, s1, GMP_RNDN);
	  mpfr_const_euler (g, GMP_RNDN);
	  mpfr_add (z_pre, z_pre, g, GMP_RNDN);
	}
      else /* Branch 2 */
        {
          size_t size;

          MPFR_TRACE (printf ("branch 2\n"));
          /* Computation of parameters n, p and working precision */
          dnep = (double) d * LOG2;
          sd = mpfr_get_d (s, GMP_RNDN);
          /* beta = dnep + 0.61 + sd * log (6.2832 / sd);
             but a larger value is ok */
#define LOG6dot2832 1.83787940484160805532
          beta = dnep + 0.61 + sd * (LOG6dot2832 - LOG2 *
                                     __gmpfr_floor_log2 (sd));
          if (beta <= 0.0)
            {
              p = 0;
              /* n = 1 + (int) (exp ((dnep - LOG2) / sd)); */
              n = 1 + (int) __gmpfr_ceil_exp2 ((d - 1.0) / sd);
            }
          else
            {
              p = 1 + (int) beta / 2;
              n = 1 + (int) ((sd + 2.0 * (double) p - 1.0) / 6.2832);
            }
          MPFR_TRACE (printf ("\nn=%d\np=%d\n",n,p));
          /* add = 4 + floor(1.5 * log(d) / log (2)).
             We should have add >= 10, which is always fulfilled since
             d = precz + 11 >= 12, thus ceil(log2(d)) >= 4 */
          add = 4 + (3 * MPFR_INT_CEIL_LOG2 (d)) / 2;
          MPFR_ASSERTD(add >= 10);
          dint = d + add;
          if (dint < precs)
            dint = precs;

          MPFR_TRACE (printf("internal precision=%d\n",dint));

          size = (p + 1) * sizeof(mpfr_t);
          tc1 = (mpfr_t*) (*__gmp_allocate_func) (size);
          for (l=1; l<=p; l++)
            mpfr_init2 (tc1[l], dint);
          mpfr_set_prec (a, dint);
          mpfr_set_prec (b, dint);
          mpfr_set_prec (c, dint);
          mpfr_set_prec (z_pre, dint);
          mpfr_set_prec (f, dint);

          MPFR_TRACE (printf ("precision of z =%d\n", precz));

          /* Computation of the coefficients c_k */
          mpfr_zeta_c (p, tc1);
          /* Computation of the 3 parts of the fonction Zeta. */
          mpfr_zeta_part_a (a, s, n);
          mpfr_zeta_part_b (b, s, n, p, tc1);
          /* s1 = s-1 is already computed above */
          /* mpfr_sub_ui (s1, s, 1, GMP_RNDN); */
          mpfr_ui_div (c, 1, s1, GMP_RNDN);
          mpfr_ui_pow (f, n, s1, GMP_RNDN);
          mpfr_div (c, c, f, GMP_RNDN);
	  MPFR_TRACE (MPFR_DUMP (c));
          mpfr_add (z_pre, a, c, GMP_RNDN);
          mpfr_add (z_pre, z_pre, b, GMP_RNDN);
          for (l=1; l<=p; l++)
            mpfr_clear (tc1[l]);
          (*__gmp_free_func) (tc1, size);
          /* End branch 2 */
        }

      MPFR_TRACE (MPFR_DUMP (z_pre));
      if (MPFR_LIKELY (mpfr_can_round (z_pre, d - 3, GMP_RNDN, GMP_RNDZ,
				       precz + (rnd_mode == GMP_RNDN))))
	break;
      MPFR_ZIV_NEXT (loop, d);
    }
  MPFR_ZIV_FREE (loop);

  inex = mpfr_set (z, z_pre, rnd_mode);
  MPFR_TRACE (MPFR_DUMP (z));

 clear_and_return:
  mpfr_clear (a);
  mpfr_clear (b);
  mpfr_clear (c);
  mpfr_clear (z_pre);
  mpfr_clear (f);
  mpfr_clear (g);
  mpfr_clear (s1);

  return inex;
}

int
mpfr_zeta (mpfr_t z, mpfr_srcptr s, mp_rnd_t rnd_mode)
{
  double sd, eps, m1, c;
  long add;
  mpfr_t z_pre, s1, s2, y, p;
  mp_prec_t precz, prec1, precs, precs1;
  int inex;
  MPFR_SAVE_EXPO_DECL (expo);

  /* Zero, Nan or Inf ? */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (s)))
    {
      if (MPFR_IS_NAN (s))
	{
	  MPFR_SET_NAN (z);
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_INF (s))
	{
	  if (MPFR_IS_POS (s))
	    return mpfr_set_ui (z, 1, GMP_RNDN); /* Zeta(+Inf) = 1 */
	  MPFR_SET_NAN (z); /* Zeta(-Inf) = NaN */
	  MPFR_RET_NAN;
	}
      else /* s iz zero */
	{
          MPFR_ASSERTD (MPFR_IS_ZERO (s));
	  mpfr_set_ui (z, 1, rnd_mode);
	  mpfr_div_2ui (z, z, 1, rnd_mode);
	  MPFR_CHANGE_SIGN (z);
	  MPFR_RET (0);
	}
    }
  MPFR_CLEAR_FLAGS(z);

  /* s is neither Nan, nor Inf, nor Zero */
  mpfr_init2 (s2, MPFR_PREC (s));
  mpfr_div_2ui (s2, s, 1, rnd_mode);
  if (MPFR_IS_NEG (s) && mpfr_floor (s2, s2) == 0) /* Case s = -2n */
    {
      mpfr_clear (s2);
      return mpfr_set_ui (z, 0, rnd_mode);
    }
  mpfr_clear (s2);

  MPFR_SAVE_EXPO_MARK (expo);

  /* Compute Zeta */
  if (MPFR_IS_POS (s) && MPFR_GET_EXP (s) >= 0) /* Case s >= 1/2 */
    inex = mpfr_zeta_pos (z, s, rnd_mode);
  else /* use reflection formula
          zeta(s) = 2^s*Pi^(s-1)*sin(Pi*s/2)*gamma(1-s)*zeta(1-s) */
    {
      MPFR_ZIV_DECL (loop);

      precz = MPFR_PREC (z);
      precs = MPFR_PREC (s);

      /* Precision precs1 needed to represent 1 - s, and s + 2,
	 without any truncation */
      precs1 = precs + 2 + MAX (0, - MPFR_GET_EXP (s));
      sd = mpfr_get_d (s, GMP_RNDN) - 1.0;
      if (sd < 0.0)
	sd = -sd; /* now sd = abs(s-1.0) */
      /* Precision prec1 is the precision on elementary computations;
	 it ensures a final precision prec1 - add for zeta(s) */
      /* eps = pow (2.0, - (double) precz - 14.0); */
      eps = __gmpfr_ceil_exp2 (- (double) precz - 14.0);
      m1 = 1.0 + MAX(1.0 / eps,  2.0 * sd) * (1.0 + eps);
      c = (1.0 + eps) * (1.0 + eps * MAX(8.0, m1));
      /* add = 1 + floor(log(c*c*c*(13 + m1))/log(2)); */
      add = __gmpfr_ceil_log2 (c * c * c * (13.0 + m1));
      prec1 = precz + add; 
      prec1 = MAX (prec1, precs1) + 10;
      
      mpfr_init2 (z_pre, prec1);
      mpfr_init2 (s1, prec1);
      mpfr_init2 (y, prec1);
      mpfr_init2 (p, prec1);
      
      MPFR_ZIV_INIT (loop, prec1);
      for (;;)
	{
          mpfr_ui_sub (s1, 1, s, GMP_RNDN); /* s1 = 1-s */
          mpfr_zeta_pos (z_pre, s1, GMP_RNDN); /* zeta(1-s)  */
          mpfr_gamma (y, s1, GMP_RNDN);        /* gamma(1-s) */
          mpfr_mul (z_pre, z_pre, y, GMP_RNDN); /* gamma(1-s)*zeta(1-s) */
          mpfr_const_pi (p, GMP_RNDD);
          mpfr_mul (y, s, p, GMP_RNDN);
          mpfr_div_2ui (y, y, 1, GMP_RNDN); /* s*Pi/2 */
          mpfr_sin (y, y, GMP_RNDN); /* sin(Pi*s/2) */
          mpfr_mul (z_pre, z_pre, y, GMP_RNDN);
          mpfr_mul_2ui (y, p, 1, GMP_RNDN); /* 2*Pi */
          mpfr_neg (s1, s1, GMP_RNDN); /* s-1 */
          mpfr_pow (y, y, s1, GMP_RNDN); /* (2*Pi)^(s-1) */
          mpfr_mul (z_pre, z_pre, y, GMP_RNDN);
          mpfr_mul_2ui (z_pre, z_pre, 1, GMP_RNDN);

          if (MPFR_LIKELY (mpfr_can_round (z_pre, prec1 - add, GMP_RNDN, 
					   GMP_RNDZ,
					   precz + (rnd_mode == GMP_RNDN))))
	    break;

	  /* Actualisation of the precision */
	  MPFR_ZIV_NEXT (loop, prec1);
	  mpfr_set_prec (z_pre, prec1);
	  mpfr_set_prec (s1, prec1);
	  mpfr_set_prec (y, prec1);
	  mpfr_set_prec (p, prec1);
        }
      MPFR_ZIV_FREE (loop);

      inex = mpfr_set (z, z_pre, rnd_mode);

      mpfr_clear(z_pre);
      mpfr_clear(s1);
      mpfr_clear(y);
      mpfr_clear(p);
    }

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (z, inex, rnd_mode);
}
