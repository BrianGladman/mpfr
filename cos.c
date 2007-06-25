/* mpfr_cos -- cosine of a floating-point number

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <limits.h>

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* f <- 1 - r/2! + r^2/4! + ... + (-1)^l r^l/(2l)! + ...
   Assumes |r| < 1/2, and f, r have the same precision.
   Returns e such that the error on f is bounded by 2^e ulps.
*/
static int
mpfr_cos2_aux (mpfr_ptr f, mpfr_srcptr r)
{
  mpz_t x, t, s;
  mp_exp_t ex, l, m;
  mp_prec_t p, q;
  unsigned long i, maxi, imax;

  MPFR_ASSERTD(mpfr_get_exp (r) <= -1);

  /* compute minimal i such that i*(i+1) does not fit in an unsigned long,
     assuming that there are no padding bits. */
  maxi = 1UL << (CHAR_BIT * sizeof(unsigned long) / 2);
  if (maxi * (maxi / 2) == 0) /* test checked at compile time */
    {
      /* can occur only when there are padding bits. */
      /* maxi * (maxi-1) is representable iff maxi * (maxi / 2) != 0 */
      do
        maxi /= 2;
      while (maxi * (maxi / 2) == 0);
    }

  mpz_init (x);
  mpz_init (s);
  mpz_init (t);
  ex = mpfr_get_z_exp (x, r); /* r = x*2^ex */

  /* remove trailing zeroes */
  l = mpz_scan1 (x, 0);
  ex += l;
  mpz_div_2exp (x, x, l);

  /* since |r| < 1, r = x*2^ex, and x is an integer, necessarily ex < 0 */

  p = mpfr_get_prec (f); /* same than r */
  /* bound for number of iterations */
  imax = p / (-mpfr_get_exp (r));
  imax += (imax == 0);
  q = 2 * MPFR_INT_CEIL_LOG2(imax) + 4; /* bound for (3l)^2 */

  mpz_set_ui (s, 1); /* initialize sum with 1 */
  mpz_mul_2exp (s, s, p + q); /* scale all values by 2^(p+q) */
  mpz_set (t, s); /* invariant: t is previous term */
  for (i = 1; (m = mpz_sizeinbase (t, 2)) >= q; i += 2)
    {
      /* adjust precision of x to that of t */
      l = mpz_sizeinbase (x, 2);
      if (l > m)
        {
          l -= m;
          mpz_div_2exp (x, x, l);
          ex += l;
        }
      /* multiply t by r */
      mpz_mul (t, t, x);
      mpz_div_2exp (t, t, -ex);
      /* divide t by i*(i+1) */
      if (i < maxi)
        mpz_div_ui (t, t, i * (i + 1));
      else
        {
          mpz_div_ui (t, t, i);
          mpz_div_ui (t, t, i + 1);
        }
      /* if m is the (current) number of bits of t, we can consider that
         all operations on t so far had precision >= m, so we can prove
         by induction that the relative error on t is of the form
         (1+u)^(3l)-1, where |u| <= 2^(-m), and l=(i+1)/2 is the # of loops.
         Since |(1+x^2)^(1/x) - 1| <= 4x/3 for |x| <= 1/2,
         for |u| <= 1/(3l)^2, the absolute error is bounded by
         4/3*(3l)*2^(-m)*t <= 4*l since |t| < 2^m.
         Therefore the error on s is bounded by 2*l*(l+1). */
      /* add or subtract to s */
      if (i % 4 == 1)
        mpz_sub (s, s, t);
      else
        mpz_add (s, s, t);
    }

  mpfr_set_z (f, s, GMP_RNDN);
  mpfr_div_2ui (f, f, p + q, GMP_RNDN);

  mpz_clear (x);
  mpz_clear (s);
  mpz_clear (t);

  l = (i - 1) / 2; /* number of iterations */
  return 2 * MPFR_INT_CEIL_LOG2 (l + 1) + 1; /* bound is 2l(l+1) */
}

int
mpfr_cos (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mp_prec_t K, precy, m, k, l, precx;
  int inexact, reduce = 0;
  mpfr_t r, s, xr, c;
  mp_exp_t exps, cancel = 0, expx;
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_GROUP_DECL (group);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", x, x, rnd_mode),
                 ("y[%#R]=%R inexact=%d", y, y, inexact));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x) || MPFR_IS_INF (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else
        {
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
          return mpfr_set_ui (y, 1, rnd_mode);
        }
    }

  MPFR_SAVE_EXPO_MARK (expo);

  /* cos(x) = 1-x^2/2 + ..., so error < 2^(2*EXP(x)-1) */
  expx = MPFR_GET_EXP (x);
  MPFR_SMALL_INPUT_AFTER_SAVE_EXPO (y, __gmpfr_one, -2 * expx,
                                    1, 0, rnd_mode, expo, {});

  /* Compute initial precision */
  precy = MPFR_PREC (y);
  precx = MPFR_PREC (x);
  m = precy + 2 * MPFR_INT_CEIL_LOG2 (precy);

  if (expx >= 3)
    {
      reduce = 1;
      mpfr_init2 (xr, m);
      mpfr_init2 (c, expx + m - 1);
    }

  MPFR_GROUP_INIT_2 (group, m, r, s);
  MPFR_ZIV_INIT (loop, m);
  for (;;)
    {
      /* If |x| >= 4, first reduce x cmod (2*Pi) into xr, using mpfr_remainder:
	 let e = EXP(x) >= 3, and m the target precision:
	 (1) c <- 2*Pi              [precision e+m-1, nearest]
	 (2) xr <- remainder (x, c) [precision m, nearest]
	 We have |c - 2*Pi| <= 1/2ulp(c) = 2^(3-e-m)
	         |xr - x - k c| <= 1/2ulp(xr) <= 2^(1-m)
		 |k| <= |x|/(2*Pi) <= 2^(e-2)
         Thus |xr - x - 2kPi| <= |k| |c - 2Pi| + 2^(1-m) <= 2^(2-m).
	 It follows |cos(xr) - cos(x)| <= 2^(2-m). */
      if (reduce)
	{
	  mpfr_const_pi (c, GMP_RNDN);
	  mpfr_mul_2ui (c, c, 1, GMP_RNDN); /* 2Pi */
	  mpfr_remainder (xr, x, c, GMP_RNDN);
	  /* now |xr| <= 4, thus r <= 16 below */
	  mpfr_mul (r, xr, xr, GMP_RNDU); /* err <= 1 ulp */
	}
      else
	mpfr_mul (r, x, x, GMP_RNDU); /* err <= 1 ulp */

      /* we need |r| < 1/2 for mpfr_cos2_aux, i.e., EXP(r) - 2K <= -1 */
      K = MAX (MPFR_GET_EXP (r) + 1, 0);
      K = (K + 1) >> 1; /* ceil(K/2) */

      MPFR_SET_EXP (r, MPFR_GET_EXP (r) - 2 * K); /* Can't overflow! */

      /* s <- 1 - r/2! + ... + (-1)^l r^l/(2l)! */
      l = mpfr_cos2_aux (s, r);
      /* l is the error bound in ulps on s */
      MPFR_SET_ONE (r);
      for (k = 0; k < K; k++)
        {
          mpfr_sqr (s, s, GMP_RNDU);            /* err <= 2*olderr */
          MPFR_SET_EXP (s, MPFR_GET_EXP (s)+1); /* Can't overflow */
          mpfr_sub (s, s, r, GMP_RNDN);         /* err <= 4*olderr */
          MPFR_ASSERTD (MPFR_GET_EXP (s) <= 1);
        }

      /* The absolute error on s is bounded by (2l+1/3)*2^(2K-m)
         2l+1/3 <= 2l+1.
	 If |x| >= 4, we need to add 2^(2-m) for the argument reduction
	 by 2Pi: if K = 0, this amounts to add 4 to 2l+1/3, i.e., to add
	 2 to l; if K >= 1, this amounts to add 1 to 2*l+1/3. */
      l = 2 * l + 1;
      if (reduce)
	l += (K == 0) ? 4 : 1;
      k = MPFR_INT_CEIL_LOG2 (l) + 2*K;
      /* now the error is bounded by 2^(k-m) = 2^(EXP(s)-err) */

      exps = MPFR_GET_EXP (s);
      if (MPFR_LIKELY (MPFR_CAN_ROUND (s, exps + m - k, precy, rnd_mode)))
        break;

      if (MPFR_UNLIKELY (exps == 1))
        /* s = 1 or -1, and except x=0 which was already checked above,
	   cos(x) cannot be 1 or -1, so we can round if the error is less
	   than 2^(-precy) for directed rounding, or 2^(-precy-1) for rounding
	   to nearest. */
        {
          if (m > k && (m - k >= precy + (rnd_mode == GMP_RNDN)))
	    {
              /* if round to nearest or away, result is s,
                 otherwise it is round(nexttoward (s, 0)) */
	      if (MPFR_IS_LIKE_RNDZ (rnd_mode, MPFR_IS_NEG (s)))
		mpfr_nexttozero (s);
	      break;
	    }
        }

      if (exps < cancel)
        {
          m += cancel - exps;
          cancel = exps;
        }

      MPFR_ZIV_NEXT (loop, m);
      MPFR_GROUP_REPREC_2 (group, m, r, s);
      if (reduce)
	{
	  mpfr_set_prec (xr, m);
	  mpfr_set_prec (c, expx + m - 1);
	}
    }
  MPFR_ZIV_FREE (loop);
  inexact = mpfr_set (y, s, rnd_mode);
  MPFR_GROUP_CLEAR (group);
  if (reduce)
    {
      mpfr_clear (xr);
      mpfr_clear (c);
    }

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}
