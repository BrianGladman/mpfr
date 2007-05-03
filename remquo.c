/* mpfr_remquo and mpfr_remainder -- argument reduction functions

Copyright 2007 Free Software Foundation, Inc.
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

#include <limits.h>  /* For CHAR_BIT */
#include "mpfr-impl.h"

/*
  Let q = x/y rounded to the nearest integer (to the nearest even number
  in case x/y = n + 1/2 with n integer).
  Put x - q*y in rem, rounded according to rnd.
  The value stored in *quo has the sign of q, and agrees with q with
  the 2^n low order bits. In other words, *quo = q (mod 2^n) and *quo q >= 0.
  If rem is zero, then it has the sign of x.
  The returned 'int' is the inexact flag giving the place of rem wrt x - q*y.

  If x or y is NaN: *quo is undefined, rem is NaN.
  If x is Inf, whatever y: *quo is undefined, rem is NaN.
  If y is Inf, x not NaN nor Inf: *quo is 0, rem is x.
  If y is 0, whatever x: *quo is undefined, rem is NaN.
  If x is 0, whatever y (not NaN nor 0): *quo is 0, rem is x.

  Otherwise if x and y are neither NaN, Inf nor 0, q is always defined,
  thus *quo is.
  Since |x - q*y| <= y/2, no overflow is possible.
  Only an underflow is possible when y is very small.
 */

/* compares the (absolute values of the) significands of x and y,
   returns a positive value if m(x) > m(y),
           zero if m(x) = m(y),
           a negative value of m(x) < m(y),
   where 1/2 <= m(x), m(y) < 1.
   Assumes x and y are neither NaN, nor Inf, nor zero. */
static int
mpfr_cmp_significand (mpfr_srcptr x, mpfr_srcptr y)
{
  mp_prec_t xn = MPFR_LIMB_SIZE (x), yn = MPFR_LIMB_SIZE (y);
  mp_ptr    xp = MPFR_MANT (x),      yp = MPFR_MANT (y);
  int i;

  if (xn <= yn)
    {
      i = mpn_cmp (xp, yp + (yn - xn), xn);
      if (i != 0)
        return i;
      /* if i = 0, either m(x) = m(y) if {yp, yn - xn} is zero,
         or m(x) < m(y) */
      for (yn = yn - xn; yn > 0 && yp[yn - 1] == 0; yn --);
      return (yn == 0) ? 0 : -1;
    }
  else /* xn > yn */
    {
      i = mpn_cmp (xp + (xn - yn), yp, yn);
      if (i != 0)
        return i;
      /* i = 0: need to check the low xn-yn limbs from x */
      for (xn = xn - yn; xn > 0 && xp[xn - 1] == 0; xn--);
      /* if xn > 0, one of the low limbs of x is non-zero */
      return (xn == 0) ? 0 : 1;
    }
}

/* we return as many bits as we can, keeping just one bit for the sign */
#define WANTED_BITS (sizeof(long) * CHAR_BIT - 1)

/* assuming q is an integer, with in addition EXP(q) <= PREC(q),
   returns |q| mod 2^WANTED_BITS */
static long
get_low_bits (mpfr_srcptr q)
{
  mp_ptr qp = MPFR_MANT(q);
  mp_size_t qn = MPFR_LIMB_SIZE(q);
  mp_size_t w = qn * BITS_PER_MP_LIMB - MPFR_EXP(q);
  long res;

  MPFR_ASSERTD(WANTED_BITS >= 3); /* to conform to C99 */
  MPFR_ASSERTN(WANTED_BITS < BITS_PER_MP_LIMB); /* required for this code
                                                   to work properly */
  /* weight of bit 0 of q is -w, with -w <= 0 since EXP(q) <= PREC(q),
     thus bit of weight 0 is w. */
  /* normally this loop should not be used, since in normal cases we have
     ulp(q)=1, i.e., EXP(q) = PREC(q), thus w is exactly the number of
     unused bits in the low significant limb of q. */
  while (w >= BITS_PER_MP_LIMB)
    {
      qp ++;
      qn --;
      w -= BITS_PER_MP_LIMB;
    }
  if ((w + WANTED_BITS <= BITS_PER_MP_LIMB) || (qn <= 1))
    /* all wanted bits in same limb */
    res = qp[0] >> w;
  else /* wanted bits are split in two consecutive limbs: BITS_PER_MP_LIMB-w
          in qp[0], and WANTED_BITS-(BITS_PER_MP_LIMB-w) in qp[1] */
    res = (qp[0] >> w) | (qp[1] << (BITS_PER_MP_LIMB - w));
  return res & ((MPFR_LIMB_ONE << WANTED_BITS) - MPFR_LIMB_ONE);
}

int
mpfr_remquo (mpfr_ptr rem, long *quo,
             mpfr_srcptr x, mpfr_srcptr y, mp_rnd_t rnd)
{
  mpfr_t q;
  int inex, compare;
  mp_exp_t ex, ey;
  mp_prec_t prec_q;

  MPFR_LOG_FUNC (("x[%#R]=%R y[%#R]=%R rnd=%d", x, x, y, y, rnd),
                 ("quo=%d rem[%#R]=%R", *quo, rem, rem));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x) || MPFR_IS_SINGULAR (y)))
    {
      if (MPFR_IS_NAN (x) || MPFR_IS_NAN (y) || MPFR_IS_INF (x)
          || MPFR_IS_ZERO (y))
        {
          MPFR_SET_NAN (rem);
          MPFR_RET_NAN;
        }
      else /* either y is Inf and x is 0 or non-special,
              or x is 0 and y is non-special,
              in both cases the quotient is zero. */
        {
          *quo = 0;
          return mpfr_set (rem, x, rnd);
        }
    }

  /* now neither x nor y is NaN, Inf or zero */

  /* first deal with the easy case where x is already reduced mod y,
     i.e., |x| <= |y|/2 */
  ex = MPFR_EXP (x);
  ey = MPFR_EXP (y);
  compare = mpfr_cmp_significand (x, y); /* compare > 0 if m(x) > m(y) */
  if ((ex + 1 < ey) || ((ex + 1 == ey) && compare <= 0))
    {
      *quo = 0;
      return mpfr_set (rem, x, rnd);
    }

  /* Now we are sure that the (true) quotient q is > 1/2 in absolute value.
     If ex = EXP(x) and ey = EXP(y), we have |x| < 2^ex and |y| >= 2^(ey-1),
     thus |q| < 2^(ex-ey+1).
     We also know that ex + 1 >= ey. */

  prec_q = ex - ey + (compare >= 0);
  mpfr_init2 (q, (prec_q < MPFR_PREC_MIN) ? MPFR_PREC_MIN : prec_q);
  inex = mpfr_div (q, x, y, GMP_RNDN);
  if (prec_q < MPFR_PREC_MIN)
    {
      int is_half_int, inex2;
      /* we might have a double-rounding problem in case q = n + 1/2 here,
         which can be obtained either from x/y = n + 1/2 + eps, or
         x/y = n + 1/2 - eps */
      mpfr_mul_2ui (q, q, 1, GMP_RNDN);
      is_half_int = mpfr_integer_p (q);
      mpfr_div_2ui (q, q, 1, GMP_RNDN);
      inex2 = mpfr_rint (q, q, GMP_RNDN);
      if (is_half_int && inex2 > 0 && inex > 0)
        /* mpfr_div rounded n + 1/2 - eps to n + 1/2, and mpfr_rint
           rounded n + 1/2 to n + 1 */
        mpfr_sub_ui (q, q, 1, GMP_RNDN);
      else if (is_half_int && inex2 < 0 && inex < 0)
        mpfr_add_ui (q, q, 1, GMP_RNDN);
    }

  /* set the low bits of the quotient */
  *quo = get_low_bits (q);
  /* quotient should have the sign of x/y */
  if (MPFR_SIGN(x) != MPFR_SIGN(y))
    *quo = -*quo;

  /* since we have no fused-multiply-and-subtract yet, we compute
     x + (-q)*y with an FMA */
  mpfr_neg (q, q, GMP_RNDN); /* exact */
  inex = mpfr_fma (rem, q, y, x, rnd);
  /* if rem is zero, it should have the sign of x */
  if (MPFR_IS_ZERO (rem))
      MPFR_SET_SAME_SIGN(rem, x);

  mpfr_clear (q);

  return inex;
}

int
mpfr_remainder (mpfr_ptr rem, mpfr_srcptr x, mpfr_srcptr y, mp_rnd_t rnd)
{
  long quo;

  return mpfr_remquo (rem, &quo, x, y, rnd);
}
