/* mpfr_sin -- sine of a floating-point number

Copyright 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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

#include <stdio.h>

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* determine the sign of sin(x) using argument reduction.
   Assumes x is not an exact multiple of Pi (this excludes x=0). */
static int
mpfr_sin_sign (mpfr_srcptr x)
{
  mpfr_t c, k;
  mp_exp_t K;
  int sign;
  mp_prec_t m;
  mpfr_srcptr y;
  MPFR_ZIV_DECL (loop);

  K = MPFR_GET_EXP(x);
  if (K < 0)  /* Trivial case if ABS(x) < 1 */
    return MPFR_SIGN (x);

  m = K + BITS_PER_MP_LIMB;
  mpfr_init2 (c, m);
  mpfr_init2 (k, m);

  MPFR_ZIV_INIT (loop, m);
  for (;;)
    {
      /* first determine round(x/Pi): does not have to be exact since
         the result is an integer */
      mpfr_const_pi (c, GMP_RNDN); /* err <= 1/2*ulp(c) = 2^(1-m) */
      /* we need that k is not-to-badly rounded to an integer,
         i.e. ulp(k) <= 1, so m >= EXP(k). */
      mpfr_div (k, x, c, GMP_RNDN);
      mpfr_round (k, k);

      sign = 1;

      if (!MPFR_IS_ZERO (k)) /* subtract k*approx(Pi) */
        {
          /* determine parity of k for sign */
          if (MPFR_GET_EXP (k) <= 0 || (mpfr_uexp_t) MPFR_GET_EXP (k) <= m)
            {
              mp_size_t j = BITS_PER_MP_LIMB * MPFR_LIMB_SIZE(k) 
		- MPFR_GET_EXP(k);
              mp_size_t l = j / BITS_PER_MP_LIMB;
              /* parity bit is j-th bit starting from least significant bits */
              if ((MPFR_MANT(k)[l] >> (j % BITS_PER_MP_LIMB)) & 1)
                sign = -1; /* k is odd */
            }
          K = MPFR_GET_EXP (k); /* k is an integer, thus K >= 1, k < 2^K */
          mpfr_mul (k, k, c, GMP_RNDN); /* err <= oldk*err(c) + 1/2*ulp(k)
                                               <= 2^(K+2-m) */
          mpfr_sub (k, x, k, GMP_RNDN);
          /* assuming |k| <= Pi, err <= 2^(1-m)+2^(K+2-m) < 2^(K+3-m) */
	  MPFR_ASSERTN (MPFR_GET_EXP (k) <= 2); 
          y = k;
        }
      else
        {
          K = 1;
          y = x;
        }
      /* sign of sign(y) is uncertain if |y| <= err < 2^(K+3-m),
         thus EXP(y) < K+4-m */
      if (MPFR_LIKELY (!MPFR_IS_ZERO (y) 
		       && MPFR_GET_EXP (y) >= K + 4 - (mp_exp_t) m))
	break;
      MPFR_ZIV_NEXT (loop, m);
      mpfr_set_prec (c, m);
      mpfr_set_prec (k, m);
    }

  if (MPFR_IS_NEG (y))
    sign = -sign;

  mpfr_clear (k);
  mpfr_clear (c);
  
  return sign;
}

int 
mpfr_sin (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode) 
{
  int precy, m, inexact, sign;
  mpfr_t c;
  mp_exp_t e;
  MPFR_ZIV_DECL (loop);
  
  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", x, x, rnd_mode),
		  ("y[%#R]=%R inexact=%d", y, y, inexact));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x) || MPFR_IS_INF (x))
	{
	  MPFR_SET_NAN (y);
	  MPFR_RET_NAN;
	}
      else /* x is zero */
	{
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
	  MPFR_SET_ZERO (y);
	  MPFR_SET_SAME_SIGN (y, x);
	  MPFR_RET (0);
	}
    }

  /* Special case when |x| is very small (and x is a good approximation) */
  e = MPFR_GET_EXP(x);
  precy = MPFR_PREC(y);
  if (e < 0)
  {
    long v;

    /*
    Assume WMLOG that x > 0 (we will compensate for negative x later); then
    x - x^3 < Sin(x) < x.  Let e be the exponent of x, px be the precision
    of x and p the precision of y.  Let L and U be such that L < x <= U
    with L and U "adjacent" representable values in precision p.  (It
    follows that U = x rounded up to precision p.)

    Then a coarse estimate establishes that if Max(p,px) + 2*e - 2 <= 0
    then L < Sin(x) < U.  If not rounding to nearest then this is sufficient,
    otherwise we need to know where it stands with respect to (L+U)/2.  A
    simple approach is to implicitly use one extra bit so that one of L or U
    is implicitly replaced by (L+U)/2.  In this case, we require that
    Max(p+1,px) + 2*e - 2 <= 0.
    */
    v = (long)MPFR_PREC(y);
    if (rnd_mode == GMP_RNDN)
      v++;
    if (v < (long)MPFR_PREC(x))
      v = (long)MPFR_PREC(x);
    v -= 2;
    v += (long)e;	/* don't use v += 2*e due to overflow paranoia */
    if (v > 0)
      v += (long)e;

    if (v <= 0)		/* small case */
    {
      /*
      If rounding to nearest there is a potential problem if x lies exactly
      halfway between two representable values: the correct end result is
      obtained by rounding towards zero in this case.  Since this cannot
      currently be specified to the setting/rounding routines, an explicit
      check and change of rounding mode is done in this case.
      */
      if (rnd_mode == GMP_RNDN && MPFR_PREC(x) > precy)
      {
	mp_limb_t *xp, ok, mask;
	mp_prec_t sh;
	int k;

	MPFR_UNSIGNED_MINUS_MODULO(sh, precy);
	k = MPFR_LIMB_SIZE(x) - MPFR_LIMB_SIZE(y);
	xp = MPFR_MANT(x) + k;
	if (MPFR_LIKELY(sh != 0))
	{
	  mask = MPFR_LIMB_ONE << (sh - 1);
	}
	else
	{
	  mask = MPFR_LIMB_HIGHBIT;
	  xp--;
	  k--;
	}
	ok = (xp[0] & mask) ^ mask;
	if (ok == 0)
	{
	  ok = xp[0] & (mask-1);
	  if (MPFR_UNLIKELY(ok == 0))
	  {
	    while (--k >= 0 && ok == 0)
	      ok = *--xp;
	  }
	}

	if (MPFR_UNLIKELY(ok == 0))		/* problem case */
	  rnd_mode = GMP_RNDZ;
      }

      sign = MPFR_INT_SIGN(x);
      inexact = mpfr_set(y, x, rnd_mode);
      if (inexact == 0)
      {
	inexact = sign;
	if (MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG_SIGN(sign)))
	{
	  inexact = -inexact;
	  mpfr_nexttozero(y);
	}
      }
      MPFR_RET(inexact);
    }
  }

  /* Compute initial precision */
  m = precy + MPFR_INT_CEIL_LOG2 (precy) + 13;
  m += (e < 0) ? -2*e : e;  

  sign = mpfr_sin_sign (x);
  mpfr_init2 (c, m);

  MPFR_ZIV_INIT (loop, m);
  for (;;)
    {
      mpfr_cos (c, x, GMP_RNDZ);    /* can't be exact */
      mpfr_nexttoinf (c);           /* now c = cos(x) rounded away */
      mpfr_mul (c, c, c, GMP_RNDU); /* away */
      mpfr_ui_sub (c, 1, c, GMP_RNDZ);
      mpfr_sqrt (c, c, GMP_RNDZ);
      if (MPFR_IS_NEG_SIGN(sign))
	MPFR_CHANGE_SIGN(c);

      /* Warning c may be 0 ! */
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (c)))
	{
	  /* Huge cancellation: increase prec a lot! */
	  m = MAX (m, MPFR_PREC (x));
	  m = 2*m;
	}
      else
	{
	  /* the absolute error on c is at most 2^(3-m-EXP(c)) */
	  e = 2 * MPFR_GET_EXP (c) + m - 3;
	  if (mpfr_can_round (c, e, GMP_RNDZ, rnd_mode, precy))
	    break;

	  /* check for huge cancellation (Near 0) */
	  if (e < (mp_exp_t) MPFR_PREC (y))
	    m += MPFR_PREC (y) - e;
	  /* Check if near 1 */
	  if (MPFR_GET_EXP (c) == 1)
	    m += m;
	}
      
      /* Else generic increase */
      MPFR_ZIV_NEXT (loop, m);
      mpfr_set_prec (c, m);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set (y, c, rnd_mode);

  /* sin(x) is exact only for x = 0, which was treated apart above;
     nevertheless, we can have inexact = 0 here if the approximation c
     is exactly representable with PREC(y) bits. Since c is an approximation
     towards zero, in that case the inexact flag should have the opposite sign
     as y. */
  if (MPFR_UNLIKELY (inexact == 0))
    inexact = -MPFR_INT_SIGN (y);

  mpfr_clear (c);

  return inexact; /* inexact */
}
