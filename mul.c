/* mpfr_mul -- multiply two floating-point numbers

Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005
  Free Software Foundation, Inc.

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

int
mpfr_mul (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mp_rnd_t rnd_mode)
{
  int sign, inexact;
  mp_exp_t  ax, ax2;
  mp_limb_t *tmp;
  mp_limb_t b1;
  mp_prec_t bq, cq;
  mp_size_t bn, cn, tn, k;
  TMP_DECL (marker);

  /* deal with special cases */
  if (MPFR_ARE_SINGULAR (b, c))
    {
      if (MPFR_IS_NAN (b) || MPFR_IS_NAN (c))
	{
	  MPFR_SET_NAN (a);
	  MPFR_RET_NAN;
	}
      sign = MPFR_MULT_SIGN (MPFR_SIGN (b), MPFR_SIGN (c));
      if (MPFR_IS_INF (b))
	{
	  if (!MPFR_IS_ZERO (c))
	    {
	      MPFR_SET_SIGN (a, sign);
	      MPFR_SET_INF (a);
	      MPFR_RET (0);
	    }
	  else
	    {
	      MPFR_SET_NAN (a);
	      MPFR_RET_NAN;
	    }
	}
      else if (MPFR_IS_INF (c))
	{
	  if (!MPFR_IS_ZERO (b))
	    {
	      MPFR_SET_SIGN (a, sign);
	      MPFR_SET_INF (a);
	      MPFR_RET(0);
	    }
	  else
	    {
	      MPFR_SET_NAN (a);
	      MPFR_RET_NAN;
	    }
	}
      else
	{
	  MPFR_ASSERTD (MPFR_IS_ZERO(b) || MPFR_IS_ZERO(c));
	  MPFR_SET_SIGN (a, sign);
	  MPFR_SET_ZERO (a);
	  MPFR_RET (0);
	}
    }
  MPFR_CLEAR_FLAGS (a);
  sign = MPFR_MULT_SIGN (MPFR_SIGN (b), MPFR_SIGN (c));

  ax = MPFR_GET_EXP (b) + MPFR_GET_EXP (c);
  /* Note: the exponent of the exact result will be e = bx + cx + ec with
     ec in {-1,0,1} and the following assumes that e is representable. */

  /* FIXME: Useful since we do an exponent check after ?
   * It is useful iff the precision is big, there is an overflow
   * and we are doing further mults...*/
#ifdef HUGE
  if (MPFR_UNLIKELY (ax > __gmpfr_emax + 1))
    return mpfr_overflow (a, rnd_mode, sign);
  if (MPFR_UNLIKELY (ax < __gmpfr_emin - 2))
  return mpfr_underflow (a, rnd_mode == GMP_RNDN ? GMP_RNDZ : rnd_mode,
			     sign);
#endif

  bq = MPFR_PREC (b);
  cq = MPFR_PREC (c);

  MPFR_ASSERTD (bq+cq > bq); /* PREC_MAX is /2 so no integer overflow */

  bn = (bq+BITS_PER_MP_LIMB-1)/BITS_PER_MP_LIMB; /* number of limbs of b */
  cn = (cq+BITS_PER_MP_LIMB-1)/BITS_PER_MP_LIMB; /* number of limbs of c */
  k = bn + cn; /* effective nb of limbs used by b*c (= tn or tn+1) below */
  tn = (bq + cq + BITS_PER_MP_LIMB - 1) / BITS_PER_MP_LIMB;
  MPFR_ASSERTD (tn <= k); /* tn <= k, thus no int overflow */

  /* Check for no size_t overflow*/
  MPFR_ASSERTD ((size_t) k <= ((size_t) ~0) / BYTES_PER_MP_LIMB);
  TMP_MARK (marker);
  tmp = (mp_limb_t *) TMP_ALLOC ((size_t) k * BYTES_PER_MP_LIMB);

  /* multiplies two mantissa in temporary allocated space */
#if 0
  b1 = MPFR_LIKELY (bn >= cn)
    ? mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn)
    : mpn_mul (tmp, MPFR_MANT (c), cn, MPFR_MANT (b), bn);

  /* now tmp[0]..tmp[k-1] contains the product of both mantissa,
     with tmp[k-1]>=2^(BITS_PER_MP_LIMB-2) */
  b1 >>= BITS_PER_MP_LIMB - 1; /* msb from the product */

  /* if the mantissas of b and c are uniformly distributed in ]1/2, 1],
     then their product is in ]1/4, 1/2] with probability 2*ln(2)-1 ~ 0.386
     and in [1/2, 1] with probability 2-2*ln(2) ~ 0.614 */
  tmp += k - tn;
  if (MPFR_UNLIKELY (b1 == 0))
    mpn_lshift (tmp, tmp, tn, 1); /* tn <= k, so no stack corruption */

#else

  if (MPFR_UNLIKELY (bn < cn)) 
    {
      mpfr_srcptr tmp = b;
      mp_size_t tn  = bn;
      b = c;
      bn = cn;
      c = tmp;
      cn = tn;
    }
  if (MPFR_UNLIKELY (bn > MPFR_MUL_THRESHOLD))
    {
      mp_size_t n;
      mp_prec_t p;

      /* Compute estimated precision of mulhigh */
      n = MPFR_LIMB_SIZE (a) + 1;
      n = MIN (n, cn);
      MPFR_ASSERTD (n >= 1 && 2*n <= k && n <= cn && n <= bn);
      p = n*BITS_PER_MP_LIMB - MPFR_INT_CEIL_LOG2 (n + (n < cn) + (n < bn) );
      if (MPFR_UNLIKELY (MPFR_PREC (a) > p - 4))
	/* MulHigh can't produce a roundable result. */
	goto full_multiply;
      /* Compute an approximation of the product of b and c */
      mpfr_mulhigh_n (tmp+k-2*n, MPFR_MANT (b) + bn - n, 
		      MPFR_MANT (c) + cn - n, n);

      /* now tmp[0]..tmp[k-1] contains the product of both mantissa,
	 with tmp[k-1]>=2^(BITS_PER_MP_LIMB-2) */
      b1 = tmp[k-1] >> (BITS_PER_MP_LIMB - 1); /* msb from the product */
      
      /* if the mantissas of b and c are uniformly distributed in ]1/2, 1],
	 then their product is in ]1/4, 1/2] with probability 2*ln(2)-1 ~ 0.386
	 and in [1/2, 1] with probability 2-2*ln(2) ~ 0.614 */
      tmp += k - tn;
      if (MPFR_UNLIKELY (b1 == 0))
	mpn_lshift (tmp, tmp, tn, 1);

      if (MPFR_UNLIKELY (!mpfr_can_round_raw (tmp, tn, sign, p + b1 - 1,
	    GMP_RNDN, GMP_RNDZ, MPFR_PREC(a)+(rnd_mode==GMP_RNDN))))
	goto full_multiply;
      }
  else
    {
    full_multiply:
      b1 = mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn);      

      /* now tmp[0]..tmp[k-1] contains the product of both mantissa,
	 with tmp[k-1]>=2^(BITS_PER_MP_LIMB-2) */
      b1 >>= BITS_PER_MP_LIMB - 1; /* msb from the product */
      
      /* if the mantissas of b and c are uniformly distributed in ]1/2, 1],
	 then their product is in ]1/4, 1/2] with probability 2*ln(2)-1 ~ 0.386
	 and in [1/2, 1] with probability 2-2*ln(2) ~ 0.614 */
      tmp += k - tn;
      if (MPFR_UNLIKELY (b1 == 0))
	mpn_lshift (tmp, tmp, tn, 1); /* tn <= k, so no stack corruption */
    }
#endif

  ax2 = ax + (mp_exp_t) (b1 - 1);
  MPFR_RNDRAW (inexact, a, tmp, bq+cq, rnd_mode, sign, ax2++);
  TMP_FREE (marker);
  MPFR_EXP  (a) = ax2; /* Can't use MPFR_SET_EXP: Expo may be out of range */
  MPFR_SET_SIGN (a, sign);
  if (MPFR_UNLIKELY (ax2 > __gmpfr_emax))
    return mpfr_overflow (a, rnd_mode, sign);
  if (MPFR_UNLIKELY (ax2 < __gmpfr_emin))
    {
      /* In the rounding to the nearest mode, if the exponent of the exact
	 result (i.e. before rounding, i.e. without taking cc into account)
	 is < __gmpfr_emin - 1 or the exact result is a power of 2 (i.e. if
	 both arguments are powers of 2), then round to zero. */
      if (rnd_mode == GMP_RNDN
	  && (ax + (mp_exp_t) b1 < __gmpfr_emin
	      || (mpfr_powerof2_raw (b) && mpfr_powerof2_raw (c))))
	rnd_mode = GMP_RNDZ;
      return mpfr_underflow (a, rnd_mode, sign);
    }
  return inexact;
}

