/* mpfr_agm -- arithmetic-geometric mean of two floating-point numbers

Copyright 1999, 2000, 2001, 2002, 2003, 2004 Free Software Foundation.

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

#include "mpfr-impl.h"

int
mpfr_agm (mpfr_ptr r, mpfr_srcptr op2, mpfr_srcptr op1, mp_rnd_t rnd_mode)
{
  int s, go_on, compare, inexact;
  mp_prec_t p, q;
#if 0
  double uo, vo;
#endif
  mp_limb_t *up, *vp, *tmpp, *tmpup, *tmpvp, *ap, *bp;
  mpfr_t u, v, tmp, tmpu, tmpv, a, b;
  TMP_DECL(marker);

  /* Deal with special values */
  if (MPFR_ARE_SINGULAR(op1, op2))
    {
      /* If a or b is NaN, the result is NaN */
      if (MPFR_IS_NAN(op1) || MPFR_IS_NAN(op2))
	{
	  MPFR_SET_NAN(r);
	  MPFR_RET_NAN;
	}
      /* now one of a or b is Inf or 0 */
      /* If a and b is +Inf, the result is +Inf.
	 Otherwise if a or b is -Inf, the result is NaN */
      else if (MPFR_IS_INF(op1) || MPFR_IS_INF(op2))
	{
          if (MPFR_IS_STRICTPOS(op1) && MPFR_IS_STRICTPOS(op2))
            {
              MPFR_SET_INF(r);
              MPFR_SET_SAME_SIGN(r, op1);
              MPFR_RET(0); /* exact */
            }
	  else
            {
              MPFR_SET_NAN(r);
              MPFR_RET_NAN;
            }
	}
      else /* a and b are neither NaN nor Inf, and one is zero */
	{  /* If a or b is 0, the result is +0 since a sqrt is positive */
          MPFR_ASSERTD(MPFR_IS_ZERO(op1) || MPFR_IS_ZERO(op2));
	  MPFR_SET_POS(r);
	  MPFR_SET_ZERO(r);
	  MPFR_RET(0); /* exact */
	}
    }
  MPFR_CLEAR_FLAGS(r);

  /* If a or b is negative (excluding -Infinity), the result is NaN */
  if (MPFR_UNLIKELY(MPFR_IS_NEG(op1)  || MPFR_IS_NEG(op2)))
    {
      MPFR_SET_NAN(r);
      MPFR_RET_NAN;
    }

  /* precision of the following calculus */
  q = MPFR_PREC(r);
  p = q + 15;

  /* Initialisations */
  go_on = 1;

  TMP_MARK(marker);
  s = (p - 1) / BITS_PER_MP_LIMB + 1;
  MPFR_TMP_INIT(ap, a, p, s);
  MPFR_TMP_INIT(bp, b, p, s);
  MPFR_TMP_INIT(up, u, p, s);
  MPFR_TMP_INIT(vp, v, p, s);
  MPFR_TMP_INIT(tmpup, tmpu, p, s);
  MPFR_TMP_INIT(tmpvp, tmpv, p, s);
  MPFR_TMP_INIT(tmpp, tmp, p, s);

  /* b and a are the 2 operands but we want b >= a */
  if ((compare = mpfr_cmp (op1, op2)) > 0)
    {
      mpfr_set (b,op1,GMP_RNDN);
      mpfr_set (a,op2,GMP_RNDN);
    }
  else if (compare == 0)
    {
      mpfr_set (r, op1, rnd_mode);
      TMP_FREE(marker);
      MPFR_RET(0); /* exact */
    }
  else
    {
      mpfr_set (b, op2, GMP_RNDN);
      mpfr_set (a, op1, GMP_RNDN);
    }

#if 0
  vo = mpfr_get_d1 (b);
  uo = mpfr_get_d1 (a);
#endif

  mpfr_set (u, a, GMP_RNDN);
  mpfr_set (v, b, GMP_RNDN);

  /* Main loop */

  while (go_on) 
    {
      int err, can_round;
      mp_prec_t eq;
      double erraux;
      
#if 0
      erraux = (vo == uo) ? 0.0
	: __gmpfr_ceil_exp2 (-2.0 * (double) p * uo / (vo - uo));
      err = 1 + (int) ((3.0 / 2.0 * (double) __gmpfr_ceil_log2 ((double) p)
			+ 1.0) * __gmpfr_ceil_exp2 (- (double) p)
		       + 3.0 * erraux);
#else
      /* since the argument of __gmpfr_ceil_exp2() is negative,
	 erraux is always bounded by 1.
	 We also have floor((3/2*ceil(log2(p))+1) * 2^(-p)) = 0,
	 since it the value f(p) inside floor() is 5/8 for p=2, 1/2 for p=3,
	 and f(2p) <= f(p)/2 + 1/2^p. So err=4.
      */
      erraux = 1.0;
      err = 4;
#endif

#if 0 /* this can not happen since p = q+15, and err=4 */
      if (p - err - 3 <= q)
	{
	  p = q + err + 4;
	  err = 1 + (int) ((3.0 / 2.0 * __gmpfr_ceil_log2 ((double) p) + 1.0)
			   * __gmpfr_ceil_exp2 (- (double) p) + 3.0 * erraux);
	}
#endif
      
      /* Calculus of un and vn */
      do
	{
	  mpfr_mul (tmp, u, v, GMP_RNDN);
	  mpfr_sqrt (tmpu, tmp, GMP_RNDN);
	  mpfr_add (tmp, u, v, GMP_RNDN);
	  mpfr_div_2ui (tmpv, tmp, 1, GMP_RNDN);
	  mpfr_set (u, tmpu, GMP_RNDN);
	  mpfr_set (v, tmpv, GMP_RNDN);
	}
      while (mpfr_cmp2 (u, v, &eq) != 0 && eq <= p - 2);
      
      /* Roundability of the result */
      can_round = mpfr_can_round (v, p - err - 3, GMP_RNDN, GMP_RNDZ,
				  q + (rnd_mode == GMP_RNDN));
      
      if (can_round)
	go_on = 0;
      
      else 
	{
	  go_on = 1;
	  p += 5;
	  s = (p - 1) / BITS_PER_MP_LIMB + 1;
	  MPFR_TMP_INIT(up, u, p, s);
	  MPFR_TMP_INIT(vp, v, p, s);
	  MPFR_TMP_INIT(tmpup, tmpu, p, s);
	  MPFR_TMP_INIT(tmpvp, tmpv, p, s);
	  MPFR_TMP_INIT(tmpp, tmp, p, s);
	  mpfr_set (u, a, GMP_RNDN);
	  mpfr_set (v, b, GMP_RNDN);
	}
    }
  /* End of while */

  /* Setting of the result */

  inexact = mpfr_set (r, v, rnd_mode);

  /* Let's clean */
  TMP_FREE(marker);

  return inexact; /* agm(u,v) can be exact for u, v rational only for u=v.
		     Proof (due to Nicolas Brisebarre): it suffices to consider
		     u=1 and v<1. Then 1/AGM(1,v) = 2F1(1/2,1/2,1;1-v^2),
		     and a theorem due to G.V. Chudnovsky states that for x a
		     non-zero algebraic number with |x|<1, then
		     2F1(1/2,1/2,1;x) and 2F1(-1/2,1/2,1;x) are algebraically
		     independent over Q. */
}
