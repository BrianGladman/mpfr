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
  int s, compare, inexact;
  mp_prec_t p, q;
  mp_limb_t *up, *vp, *tmpup, *tmpvp;
  mpfr_t u, v, tmpu, tmpv;
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
  if (MPFR_UNLIKELY(MPFR_IS_NEG(op1) || MPFR_IS_NEG(op2)))
    {
      MPFR_SET_NAN(r);
      MPFR_RET_NAN;
    }

  /* precision of the following calculus */
  q = MPFR_PREC(r);
  p = q + 15;
  s = (p - 1) / BITS_PER_MP_LIMB + 1;

  /* b (op2) and a (op1) are the 2 operands but we want b >= a */
  compare = mpfr_cmp (op1, op2);
  if (MPFR_UNLIKELY( compare == 0 ))
    {
      mpfr_set (r, op1, rnd_mode);
      MPFR_RET (0); /* exact */
    }
  else if (compare > 0)
    {
      mpfr_srcptr t = op1;
      op1 = op2;
      op2 = t;
    }
  /* Now b(=op2) >= a (=op1) */

  TMP_MARK(marker);

  /* Main loop */
  for (;;)
    {
      mp_prec_t eq;

      /* Init temporary vars */
      MPFR_TMP_INIT (up, u, p, s);
      MPFR_TMP_INIT (vp, v, p, s);
      MPFR_TMP_INIT (tmpup, tmpu, p, s);
      MPFR_TMP_INIT (tmpvp, tmpv, p, s);

      /* Set initial values */
      mpfr_set (u, op1, GMP_RNDN);
      mpfr_set (v, op2, GMP_RNDN);

      /* Calculus of un and vn */
      do
	{
	  mpfr_mul (tmpu, u, v, GMP_RNDN);
	  mpfr_sqrt (tmpu, tmpu, GMP_RNDN);
	  mpfr_add (tmpv, u, v, GMP_RNDN);
	  mpfr_div_2ui (tmpv, tmpv, 1, GMP_RNDN);
	  mpfr_set (u, tmpu, GMP_RNDN);
	  mpfr_set (v, tmpv, GMP_RNDN);
	}
      while (mpfr_cmp2 (u, v, &eq) != 0 && eq <= p - 2);
      
      /* Roundability of the result */
      if (mpfr_can_round (v, p - 4 - 3, GMP_RNDN, GMP_RNDZ,
			  q + (rnd_mode == GMP_RNDN)))
	break; /* Stop the loop */
  
      /* Next iteration */
      p += 5;
      s = (p - 1) / BITS_PER_MP_LIMB + 1;
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
