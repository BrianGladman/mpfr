/* mpfr_div_ui -- divide a floating-point number by a machine integer

Copyright 1999, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

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
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"
#include "mpfr-impl.h"

/* returns 0 if result exact, non-zero otherwise */
int
mpfr_div_ui (mpfr_ptr y, mpfr_srcptr x, unsigned long int u, mp_rnd_t rnd_mode)
{
  long int xn, yn, dif, sh, i;
  mp_limb_t *xp, *yp, *tmp, c, d;
  mp_exp_t exp;
  int inexact, middle = 1;
  TMP_DECL(marker);

  if (MPFR_UNLIKELY( MPFR_IS_SINGULAR(x) ))
    {
      if (MPFR_IS_NAN(x))
	{
	  MPFR_SET_NAN(y);
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_INF(x)) 
	{ 
	  MPFR_SET_INF(y);
	  MPFR_SET_SAME_SIGN(y, x);
	  MPFR_RET(0);
	}
      else if (MPFR_IS_ZERO(x))
	{
	  if (u == 0)/* 0/0 is NaN */
	    {
	      MPFR_SET_NAN(y);
	      MPFR_RET_NAN;
	    }
	  else
	    {
	      MPFR_SET_ZERO(y);
	      MPFR_RET(0);
	    }
	}
      else
	MPFR_ASSERTN(1);
    }

  if (MPFR_UNLIKELY(u == 0))
    {
      /* x/0 is Inf */
      MPFR_SET_INF(y);
      MPFR_SET_SAME_SIGN(y, x);
      MPFR_RET(0);
    }

  MPFR_CLEAR_FLAGS(y);

  MPFR_SET_SAME_SIGN(y, x);

  TMP_MARK(marker);
  xn = MPFR_LIMB_SIZE(x);
  /*(MPFR_PREC(x) - 1)/BITS_PER_MP_LIMB + 1;*/
  yn = MPFR_LIMB_SIZE(y);
  /*(MPFR_PREC(y) - 1)/BITS_PER_MP_LIMB + 1;*/

  xp = MPFR_MANT(x);
  yp = MPFR_MANT(y);
  exp = MPFR_GET_EXP (x);

  dif = yn + 1 - xn;

  /* we need to store yn+1 = xn + dif limbs of the quotient */
  /* don't use tmp=yp since the mpn_lshift call below requires yp >= tmp+1 */
  tmp = (mp_limb_t*) TMP_ALLOC((yn + 1) * BYTES_PER_MP_LIMB);

  c = (mp_limb_t) u;
  MPFR_ASSERTN(u == c);
  if (dif >= 0)
    c = mpn_divrem_1 (tmp, dif, xp, xn, c); /* used all the dividend */
  else /* dif < 0 i.e. xn > yn, don't use the (-dif) low limbs from x */
    c = mpn_divrem_1 (tmp, 0, xp - dif, yn + 1, c);

  inexact = (c != 0);
  if (rnd_mode == GMP_RNDN)
    {
      if (2 * c < (mp_limb_t) u)
	middle = -1;
      else if (2 * c > (mp_limb_t) u)
	middle = 1;
      else
	middle = 0; /* exactly in the middle */
    }
  for (i=0; ((inexact == 0) || (middle == 0)) && (i < -dif); i++)
    if (xp[i])
      inexact = middle = 1; /* larger than middle */

  if (tmp[yn] == 0) /* high limb is zero */
    {
      tmp--;
      sh = 0;
      exp -= BITS_PER_MP_LIMB;
    }

  /* now we have yn limbs starting from tmp[1], with tmp[yn]<>0 */

  /* shift left to normalize */
  count_leading_zeros(sh, tmp[yn]);
  if (sh)
    {
      mpn_lshift (yp, tmp + 1, yn, sh);
      yp[0] += tmp[0] >> (BITS_PER_MP_LIMB - sh);
      middle = middle || ((tmp[0] << sh) != 0);
      inexact = inexact || ((tmp[0] << sh) != 0);
      exp -= sh;
    }
  else
    MPN_COPY(yp, tmp + 1, yn);

  sh = yn * BITS_PER_MP_LIMB - MPFR_PREC(y);
  /* it remains sh bits in less significant limb of y */

  d = *yp & ((MP_LIMB_T_ONE << sh) - MP_LIMB_T_ONE);
  *yp ^= d; /* set to zero lowest sh bits */

  MPFR_SET_EXP (y, exp);

  TMP_FREE(marker);
  if ((d == 0) && (inexact == 0))
    return 0; /* result is exact */

  switch (rnd_mode)
    {
    case GMP_RNDZ:
      MPFR_RET(-MPFR_INT_SIGN(x)); /* result is inexact */

    case GMP_RNDU:
      if (MPFR_IS_POS(y))
	mpfr_add_one_ulp (y, rnd_mode);
      MPFR_RET(1); /* result is inexact */

    case GMP_RNDD:
      if (MPFR_IS_NEG(y))
	mpfr_add_one_ulp (y, rnd_mode);
      MPFR_RET(-1); /* result is inexact */

    case GMP_RNDN:
      if (sh && d < (MP_LIMB_T_ONE << (sh - 1)))
	MPFR_RET(-MPFR_INT_SIGN(x));
      else if (sh && d > (MP_LIMB_T_ONE << (sh - 1)))
	{
	  mpfr_add_one_ulp (y, rnd_mode);
	  MPFR_RET(MPFR_INT_SIGN(x));
	}
      else /* sh = 0 or d = 1 << (sh-1) */
	{
	  /* we are in the middle if:
	     (a) sh > 0 and inexact == 0
	     (b) sh=0 and middle=1
	  */
	  if ((sh && inexact) || (!sh && (middle > 0)) ||
	      (*yp & (MP_LIMB_T_ONE << sh)))
	    {
	      mpfr_add_one_ulp (y, rnd_mode);
	      MPFR_RET(MPFR_INT_SIGN(x));
	    }
	  else
	    MPFR_RET(-MPFR_INT_SIGN(x));
	}
    }
  MPFR_ASSERTN(1); /* should never go here */
  return 0; /* To avoid warning*/
}
