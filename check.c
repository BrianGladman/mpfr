/* mpfr_check -- Check if a floating-point number has not been corrupted.

Copyright 2003, 2004 Free Software Foundation, Inc.

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

/*
 * Check if x is a valid mpfr_t initializes by mpfr_init
 * Returns 0 if isn't valid 
 */
int
mpfr_check (mpfr_srcptr x)
{
  mp_size_t s, i;
  mp_limb_t tmp;
  volatile mp_limb_t *xm;
  int rw;

  /* Check Sign */
  if (MPFR_SIGN(x) != MPFR_SIGN_POS && MPFR_SIGN(x) != MPFR_SIGN_NEG)
    return 0;
  /* Check Precision */
  if ( (MPFR_PREC(x) < MPFR_PREC_MIN) || (MPFR_PREC(x) > MPFR_PREC_MAX))
    return 0;
  /* Check Mantissa */
  xm = MPFR_MANT(x);
  if (!xm)
    return 0;
  /* Check size of mantissa */
  s = MPFR_GET_ALLOC_SIZE(x);
  if (s<=0 || s > MP_SIZE_T_MAX || 
      MPFR_PREC(x) > ((mp_prec_t)s*BITS_PER_MP_LIMB))
    return 0;
  /* Acces all the mp_limb of the mantissa: may do a seg fault */
  for(i = 0 ; i < s ; i++)
    tmp = xm[i];
  /* Check if it isn't singular*/
  if (MPFR_IS_PURE_FP(x))
    {
      /* Check first mp_limb of mantissa (Must start with a 1 bit) */
      if ( ((xm[MPFR_LIMB_SIZE(x)-1])>>(BITS_PER_MP_LIMB-1)) == 0)
	return 0;
      /* Check last mp_limb of mantissa */
      rw = (MPFR_PREC(x) % BITS_PER_MP_LIMB);
      if (rw != 0)
	{
	  tmp = MPFR_LIMB_MASK (BITS_PER_MP_LIMB - rw);
	  if ((xm[0] & tmp) != 0)
	    return 0;
	}
      /* Check exponent range */
      if ((MPFR_EXP (x) < __gmpfr_emin) || (MPFR_EXP (x) > __gmpfr_emax))
	return 0;
    }
  else 
    {
      /* Singular value is zero, inf or nan */
      MPFR_ASSERTD(MPFR_IS_ZERO(x) || MPFR_IS_NAN(x) || MPFR_IS_INF(x));
    }
  return 1;
}

