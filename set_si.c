/* mpfr_set_si, mpfr_set_ui -- set a MPFR number from a machine integer

Copyright (C) 1999 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

void
#if __STDC__
mpfr_set_si(mpfr_ptr x, long int i, mp_rnd_t rnd_mode)
#else
mpfr_set_si(x, i, rnd_mode)
     mpfr_ptr x;
     long int i;
     mp_rnd_t rnd_mode;
#endif
{
  unsigned long xn, cnt; mp_limb_t ai, *xp;

  if (i==0) { MPFR_SET_ZERO(x); return; }
  xn = (MPFR_PREC(x)-1)/BITS_PER_MP_LIMB;
  ai = ABS(i); 

  count_leading_zeros(cnt, ai); 

  xp = MPFR_MANT(x);
  xp[xn] = ai << cnt;
  /* don't forget to put zero in lower limbs */
  MPN_ZERO(xp, xn);

  MPFR_EXP(x) = BITS_PER_MP_LIMB - cnt;

  /* round if MPFR_PREC(x) smaller than length of i */
  if (MPFR_PREC(x) < BITS_PER_MP_LIMB-cnt) {
    cnt = mpfr_round_raw(xp+xn, xp+xn, BITS_PER_MP_LIMB-cnt, (ai<0), MPFR_PREC(x), 
		   rnd_mode);
    if (cnt) { /* special case 1.000...000 */
      MPFR_EXP(x)++;
      xp[xn] = ((mp_limb_t) 1) << (BITS_PER_MP_LIMB-1);
    }
  }

  /* warning: don't change the precision of x! */
  if (i*MPFR_SIGN(x) < 0) MPFR_CHANGE_SIGN(x);

  return; 
}

void
#if __STDC__
mpfr_set_ui(mpfr_ptr x, unsigned long int i, mp_rnd_t rnd_mode)
#else
mpfr_set_ui(x, i, rnd_mode)
     mpfr_ptr x;
     long int i;
     mp_rnd_t rnd_mode;
#endif  
{
  unsigned int xn, cnt; mp_limb_t *xp;

  if (i==0) { MPFR_SET_ZERO(x); return; }
  xn = (MPFR_PREC(x)-1)/BITS_PER_MP_LIMB;
  count_leading_zeros(cnt, (mp_limb_t) i); 

  xp = MPFR_MANT(x);
  xp[xn] = ((mp_limb_t) i) << cnt; 
  /* don't forget to put zero in lower limbs */
  MPN_ZERO(xp, xn);

  MPFR_EXP(x) = BITS_PER_MP_LIMB - cnt;

  /* round if MPFR_PREC(x) smaller than length of i */
  if (MPFR_PREC(x) < BITS_PER_MP_LIMB-cnt) {
    cnt = mpfr_round_raw(xp+xn, xp+xn, BITS_PER_MP_LIMB-cnt, 0, MPFR_PREC(x), 
			 rnd_mode);
    if (cnt) { /* special case 1.000...000 */
      MPFR_EXP(x)++;
      xp[xn] = ((mp_limb_t) 1) << (BITS_PER_MP_LIMB-1);
    }
  }

  /* warning: don't change the precision of x! */
  if (MPFR_SIGN(x) < 0) MPFR_CHANGE_SIGN(x);

  return; 
}

