/* mpfr_ui_sub -- subtract a floating-point number from an integer

Copyright 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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
mpfr_ui_sub (mpfr_ptr y, unsigned long int u, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t uu;
  mp_limb_t up[1];
  unsigned long cnt;

  if (MPFR_UNLIKELY (u == 0))
    return mpfr_neg (y, x, rnd_mode);

  if (MPFR_UNLIKELY(MPFR_IS_SINGULAR(x)))
    {
      if (MPFR_IS_NAN(x))
        {
          MPFR_SET_NAN(y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF(x))
        {
          /*  u - Inf = -Inf and u - -Inf = +Inf  */
          MPFR_SET_INF(y);
          MPFR_SET_OPPOSITE_SIGN(y,x);
          MPFR_RET(0); /* +/-infinity is exact */
        }
      else /* x is zero */
        /* u - 0 = u */
        return mpfr_set_ui(y, u, rnd_mode);
    }
  else
    {
      MPFR_TMP_INIT1 (up, uu, BITS_PER_MP_LIMB);
      MPFR_ASSERTN(u == (mp_limb_t) u);
      count_leading_zeros (cnt, (mp_limb_t) u);
      *up = (mp_limb_t) u << cnt;
      MPFR_SET_EXP (uu, BITS_PER_MP_LIMB - cnt);
      return mpfr_sub (y, uu, x, rnd_mode);
    }
}
