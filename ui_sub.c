/* mpfr_ui_sub -- divide a machine integer by a floating-point number

Copyright (C) 2000 Free Software Foundation.

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
#include "mpfr-impl.h"

void
#if __STDC__
mpfr_ui_sub (mpfr_ptr y, unsigned long int u, mpfr_srcptr x, mp_rnd_t rnd_mode)
#else
mpfr_ui_sub (y, u, x, rnd_mode)
     mpfr_ptr y;
     unsigned long int u;
     mpfr_srcptr x;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t uu;
  mp_limb_t up[1];
  unsigned long cnt;

  if (MPFR_IS_NAN(x)) 
    {
      MPFR_SET_NAN(y);
      return;
    }

  MPFR_CLEAR_NAN(y);

  if (MPFR_IS_INF(x))
    {
      MPFR_SET_INF(y);
      if (MPFR_SIGN(x) == MPFR_SIGN(y)) MPFR_CHANGE_SIGN(y);
      return;
    }

  if (u) {
    MPFR_INIT1(up, uu, BITS_PER_MP_LIMB, 1);
    count_leading_zeros(cnt, (mp_limb_t) u);
    *up = (mp_limb_t) u << cnt;
    MPFR_EXP(uu) = BITS_PER_MP_LIMB-cnt;
  
    mpfr_sub (y, uu, x, rnd_mode);
  }
  else mpfr_neg (y, x, rnd_mode); /* if u=0, then set y to -x */
}
