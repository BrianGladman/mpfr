/* mpfr_sqrt_ui -- square root of a machine integer

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

#define MON_INIT(xp, x, p, s) xp = (mp_ptr) TMP_ALLOC(s*BYTES_PER_MP_LIMB); \
   MPFR_PREC(x) = p; MPFR_MANT(x) = xp; x -> _mp_size = s;

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "longlong.h"

int
#if __STDC__
mpfr_sqrt_ui (mpfr_ptr r, unsigned long u, mp_rnd_t rnd_mode)
#else
mpfr_sqrt_ui (r, u, rnd_mode)
     mpfr_ptr r;
     unsigned long u;
     mp_rnd_t rnd_mode;
#endif
{
  int error = 0;
  mpfr_t uu;
  mp_limb_t *up;
  unsigned long cnt;
  TMP_DECL(marker);

  if (u) { /* if u=0, do nothing */
    TMP_MARK(marker);
    MON_INIT(up, uu, BITS_PER_MP_LIMB, 1);
    count_leading_zeros(cnt, (mp_limb_t) u);
    *up = (mp_limb_t) u << cnt;
    MPFR_EXP(uu) = BITS_PER_MP_LIMB-cnt;

    error = mpfr_sqrt(r, uu, rnd_mode);

    TMP_FREE(marker);
  }
  else {
    MPFR_CLEAR_FLAGS(r);
    MPFR_SET_ZERO(r);
  }
  return error;
}
