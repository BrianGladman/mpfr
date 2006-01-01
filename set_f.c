/* mpfr_set_f -- set a MPFR number from a GNU MPF number

Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

int
mpfr_set_f (mpfr_ptr y, mpf_srcptr x, mp_rnd_t rnd_mode)
{
  mp_limb_t *my, *mx, *tmp;
  unsigned long cnt, sx, sy;
  int inexact, carry = 0;
  MPFR_TMP_DECL(marker);

  sx = ABS(SIZ(x)); /* number of limbs of the mantissa of x */

  if (sx == 0) /* x is zero */
    {
      MPFR_CLEAR_FLAGS (y);
      MPFR_SET_ZERO(y);
      MPFR_SET_POS(y);
      return 0; /* 0 is exact */
    }

  if (SIZ(x) * MPFR_FROM_SIGN_TO_INT(MPFR_SIGN(y)) < 0)
    MPFR_CHANGE_SIGN (y);

  MPFR_CLEAR_FLAGS (y);

  sy = 1 + (MPFR_PREC(y) - 1) / BITS_PER_MP_LIMB;
  my = MPFR_MANT(y);
  mx = PTR(x);

  count_leading_zeros(cnt, mx[sx - 1]);

  if (sy <= sx) /* we may have to round even when sy = sx */
    {
      unsigned long xprec = sx * BITS_PER_MP_LIMB;

      MPFR_TMP_MARK(marker);
      tmp = (mp_limb_t*) MPFR_TMP_ALLOC(sx * BYTES_PER_MP_LIMB);
      if (cnt)
        mpn_lshift (tmp, mx, sx, cnt);
      else
        /* FIXME: we may avoid the copy here, and directly call mpfr_round_raw
           on mx instead of tmp */
        MPN_COPY (tmp, mx, sx);
      carry = mpfr_round_raw (my, tmp, xprec, (SIZ(x) < 0), MPFR_PREC(y),
                              rnd_mode, &inexact);
      if (MPFR_UNLIKELY(carry)) /* result is a power of two */
        my[sy - 1] = MPFR_LIMB_HIGHBIT;
      MPFR_TMP_FREE(marker);
    }
  else
    {
      if (cnt)
        mpn_lshift (my + sy - sx, mx, sx, cnt);
      else
        MPN_COPY (my + sy - sx, mx, sx);
      MPN_ZERO(my, sy - sx);
      /* no rounding necessary, since y has a larger mantissa */
      inexact = 0;
    }

  /* warning: EXP(x) * BITS_PER_MP_LIMB may exceed the maximal exponent */
  if (EXP(x) > 1 + (__mpfr_emax - 1) / BITS_PER_MP_LIMB)
    {
      /* EXP(x) >= 2 + floor((__mpfr_emax-1)/BITS_PER_MP_LIMB)
	 EXP(x) >= 2 + (__mpfr_emax - BITS_PER_MP_LIMB) / BITS_PER_MP_LIMB
	        >= 1 + __mpfr_emax / BITS_PER_MP_LIMB
         EXP(x) * BITS_PER_MP_LIMB >= __mpfr_emax + BITS_PER_MP_LIMB
	 Since 0 <= cnt <= BITS_PER_MP_LIMB-1, and 0 <= carry <= 1,
	 we have then EXP(x) * BITS_PER_MP_LIMB - cnt + carry > __mpfr_emax */
      return mpfr_overflow (y, rnd_mode, MPFR_SIGN (y));
    }
  else
    MPFR_SET_EXP(y, EXP(x) * BITS_PER_MP_LIMB - (mp_exp_t) cnt + carry);

  return mpfr_check_range (y, inexact, rnd_mode);
}
