/* mpfr_get_f -- convert a MPFR number to a GNU MPF number

Copyright 2005 Free Software Foundation, Inc.

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

#include <stdio.h>
#include "mpfr-impl.h"

/* return value is 0 iff no error occurred in the conversion
   (1 for NaN, +Inf, -Inf that have no equivalent in mpf)
*/
int
mpfr_get_f (mpf_ptr x, mpfr_srcptr y, mp_rnd_t rnd_mode)
{
  unsigned long sx, sy, precx, precy, sh;
  mp_exp_t ey;

  if (MPFR_UNLIKELY(MPFR_IS_SINGULAR(y)))
    {
      if (MPFR_IS_ZERO(y))
        {
          mpf_set_ui (x, 0);
          return 0;
        }
      else /* NaN or Inf */
        return 1;
    }

  sx = PREC(x); /* number of limbs of the mantissa of x */

  precy = MPFR_PREC(y);
  precx = sx * BITS_PER_MP_LIMB;
  sy = 1 + (MPFR_PREC(y) - 1) / BITS_PER_MP_LIMB;

  /* since mpf numbers are represented in base 2^BITS_PER_MP_LIMB,
     we loose -EXP(y) % BITS_PER_MP_LIMB bits in the most significant limb */
  ey = MPFR_GET_EXP(y) % BITS_PER_MP_LIMB;
  if (ey <= 0)
    sh = (unsigned long) (-ey);
  else /* 0 < ey < BITS_PER_MP_LIMB */
    sh = BITS_PER_MP_LIMB - (unsigned long) ey;
  if (precy + sh <= precx) /* we can copy directly */
    {
      /* necessarily sy <= sx */
      if (sh)
        mpn_rshift (PTR(x) + sx - sy, MPFR_MANT(y), sy, sh);
      else
        MPN_COPY (PTR(x) + sx - sy, MPFR_MANT(y), sy);
      if (sx > sy)
        MPN_ZERO (PTR(x), sx - sy);
      EXP(x) = (MPFR_GET_EXP(y) + sh) / BITS_PER_MP_LIMB;
    }
  else /* we have to round to precx - sh bits */
    {
      mpfr_t z;
      unsigned long sz;

      mpfr_init2 (z, precx - sh);
      sz = 1 + (MPFR_PREC(z) - 1) / BITS_PER_MP_LIMB;
      mpfr_set (z, y, rnd_mode);
      /* warning, sh may change due to rounding, but then z is a power of two,
         thus we can safely ignore its last bit which is 0 */
      ey = MPFR_GET_EXP(z) % BITS_PER_MP_LIMB;
      sh = (ey <= 0) ? (unsigned long) (-ey)
        : BITS_PER_MP_LIMB - (unsigned long) ey;
      if (sh)
        mpn_rshift (PTR(x) + sx - sz, MPFR_MANT(z), sz, sh);
      else
        MPN_COPY (PTR(x) + sx - sz, MPFR_MANT(z), sz);
      if (sx > sz)
        MPN_ZERO (PTR(x), sx - sz);
      EXP(x) = (MPFR_GET_EXP(z) + sh) / BITS_PER_MP_LIMB;
      mpfr_clear (z);
    }

  /* set size and sign */
  SIZ(x) = (MPFR_FROM_SIGN_TO_INT(MPFR_SIGN(y)) < 0) ? -sx : sx;

  return 0;
}
