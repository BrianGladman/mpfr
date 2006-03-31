/* mpfr_get_f -- convert a MPFR number to a GNU MPF number

Copyright 2005, 2006 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include "mpfr-impl.h"

/* return value is 0 iff no error occurred in the conversion
   (1 for NaN, +Inf, -Inf that have no equivalent in mpf)
*/
int
mpfr_get_f (mpf_ptr x, mpfr_srcptr y, mp_rnd_t rnd_mode)
{
  mp_size_t sx, sy;
  mp_prec_t precx, precy;
  int sh;

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
  precx = (mp_prec_t) sx * BITS_PER_MP_LIMB;
  sy = MPFR_LIMB_SIZE (y);

  /* since mpf numbers are represented in base 2^BITS_PER_MP_LIMB,
     we loose -EXP(y) % BITS_PER_MP_LIMB bits in the most significant limb */
  sh = MPFR_GET_EXP(y) % BITS_PER_MP_LIMB;
  sh = sh <= 0 ? - sh : BITS_PER_MP_LIMB - sh;
  MPFR_ASSERTD (sh >= 0);
  if (precy + sh <= precx) /* we can copy directly */
    {
      MPFR_ASSERTN (sx >= sy);
      if (sh != 0)
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
      mp_size_t sz;

      mpfr_init2 (z, precx - sh);
      sz = MPFR_LIMB_SIZE (z);
      mpfr_set (z, y, rnd_mode);
      /* warning, sh may change due to rounding, but then z is a power of two,
         thus we can safely ignore its last bit which is 0 */
      sh = MPFR_GET_EXP(z) % BITS_PER_MP_LIMB;
      sh = sh <= 0 ? - sh : BITS_PER_MP_LIMB - sh;
      MPFR_ASSERTD (sh >= 0);
      if (sh != 0)
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
