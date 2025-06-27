/* mpfr_get_float16 -- convert a multiple precision floating-point
                       number to a _Float16 number

Copyright 2012-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

#include "mpfr-impl.h"

#ifdef MPFR_WANT_FLOAT16

#include <stdint.h>

typedef union { _Float16 x; uint16_t n; } b16u16;

_Float16
mpfr_get_float16 (mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{

  mpfr_t y;
  mpfr_exp_t e;
  int16_t m;
  b16u16 v;
  MPFR_SAVE_EXPO_DECL (expo);

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    return (_Float16) mpfr_get_d (x, rnd_mode);

  e = mpfr_get_exp (x); /* 2^(e-1) <= |x| < 2^e */

  if (e > 16) /* |x| >= 2^16 */
    {
      static const uint16_t s[2] = {0x7c00, 0xfc00};
      int neg = mpfr_signbit (x);
      v.n = s[neg];
      if (MPFR_IS_LIKE_RNDZ(rnd_mode,neg))
        v.n--;
      return v.x;
    }

  /* now x is a normal non-zero number, with |x| < 2^16 */
  MPFR_SAVE_EXPO_MARK (expo);

  mpfr_init2 (y, MPFR_PREC(x));

  /* we round x*2^(11-e) to an integer to get the significand of the result,
     except when x is in the subnormal range */
  if (e <= -14) /* subnormal range */
    {
      /* divide x by 2^-24 which is the smallest positive subnormal */
      mpfr_mul_2si (y, x, 24, MPFR_RNDN); /* exact */
      m = mpfr_get_si (y, rnd_mode);
      /* the result is m*2^-24 */
      MPFR_ASSERTD(-0x400 <= m && m <= 0x400);
      // the code below also works in the case where |m| = 0x400
      v.n = (mpfr_signbit (y)) ? 0x8000 + (-m) : m;
    }
  else
    {
      /* x is in the normal range */

      mpfr_mul_2si (y, x, 11 - e, MPFR_RNDN); /* exact */
      /* 2^10 <= |y| < 2^11 */
      m = mpfr_get_si (y, rnd_mode);
      /* 2^10 <= |m| <= 2^11 with 1 <= 14 + e <= 30 */
      v.n = ((14 + e) << 10) + ((m < 0) ? 0x7c00 - m : m - 0x400);
      /* Note: using + instead of | above allows the code to also work in case of
         overflow: when e=16 and m=0x800 for example, the exponent part is
         30 << 10 while the significand part is 0x400, which adds to 0x7c00,
         which is the encoding of +Inf. When e=16 and m=-0x800, the significant
         part is 0x7c00 + 0x800 = 0x8400, which added to 30 << 10 yields 0xfc00,
         which is the encoding of -Inf. */
  }

  mpfr_clear (y);
  MPFR_SAVE_EXPO_FREE (expo);
  return v.x;
}

#endif /* MPFR_WANT_FLOAT16 */
