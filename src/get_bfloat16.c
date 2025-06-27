/* mpfr_get_bfloat16 -- convert a multiple precision floating-point
                       number to a bfloat16 number

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

#ifdef MPFR_WANT_BFLOAT16

#include <stdint.h>

typedef union { __bf16 x; uint16_t n; } b16u16;

__bf16
mpfr_get_bfloat16 (mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{

  mpfr_t y;
  mpfr_exp_t e;
  int16_t m;
  b16u16 v;
  MPFR_SAVE_EXPO_DECL (expo);

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    return (__bf16) mpfr_get_d (x, rnd_mode);

  e = mpfr_get_exp (x); /* 2^(e-1) <= |x| < 2^e */

  if (e > 128) /* |x| >= 2^128 */
    {
      static const uint16_t s[2] = {0x7f80, 0xff80}; /* +Inf, -Inf */
      int neg = mpfr_signbit (x);
      v.n = s[neg];
      if (MPFR_IS_LIKE_RNDZ(rnd_mode,neg))
        v.n--;
      return v.x;
    }

  /* now x is a normal non-zero number, with |x| < 2^128 */
  MPFR_SAVE_EXPO_MARK (expo);

  mpfr_init2 (y, MPFR_PREC(x));

  /* we round x*2^(8-e) to an integer to get the significand of the result,
     except when x is in the subnormal range */
  if (e <= -126) /* subnormal range */
    {
      /* divide x by 2^-133 which is the smallest positive subnormal */
      mpfr_mul_2si (y, x, 133, MPFR_RNDN); /* exact */
      m = mpfr_get_si (y, rnd_mode);
      /* the result is m*2^-133 */
      MPFR_ASSERTD(-0x80 <= m && m <= 0x80);
      /* the code below also works in the case where |m| = 0x80 */
      v.n = (mpfr_signbit (y)) ? 0x8000 + (-m) : m;
    }
  else
    {
      /* x is in the normal range */

      mpfr_mul_2si (y, x, 8 - e, MPFR_RNDN); /* exact */
      /* 2^7 <= |y| < 2^8 */
      m = mpfr_get_si (y, rnd_mode);
      /* 2^7 <= |m| <= 2^8 with 1 <= 126 + e <= 126 */
      v.n = ((126 + e) << 7) + ((m < 0) ? 0x7f80 - m : m - 0x80);
      /* Note: using + instead of | above allows the code to also work in case
         of overflow: when e=128 and m=0x100 for example, the exponent part is
         254 << 7 while the significand part is 0x80, which adds to 0x7f80,
         which is the encoding of +Inf. When e=128 and m=-0x100, the
         significand part is 0x7c80 + 0x100 = 0x7d80, which added to
         254 << 7 yields 0xfc80, which is the encoding of -Inf. */
  }

  mpfr_clear (y);
  MPFR_SAVE_EXPO_FREE (expo);
  return v.x;
}

#endif /* MPFR_WANT_BFLOAT16 */
