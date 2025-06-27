/* mpfr_set_float16 -- convert a machine _Float116 number to
                       a multiple precision floating-point number

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

int
mpfr_set_float16 (mpfr_ptr r, _Float16 d, mpfr_rnd_t rnd_mode)
{
  b16u16 v;
  int e, sbit;
  int16_t m;

  v.x = d;
  e = (v.n >> 10) & 0x1f;
  sbit = v.n >> 15;
  m = v.n & 0x3ff;

  /*
    NaN is encoded by e=31 and m!=0
    +Inf is encoded by e=31 and m=0
    the largest number is 0x1.ffcp+15 (e=30, m=0x3ff)
    1.0 is encoded by e=15 and m=0
    the smallest positive normal number is 0x1p-14 (e=1, m=0)
    the largest subnormal number is 0x1.ff8p-15 (e=0, m=0x3ff)
    the smallest positive subnormal number is 0x1p-24 (e=0, m=1)
   */

  /* Check for NaN or INF */
  if (MPFR_UNLIKELY (e == 0x1f))
    {
      if (m != 0) /* NaN */
        {
          MPFR_SET_NAN(r);
          MPFR_RET_NAN;
        }
      /* INF case */
      MPFR_SET_INF (r);
      if (sbit) /* sign bit is set */
        MPFR_SET_NEG (r);
      else
        MPFR_SET_POS (r);
      return 0;
    }
  else if (MPFR_UNLIKELY (e == 0)) /* subnormal case */
    {
      if (m == 0) /* case +/-0 */
        return mpfr_set_d (r, (double) d, rnd_mode);
      e ++;
    }
  else
    m += 0x400; /* add implicit bit */

  if (sbit)
    m = -m;

  /* d = m * 2^(e-25) */
  return mpfr_set_si_2exp (r, m, e - 25, rnd_mode);
}

#endif /* MPFR_WANT_FLOAT16 */
