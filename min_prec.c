/* mpfr_min_prec -- minimal size in bits to hold the mantissa

Copyright 2009 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "mpfr-impl.h"

int
mpfr_min_prec (mpfr_srcptr x)
{
  mp_limb_t *mx;
  mp_prec_t px;
  mp_size_t n;
  int res, i;


  if (MPFR_IS_NAN (x) || MPFR_IS_INF (x))
    return 0;

  if (MPFR_IS_ZERO (x))
    return 0;

  mx = MPFR_MANT (x);
  px = MPFR_PREC (x);

  res = 0;
  /* Count full limbs set to zero */
  for (n = (px - 1) / BITS_PER_MP_LIMB; mx[n] == 0; n--)
    {
      res += BITS_PER_MP_LIMB;
    }

  i = 0;
  /* mx[n] is now the first limb which is not null. Count number
   * of null bits in mx[n], from the right */
  while ((i < BITS_PER_MP_LIMB) && (mx[n] & (1 << i) == 0))
    i++;

  res += i;
  /* If we have trailing zero bits because the precision
   * is not a multiple of BITS_PER_MP_LIMB, we must not count
   * those. */
  i = px % BITS_PER_MP_LIMB;
  if (px != 0)
    res -= BITS_PER_MP_LIMB - i;

  return px - res;
}
