/* mpfr_total_order -- total order of two floating-point numbers

Copyright 2018 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "mpfr-impl.h"

/* follows IEEE 754-2008, Section 5.10 */
int
mpfr_total_order (mpfr_srcptr x, mpfr_srcptr y)
{
  /* FIXME: test the sign bit first? */

  if (MPFR_IS_NAN(x))
    {
      if (!MPFR_IS_NAN(y)) /* -NaN < y */
        return MPFR_IS_POS(x) ? 0 : 1; /* +NaN > y, -NaN < y */
      /* both x and y are NaN */
      if (MPFR_SIGN(x) != MPFR_SIGN(y))
        return MPFR_IS_POS(x) ?	0 : 1; /* +NaN > -NaN, -NaN < +NaN */
      return 1;
    }

  /* now x is not NaN */
  if (MPFR_IS_NAN(y))
    return MPFR_IS_POS(y) ? 1 : 0; /* x < +NaN, x > -NaN */

  /* now neither x nor y are NaN */
  if (MPFR_IS_ZERO(x) && MPFR_IS_ZERO(y))
    {
      if (MPFR_SIGN(x) != MPFR_SIGN(y))
        return MPFR_IS_POS(x) ? 0 : 1; /* +0 > -0, -0 < +0 */
      else
        return 1;
    }

  return mpfr_cmp (x, y) <= 0 ? 1 : 0;
}
