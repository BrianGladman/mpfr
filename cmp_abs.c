/* mpfr_cmpabs -- compare the absolute values of two FP numbers

Copyright 1999, 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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

#include "mpfr-impl.h"

/* Return a positive value if abs(b) > abs(c), 0 if abs(b) = abs(c), and
   a negative value if abs(b) < abs(c). Neither b nor c may be NaN. */

int
mpfr_cmpabs (mpfr_srcptr b, mpfr_srcptr c)
{
  mp_exp_t be, ce;
  mp_size_t bn, cn;
  mp_limb_t *bp, *cp;

  if (MPFR_ARE_SINGULAR (b, c))
    {
      if (MPFR_IS_NAN (b) || MPFR_IS_NAN (c))
	{
	  MPFR_SET_ERANGE ();
	  return 0;
	}
      else if (MPFR_IS_INF (b))
	return ! MPFR_IS_INF (c);
      else if (MPFR_IS_INF (c))
	return -1;
      else if (MPFR_IS_ZERO (c))
	return ! MPFR_IS_ZERO (b);
      else /* b == 0 */
	return -1;
    }

  be = MPFR_GET_EXP (b);
  ce = MPFR_GET_EXP (c);
  if (be > ce)
    return 1;
  if (be < ce)
    return -1;

  /* exponents are equal */

  bn = MPFR_LIMB_SIZE(b)-1;
  cn = MPFR_LIMB_SIZE(c)-1;

  bp = MPFR_MANT(b);
  cp = MPFR_MANT(c);

  for ( ; bn >= 0 && cn >= 0; bn--, cn--)
    {
      if (bp[bn] > cp[cn])
        return 1;
      if (bp[bn] < cp[cn])
        return -1;
    }

  for ( ; bn >= 0; bn--)
    if (bp[bn])
      return 1;

  for ( ; cn >= 0; cn--)
    if (cp[cn])
      return -1;

   return 0;
}
