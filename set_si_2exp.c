/* mpfr_set_si_2exp -- set a MPFR number from a machine signed integer with
   a shift

Copyright 2004, 2006, 2007 Free Software Foundation, Inc.

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

#include <limits.h>
#include "mpfr-impl.h"

int
mpfr_set_si_2exp (mpfr_ptr x, long i, mp_exp_t e, mp_rnd_t rnd_mode)
{
  mpfr_t ii;
  int res;
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_SAVE_EXPO_MARK (expo);
  mpfr_init2 (ii, sizeof (long) * CHAR_BIT);
  res = mpfr_set_si (ii, i, rnd_mode);  /* exact, no exceptions */
  MPFR_ASSERTN (res == 0);
  MPFR_ASSERTN (e >= LONG_MIN && e <= LONG_MAX);
  /* FIXME: this may no longer be the case in the future. */
  res = mpfr_mul_2si (x, ii, e, rnd_mode);
  mpfr_clear (ii);
  MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, __gmpfr_flags);
  MPFR_SAVE_EXPO_FREE (expo);
  res = mpfr_check_range(x, res, rnd_mode);
  return res;
}
