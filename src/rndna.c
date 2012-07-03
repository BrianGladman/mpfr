/* mpfr_round_nearest_away -- round to nearest away

Copyright 2012 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

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

/* Note: this doesn't work for 2^(emin-2). Currently, the best that can be
   done is to extend the exponent range internally, and do not support the
   case emin = MPFR_EMIN_MIN from the caller. */

/* put in rop the value of foo(op), rounded to nearest-away,
   and return the ternary value */
int
mpfr_round_nearest_away (mpfr_t rop, mpfr_srcptr op,
                         int foo(mpfr_t, mpfr_srcptr, mpfr_rnd_t))
{
  mpfr_t tmp;
  int inex;
  mpfr_prec_t n = mpfr_get_prec (rop);
  MPFR_SAVE_EXPO_DECL (expo);

  /* we check emin has not the smallest possible value, otherwise we cannot
     determine the correct round-to-nearest-away rounding for
     0.25*2^emin_min, which gets rounded to 0 with nearest-even,
     like 0.24*2^emin_min */
  MPFR_ASSERTN(mpfr_get_emin () > mpfr_get_emin_min ());

  MPFR_SAVE_EXPO_MARK (expo);

  mpfr_init2 (tmp, n + 1);

  /* first round to n+1 bits with rounding to nearest-even */
  inex = foo (tmp, op, MPFR_RNDN);

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (tmp)))
    mpfr_set (rop, tmp, MPFR_RNDN); /* inex unchanged */
  else
    {
      int lastbit, sh;

      MPFR_UNSIGNED_MINUS_MODULO(sh, n + 1);
      lastbit = (MPFR_MANT(tmp)[0] >> sh) & 1;

      if (lastbit == 0)
        mpfr_set (rop, tmp, MPFR_RNDN); /* exact, inex unchanged */
      else if (inex == 0)  /* midpoint: round away from zero */
        inex = mpfr_set (rop, tmp, MPFR_RNDA);
      else  /* lastbit == 1, inex != 0: double rounding */
        inex = mpfr_set (rop, tmp, (inex > 0) ? MPFR_RNDD : MPFR_RNDU);
    }

  mpfr_clear (tmp);
  MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, __gmpfr_flags);
  MPFR_SAVE_EXPO_FREE (expo);
  return inex;
}
