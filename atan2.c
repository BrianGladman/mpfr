/* mpfr_atan2 -- arc-tan 2 of a floating-point number

Copyright 2005 Free Software Foundation.

This file is part of the MPFR Library, and was contributed by Mathieu Dutour.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* Returns arctan (y/x) */
int
mpfr_atan2 (mpfr_ptr dest, mpfr_srcptr y, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t tmp;
  int inexact;
  mp_prec_t prec;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("y[%#R]=%R x[%#R]=%R rnd=%d", y, y, x, x, rnd_mode),
		 ("atan[%#R]=%R inexact=%d", dest, dest, inexact));

  /* Special cases */
  if (MPFR_ARE_SINGULAR (x, y))
    {
      /* Don't want to handle all the cases. 
	 Just call div and atan is simplier, and still fast */
      mpfr_div (dest, y, x, rnd_mode);
      return mpfr_atan (dest, dest, rnd_mode);
    }
  MPFR_SAVE_EXPO_MARK (expo);

  /* Set up initial prec */
  prec = MPFR_PREC (dest) + 3 + MPFR_INT_CEIL_LOG2 (MPFR_PREC (dest));
  mpfr_init2 (tmp, prec);

  /* use atan2(y,x) = atan(y/x) */
  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      mpfr_div (tmp, y, x, GMP_RNDN);   /* Error <= ulp (tmp) */
      mpfr_atan (tmp, tmp, GMP_RNDN);   /* Error <= 2*ulp (tmp) since
					   abs(D(arctan)) <= 1 */
      /*FIXME: Error <= ulp(tmp) ? */
      if (MPFR_LIKELY (MPFR_CAN_ROUND (tmp, prec - 2, MPFR_PREC (dest),
				       rnd_mode)))
	break;
      MPFR_ZIV_NEXT (loop, prec);
      mpfr_set_prec (tmp, prec);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set (dest, tmp, rnd_mode);
  mpfr_clear (tmp);
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (dest, inexact, rnd_mode);
}
