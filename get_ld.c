/* mpfr_get_ld -- convert a multiple precision floating-point number
                  to a machine long double

Copyright 2002, 2003, 2004 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <float.h>


#include "mpfr-impl.h"

long double
mpfr_get_ld (mpfr_srcptr x, mp_rnd_t rnd_mode)
{

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    return (long double) mpfr_get_d (x, rnd_mode);
  else /* now x is a normal non-zero number */
    {
      long double r; /* result */
      long double m;
      double s; /* part of result */
      mp_exp_t sh; /* exponent shift, so that x/2^sh is in the double range */
      mpfr_t y, z;
      int sign;

      /* first round x to the target long double precision, so that
         all subsequent operations are exact (this avoids double rounding
         problems) */
      mpfr_init2 (y, MPFR_LDBL_MANT_DIG);
      mpfr_init2 (z, IEEE_DBL_MANT_DIG);
 
      mpfr_set (y, x, rnd_mode);
      sh = MPFR_GET_EXP (y);
      sign = MPFR_SIGN (y);
      MPFR_SET_EXP (y, 0);
      MPFR_SET_POS (y);

      r = 0.0;
      do {
	s = mpfr_get_d (y, GMP_RNDN); /* high part of y */
	r += (long double) s;
	mpfr_set_d (z, s, GMP_RNDN);  /* exact */
	mpfr_sub (y, y, z, GMP_RNDN); /* exact */
      } while (!MPFR_IS_ZERO (y));

      mpfr_clear (z);
      mpfr_clear (y);

      /* we now have to multiply back by 2^sh */
      MPFR_ASSERTD (r > 0);
      if (sh != 0)
        {
          /* An overflow may occurs (example: 0.5*2^1024) */
          while (r < 1.0)
            {
              r += r;
              sh--;
            }

          if (sh > 0)
            m = 2.0;
          else
            {
              m = 0.5;
              sh = -sh;
            }

          for (;;)
            {
              if (sh % 2)
                r = r * m;
              sh >>= 1;
              if (sh == 0)
                break;
              m = m * m;
            }
        }
      if (sign < 0)
	r = -r;
      return r;
    }
}
