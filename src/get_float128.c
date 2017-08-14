/* mpfr_get_float128 -- convert a multiple precision floating-point
                        number to a __float128 number

Copyright 2012-2017 Free Software Foundation, Inc.
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

/* FIXME: The current code depends on double's.
 * mpfr_get_ui should be used instead of mpfr_get_d.
 */

#ifdef MPFR_WANT_FLOAT128

/* generic code */
__float128
mpfr_get_float128 (mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    return (__float128) mpfr_get_d (x, rnd_mode);
  else /* now x is a normal non-zero number */
    {
      __float128 r; /* result */
      __float128 m;
      double s; /* part of result */
      mpfr_exp_t e;  /* exponent of x (before rounding) */
      mpfr_exp_t sh; /* exponent shift, so that x/2^sh is in the double range */
      mpfr_t y, z;
      const int emin = -16381;
      const int esub = emin - IEEE_FLOAT128_MANT_DIG;
      int prec;
      int sign;

      sign = MPFR_SIGN (x);
      e = MPFR_GET_EXP (x);

      if (MPFR_UNLIKELY (e <= esub))
        {
          if (MPFR_IS_LIKE_RNDZ (rnd_mode, sign < 0) ||
              (rnd_mode == MPFR_RNDN && (e < esub || mpfr_powerof2_raw (x))))
            return sign < 0 ? -0.0 : 0.0;
          r = 1.0;
          sh = esub;
        }
      else
        {
          /* First round x to the target __float128 precision, taking the
             reduced precision of the subnormals into account, so that all
             subsequent operations are exact (this avoids double rounding
             problems). */
          prec = e < emin ? e - esub : IEEE_FLOAT128_MANT_DIG;
          MPFR_ASSERTD (prec >= MPFR_PREC_MIN);
          mpfr_init2 (y, prec);
          mpfr_init2 (z, prec);

          mpfr_set (y, x, rnd_mode);
          sh = MPFR_GET_EXP (y);
          MPFR_SET_EXP (y, 0);
          MPFR_SET_POS (y);

          r = 0.0;
          do
            {
              s = mpfr_get_d (y, MPFR_RNDN); /* high part of y */
              r += (__float128) s;
              mpfr_set_d (z, s, MPFR_RNDN);  /* exact */
              mpfr_sub (y, y, z, MPFR_RNDN); /* exact */
            }
          while (MPFR_NOTZERO (y));

          mpfr_clear (z);
          mpfr_clear (y);
        }

      /* we now have to multiply r by 2^sh */
      MPFR_ASSERTD (r > 0);
      if (sh != 0)
        {
          /* An overflow may occur (example: 0.5*2^1024) */
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

#endif /* MPFR_WANT_FLOAT128 */
