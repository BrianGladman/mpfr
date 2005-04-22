/* generic inverse of a function.

Copyright 2005 Free Software Foundation, Inc.

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

#include "mpfr-impl.h"

/* example of use:
   #define FUNCTION mpfr_sec
   #define INVERSE  mpfr_cos
   #define ACTION_NAN(y) do { MPFR_SET_NAN(y); MPFR_RET_NAN; } while (1)
   #define ACTION_INF(y) do { MPFR_SET_NAN(y); MPFR_RET_NAN; } while (1)
   #define ACTION_ZERO(y) return mpfr_set_ui (y, 1, GMP_RNDN)
   #include "gen_inverse.h"
*/

int
FUNCTION (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  if (MPFR_UNLIKELY(MPFR_IS_SINGULAR(x)))
    {
      if (MPFR_IS_NAN(x))
        ACTION_NAN(y);
      else if (MPFR_IS_INF(x))
        ACTION_INF(y);
      else /* x = 0 */
        ACTION_ZERO(y,x);
    }
  else /* x is neither NaN, Inf nor zero */
    {
      mp_prec_t precy; /* target precision */
      mp_prec_t m;     /* working precision */
      mpfr_t z;        /* temporary variable to store INVERSE(x) */
      int inexact;     /* inexact flag */
      MPFR_ZIV_DECL (loop);
      MPFR_SAVE_EXPO_DECL (expo);
      MPFR_SAVE_EXPO_MARK (expo);
      precy = MPFR_PREC(y);
      m = precy + __gmpfr_ceil_log2 ((double) precy) + 3;
      MPFR_ZIV_INIT (loop, m);
      mpfr_init2 (z, MPFR_PREC_MIN);
      for(;;)
        {
          mpfr_set_prec (z, m);
          INVERSE (z, x, GMP_RNDZ); /* error k_u < 1 ulp */
          mpfr_ui_div (z, 1, z, GMP_RNDN);
          /* the error is less than c_w + 2*c_u*k_u (see algorithms.tex),
             where c_w = 1/2, c_u = 1 since z was rounded towards zero,
             thus 1/2 + 2 < 4 */
          if (MPFR_LIKELY (MPFR_CAN_ROUND (z, m - 2, precy, rnd_mode)))
            break;

          MPFR_ZIV_NEXT (loop, m);
        }
      MPFR_ZIV_FREE (loop);
      inexact = mpfr_set (y, z, rnd_mode);
      mpfr_clear (z);
      MPFR_SAVE_EXPO_FREE (expo);
      MPFR_RET (mpfr_check_range (y, inexact, rnd_mode));
    }
}

