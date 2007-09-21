/* mpfr_sinh_cosh -- hyperbolic sine and cosine

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

 /* The computations are done by
    cosh(x) = 1/2 [e^(x)+e^(-x)]
    sinh(x) = 1/2 [e^(x)-e^(-x)]
    Adapted from mpfr_sinh.c     */

int
mpfr_sinh_cosh (mpfr_ptr sh, mpfr_ptr ch, mpfr_srcptr xt, mp_rnd_t rnd_mode)
{
  mpfr_t x;
  int inexact, inexact_sh, inexact_ch;

  MPFR_ASSERTN (sh != ch);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", xt, xt, rnd_mode),
                 ("sh[%#R]=%R ch[%#R]=%R inexact=%d", sh, sh, ch, ch, inexact));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (xt)))
    {
      if (MPFR_IS_NAN (xt))
        {
          MPFR_SET_NAN (ch);
          MPFR_SET_NAN (sh);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (xt))
        {
          MPFR_SET_INF (sh);
          MPFR_SET_SAME_SIGN (sh, xt);
          MPFR_SET_INF (ch);
          MPFR_SET_POS (ch);
          MPFR_RET (0);
        }
      else /* xt is zero */
        {
          MPFR_ASSERTD (MPFR_IS_ZERO (xt));
          MPFR_SET_ZERO (sh);                   /* sinh(0) = 0 */
          MPFR_SET_SAME_SIGN (sh, xt);
          return mpfr_set_ui (ch, 1, rnd_mode); /* cosh(0) = 1 */
        }
    }

  MPFR_TMP_INIT_ABS (x, xt);

  {
    mpfr_t s, c, ti;
    mp_exp_t d;
    mp_prec_t N;    /* Precision of the intermediary variables */
    long int err;    /* Precision of error */
    MPFR_ZIV_DECL (loop);
    MPFR_SAVE_EXPO_DECL (expo);
    MPFR_GROUP_DECL (group);

    MPFR_SAVE_EXPO_MARK (expo);

    /* compute the precision of intermediary variable */
    N = MPFR_PREC (ch);
    N = MAX (N, MPFR_PREC (sh));
    N = MAX (N, MPFR_PREC (x));
    /* the optimal number of bits : see algorithms.ps */
    N = N + MPFR_INT_CEIL_LOG2 (N) + 4;

    /* initialise of intermediary variables */
    MPFR_GROUP_INIT_3 (group, N, s, c, ti);

    /* First computation of sinh_cosh */
    MPFR_ZIV_INIT (loop, N);
    for (;;) {
      /* compute sinh_cosh */
      mpfr_clear_flags ();
      mpfr_exp (s, x, GMP_RNDD);        /* exp(x) */
      /* exp(x) can overflow! */
      /* BUG/TODO/FIXME: exp can overflow but sinh or cosh may be 
	 representable! */
      if (MPFR_UNLIKELY (mpfr_overflow_p ())) {
        inexact_ch = mpfr_overflow (ch, rnd_mode, MPFR_SIGN_POS);
        inexact_sh = mpfr_overflow (sh, rnd_mode, MPFR_SIGN (xt));
        MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, MPFR_FLAGS_OVERFLOW);
        break;
      }
      d = MPFR_GET_EXP (s);
      mpfr_ui_div (ti, 1, s, GMP_RNDU);  /* 1/exp(x) */
      mpfr_add (c, s, ti, GMP_RNDU);     /* exp(x) + 1/exp(x) */
      mpfr_sub (s, s, ti, GMP_RNDN);     /* exp(x) - 1/exp(x) */
      mpfr_div_2ui (c, c, 1, GMP_RNDN);  /* 1/2(exp(x) + 1/exp(x)) */
      mpfr_div_2ui (s, s, 1, GMP_RNDN);  /* 1/2(exp(x) - 1/exp(x)) */

      /* it may be that s is zero (in fact, it can only occur when exp(x)=1,
         and thus ti=1 too) */
      if (MPFR_IS_ZERO (s))
        err = N; /* double the precision */
      else
        {
          /* calculation of the error */
          d = d - MPFR_GET_EXP (s) + 2;
          /* error estimate: err = N-(__gmpfr_ceil_log2(1+pow(2,d)));*/
          err = N - (MAX (d, 0) + 1);
          if (MPFR_LIKELY (MPFR_CAN_ROUND (s, err, MPFR_PREC (sh), rnd_mode) \
		       && (MPFR_CAN_ROUND (c, err, MPFR_PREC (ch), rnd_mode))))
            {
              inexact_sh = mpfr_set4 (sh, s, rnd_mode, MPFR_SIGN (xt));
	      inexact_ch = mpfr_set (ch, c, rnd_mode);
              break;
            }
        }
      /* actualisation of the precision */
      N += err;
      MPFR_ZIV_NEXT (loop, N);
      MPFR_GROUP_REPREC_3 (group, N, s, c, ti);
    }
    MPFR_ZIV_FREE (loop);
    MPFR_GROUP_CLEAR (group);
    MPFR_SAVE_EXPO_FREE (expo);
  }
    /* now, let's raise the flags if needed */
    inexact = mpfr_check_range (sh, inexact_sh, rnd_mode);
    inexact = mpfr_check_range (ch, inexact_ch, rnd_mode) || inexact;

    return inexact;
}
