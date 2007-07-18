/* mpfr_atan2 -- arc-tan 2 of a floating-point number

Copyright 2005, 2006, 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

int
mpfr_atan2 (mpfr_ptr dest, mpfr_srcptr y, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t tmp, pi;
  int inexact;
  mp_prec_t prec;
  mp_exp_t e;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("y[%#R]=%R x[%#R]=%R rnd=%d", y, y, x, x, rnd_mode),
                 ("atan[%#R]=%R inexact=%d", dest, dest, inexact));

  /* Special cases */
  if (MPFR_ARE_SINGULAR (x, y))
    {
      /* atan2(0, 0) does not raise the "invalid" floating-point
         exception, nor does atan2(y, 0) raise the "divide-by-zero"
         floating-point exception.
         -- atan2(�0, -0) returns �pi.313)
         -- atan2(�0, +0) returns �0.
         -- atan2(�0, x) returns �pi, for x < 0.
         -- atan2(�0, x) returns �0, for x > 0.
         -- atan2(y, �0) returns -pi/2 for y < 0.
         -- atan2(y, �0) returns pi/2 for y > 0.
         -- atan2(�oo, -oo) returns �3pi/4.
         -- atan2(�oo, +oo) returns �pi/4.
         -- atan2(�oo, x) returns �pi/2, for finite x.
         -- atan2(�y, -oo) returns �pi, for finite y > 0.
         -- atan2(�y, +oo) returns �0, for finite y > 0.
      */
      if (MPFR_IS_NAN (x) || MPFR_IS_NAN (y))
        {
          MPFR_SET_NAN (dest);
          MPFR_RET_NAN;
        }
      if (MPFR_IS_ZERO (y))
        {
          if (MPFR_IS_NEG (x)) /* +/- PI */
            {
            set_pi:
              if (MPFR_IS_NEG (y))
                {
                  inexact =  mpfr_const_pi (dest, MPFR_INVERT_RND (rnd_mode));
                  MPFR_CHANGE_SIGN (dest);
                  return -inexact;
                }
              else
                return mpfr_const_pi (dest, rnd_mode);
            }
          else /* +/- 0 */
            {
            set_zero:
              MPFR_SET_ZERO (dest);
              MPFR_SET_SAME_SIGN (dest, y);
              return 0;
            }
        }
      if (MPFR_IS_ZERO (x))
        {
        set_pi_2:
          if (MPFR_IS_NEG (y)) /* -PI/2 */
            {
              inexact = mpfr_const_pi (dest, MPFR_INVERT_RND(rnd_mode));
              MPFR_CHANGE_SIGN (dest);
              mpfr_div_2ui (dest, dest, 1, rnd_mode);
              return -inexact;
            }
          else /* PI/2 */
            {
              inexact = mpfr_const_pi (dest, rnd_mode);
              mpfr_div_2ui (dest, dest, 1, rnd_mode);
              return inexact;
            }
        }
      if (MPFR_IS_INF (y))
        {
          if (!MPFR_IS_INF (x)) /* +/- PI/2 */
            goto set_pi_2;
          else if (MPFR_IS_POS (x)) /* +/- PI/4 */
            {
              if (MPFR_IS_NEG (y))
                {
                  rnd_mode = MPFR_INVERT_RND (rnd_mode);
                  inexact = mpfr_const_pi (dest, rnd_mode);
                  MPFR_CHANGE_SIGN (dest);
                  mpfr_div_2ui (dest, dest, 2, rnd_mode);
                  return -inexact;
                }
              else
                {
                  inexact = mpfr_const_pi (dest, rnd_mode);
                  mpfr_div_2ui (dest, dest, 2, rnd_mode);
                  return inexact;
                }
            }
          else /* +/- 3*PI/4: Ugly since we have to round properly */
            {
              mpfr_t tmp;
              MPFR_ZIV_DECL (loop);
              mp_prec_t prec = MPFR_PREC (dest) + BITS_PER_MP_LIMB;

              mpfr_init2 (tmp, prec);
              MPFR_ZIV_INIT (loop, prec);
              for (;;)
                {
                  mpfr_const_pi (tmp, GMP_RNDN);
                  mpfr_mul_ui (tmp, tmp, 3, GMP_RNDN); /* Error <= 2  */
                  mpfr_div_2ui (tmp, tmp, 2, GMP_RNDN);
                  if (mpfr_round_p (MPFR_MANT (tmp), MPFR_LIMB_SIZE (tmp),
                                    MPFR_PREC (tmp)-2,
                                    MPFR_PREC (dest) + (rnd_mode == GMP_RNDN)))
                    break;
                  MPFR_ZIV_NEXT (loop, prec);
                  mpfr_set_prec (tmp, prec);
                }
              MPFR_ZIV_FREE (loop);
              if (MPFR_IS_NEG (y))
                MPFR_CHANGE_SIGN (tmp);
              inexact = mpfr_set (dest, tmp, rnd_mode);
              mpfr_clear (tmp);
              return inexact;
            }
        }
      MPFR_ASSERTD (MPFR_IS_INF (x));
      if (MPFR_IS_NEG (x))
        goto set_pi;
      else
        goto set_zero;
    }

  /* When x=1, atan2(y,x) = atan(y). FIXME: more generally, if x is a power
     of two, we could call directly atan(y/x) since y/x is exact. */
  if (mpfr_cmp_ui (x, 1) == 0)
    return mpfr_atan (dest, y, rnd_mode);

  MPFR_SAVE_EXPO_MARK (expo);

  /* Set up initial prec */
  prec = MPFR_PREC (dest) + 3 + MPFR_INT_CEIL_LOG2 (MPFR_PREC (dest));
  mpfr_init2 (tmp, prec);

  MPFR_ZIV_INIT (loop, prec);
  if (MPFR_IS_POS (x))
    /* use atan2(y,x) = atan(y/x) */
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
  else /* x < 0 */
    /*  Use sign(y)*(PI - atan (|y/x|)) */
    {
      mpfr_init2 (pi, prec);
      for (;;)
        {
          mpfr_div (tmp, y, x, GMP_RNDN);   /* Error <= ulp (tmp) */
          MPFR_SET_POS (tmp);               /* no error */
          mpfr_atan (tmp, tmp, GMP_RNDN);   /* Error <= 2*ulp (tmp) since
                                               abs(D(arctan)) <= 1 */
          mpfr_const_pi (pi, GMP_RNDN);     /* Error <= ulp(pi) /2 */
          e = MPFR_GET_EXP (tmp);
          mpfr_sub (tmp, pi, tmp, GMP_RNDN);          /* see above */
          if (MPFR_IS_NEG (y))
            MPFR_CHANGE_SIGN (tmp);
          /* Error(tmp) <= (1/2+2^(EXP(pi)-EXP(tmp)-1)+2^(e-EXP(tmp)+1))*ulp
                        <= 2^(MAX (MAX (EXP(PI)-EXP(tmp)-1, e-EXP(tmp)+1),
                                        -1)+2)*ulp(tmp) */
          e = MAX (MAX (MPFR_GET_EXP (pi)-MPFR_GET_EXP (tmp) - 1,
                        e - MPFR_GET_EXP (tmp) + 1), -1) + 2;
          if (MPFR_LIKELY (MPFR_CAN_ROUND (tmp, prec - e, MPFR_PREC (dest),
                                           rnd_mode)))
            break;
          MPFR_ZIV_NEXT (loop, prec);
          mpfr_set_prec (tmp, prec);
          mpfr_set_prec (pi, prec);
        }
      mpfr_clear (pi);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set (dest, tmp, rnd_mode);
  mpfr_clear (tmp);
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (dest, inexact, rnd_mode);
}
