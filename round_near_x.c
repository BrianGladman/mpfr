/* mpfr_round_near_x -- Round a floating point number nears another one.

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

#include "mpfr-impl.h"

/* Use MPFR_FAST_COMPUTE_IF_SMALL_INPUT instead (a simple wrapper) */

/* int mpfr_round_near_x (mpfr_ptr y, mpfr_srcptr x, mpfr_uexp_t err, int dir,
                          mp_rnd_t rnd)

   Assuming y = o(f(x)) = o(x + g(x)) with |g(x)| < 2^(EXP(x)-error)
   If x is small enough, y ~= x. This function checks and does this.

   It assumes that f(x) is not representable exactly as a FP number.
   x must not be a singular value (NAN, INF or ZERO).

   y is the destination (a mpfr_t), x the value to set (a mpfr_t),
   err the error term (a mpfr_uexp_t) such that |g(x)| < 2^(EXP(x)-err),
   dir (an int) is the direction of the error (if dir = 0,
   it rounds towards 0, if dir=1, it rounds away from 0),
   rnd the rounding mode.

   It returns 0 if it can't round.
   Otherwise it returns the ternary flag (It can't return an exact value).
*/

/* What "small enough" means?

   We work with the positive values.
   Assuming err > Prec (y)+1

   i = [ y = o(x)]   // i = inexact flag
   If i == 0
       Setting x in y is exact. We have:
       y = [XXXXXXXXX[...]]0[...] + error where [..] are optional zeros
      if dirError = ToInf,
        x < f(x) < x + 2^(EXP(x)-err)
        since x=y, and ulp (y)/2 > 2^(EXP(x)-err), we have:
        y < f(x) < y+ulp(y) and |y-f(x)| < ulp(y)/2
       if rnd = RNDN, nothing
       if rnd = RNDZ, nothing
       if rnd = RNDA, addoneulp
      elif dirError = ToZero
        x -2^(EXP(x)-err) < f(x) < x
        since x=y, and ulp (y)/2 > 2^(EXP(x)-err), we have:
        y-ulp(y) < f(x) < y and |y-f(x)| < ulp(y)/2
       if rnd = RNDN, nothing
       if rnd = RNDZ, nexttozero
       if rnd = RNDA, nothing
     NOTE: err > prec (y)+1 is needed only for RNDN.
   elif i > 0 and i = EVEN_ROUNDING
      So rnd = RNDN and we have y = x + ulp(y)/2
       if dirError = ToZero,
         we have x -2^(EXP(x)-err) < f(x) < x
         so y - ulp(y)/2 - 2^(EXP(x)-err) < f(x) < y-ulp(y)/2
         so y -ulp(y) < f(x) < y-ulp(y)/2
         => nexttozero(y)
       elif dirError = ToInf
         we have x < f(x) < x + 2^(EXP(x)-err)
         so y - ulp(y)/2 < f(x) < y+ulp(y)/2-ulp(y)/2
         so y - ulp(y)/2 < f(x) < y
         => do nothing
   elif i < 0 and i = -EVEN_ROUNDING
      So rnd = RNDN and we have y = x - ulp(y)/2
      if dirError = ToZero,
        y < f(x) < y + ulp(y)/2 => do nothing
      if dirError = ToInf
        y + ulp(y)/2 < f(x) < y + ulp(y) => AddOneUlp
   elif i > 0
     we can't have rnd = RNDZ, and prec(x) > prec(y), so ulp(x) < ulp(y)
     we have y - ulp (y) < x < y
     or more exactly y - ulp(y) + ulp(x)/2 <= x <= y - ulp(x)/2
     if rnd = RNDA,
      if dirError = ToInf,
       we have x < f(x) < x + 2^(EXP(x)-err)
       if err > prec (x),
         we have 2^(EXP(x)-err) < ulp(x), so 2^(EXP(x)-err) <= ulp(x)/2
         so f(x) <= y - ulp(x)/2+ulp(x)/2 <= y
         and y - ulp(y) < x < f(x)
         so we have y - ulp(y) < f(x) < y
         so do nothing.
       elif we can round, ie y - ulp(y) < x + 2^(EXP(x)-err) < y
         we have y - ulp(y) < x <  f(x) < x + 2^(EXP(x)-err) < y
         so do nothing
       otherwise
         Wrong. Example X=[0.11101]111111110000
                         +             1111111111111111111....
      elif dirError = ToZero
       we have x - 2^(EXP(x)-err) < f(x) < x
       so f(x) < x < y
       if err > prec (x)
         x-2^(EXP(x)-err) >= x-ulp(x)/2 >= y - ulp(y) + ulp(x)/2-ulp(x)/2
         so y - ulp(y) < f(x) < y
         so do nothing
       elif we can round, ie y - ulp(y) < x - 2^(EXP(x)-err) < y
         y - ulp(y) < x - 2^(EXP(x)-err) < f(x) < y
         so do nothing
       otherwise
        Wrong. Example: X=[1.111010]00000010
                         -             10000001000000000000100....
     elif rnd = RNDN,
      y - ulp(y)/2 < x < y and we can't have x = y-ulp(y)/2:
      so we have:
       y - ulp(y)/2 + ulp(x)/2 <= x <= y - ulp(x)/2
      if dirError = ToInf
        we have x < f(x) < x+2^(EXP(x)-err) and ulp(y) > 2^(EXP(x)-err)
        so y - ulp(y)/2 + ulp (x)/2 < f(x) < y + ulp (y)/2 - ulp (x)/2
        we can round but we can't compute inexact flag.
        if err > prec (x)
          y - ulp(y)/2 + ulp (x)/2 < f(x) < y + ulp(x)/2 - ulp(x)/2
          so y - ulp(y)/2 + ulp (x)/2 < f(x) < y
          we can round and compute inexact flag. do nothing
        elif we can round, ie y - ulp(y)/2 < x + 2^(EXP(x)-err) < y
          we have  y - ulp(y)/2 + ulp (x)/2 < f(x) < y
          so do nothing
        otherwise
          Wrong
      elif dirError = ToZero
        we have x -2^(EXP(x)-err) < f(x) < x and ulp(y)/2 > 2^(EXP(x)-err)
        so y-ulp(y)+ulp(x)/2 < f(x) < y - ulp(x)/2
        if err > prec (x)
           x- ulp(x)/2 < f(x) < x
           so y - ulp(y)/2+ulp(x)/2 - ulp(x)/2 < f(x) < x <= y - ulp(x)/2 < y
           do nothing
        elif we can round, ie y-ulp(y)/2 < x-2^(EXP(x)-err) < y
           we have y-ulp(y)/2 < x-2^(EXP(x)-err) < f(x) < x < y
           do nothing
        otherwise
          Wrong
   elif i < 0
     same thing?
 */

int
mpfr_round_near_x (mpfr_ptr y, mpfr_srcptr x, mpfr_uexp_t err, int dir,
                   mp_rnd_t rnd)
{
  int inexact, sign;
  unsigned int old_flags = __gmpfr_flags;

  MPFR_ASSERTD (!MPFR_IS_SINGULAR (x));
  MPFR_ASSERTD (dir == 0 || dir == 1);

  /* First check if we can round. The test is more restrictive than
     necessary. Note that if err is not representable in an mp_exp_t,
     then err > MPFR_PREC (x) and the conversion to mp_exp_t will not
     occur. */
  if (!(err > MPFR_PREC (y) + 1
        && (err > MPFR_PREC (x)
            || mpfr_round_p (MPFR_MANT (x), MPFR_LIMB_SIZE (x),
                             (mp_exp_t) err,
                             MPFR_PREC (y) + (rnd == GMP_RNDN)))))
    /* If we assume we can not round, return 0 */
    return 0;

  /* First round x in y */
  sign = MPFR_SIGN (x);
  MPFR_SET_EXP (y, MPFR_GET_EXP (x));
  MPFR_SET_SIGN (y, sign);
  MPFR_RNDRAW_GEN (inexact, y, MPFR_MANT (x), MPFR_PREC (x), rnd, sign,
                   if (dir == 0)
                     {
                       inexact = -sign;
                       goto trunc_doit;
                     }
                   else
                     goto addoneulp;
                   , if (MPFR_UNLIKELY (++MPFR_EXP (y) > __gmpfr_emax))
                       mpfr_overflow (y, rnd, sign)
                  );

  /* Fix it in some cases */
  MPFR_ASSERTD (!MPFR_IS_NAN (y) && !MPFR_IS_ZERO (y));
  /* If inexact == 0, setting y from x is exact but we haven't
     take into account yet the error term */
  if (inexact == 0)
    {
      if (dir == 0) /* The error term is negative for x positive */
        {
          inexact = sign;
          if (MPFR_IS_LIKE_RNDZ (rnd, MPFR_IS_NEG_SIGN (sign)))
            {
            nexttozero:
              /* The underflow flag should be set if the result is zero */
              __gmpfr_flags = old_flags;
              inexact = -sign;
              mpfr_nexttozero (y);
              if (MPFR_UNLIKELY (MPFR_IS_ZERO (y)))
                mpfr_set_underflow ();
            }
        }
      else /* The error term is positive for x positive */
        {
          inexact = -sign;
          /* Round Away */
          if (rnd != GMP_RNDN && rnd != GMP_RNDZ
              && MPFR_IS_RNDUTEST_OR_RNDDNOTTEST (rnd, MPFR_IS_POS_SIGN(sign)))
            {
            nexttoinf:
              /* The overflow flag should be set if the result is infinity */
              inexact = sign;
              mpfr_nexttoinf (y);
              if (MPFR_UNLIKELY (MPFR_IS_INF (y)))
                mpfr_set_overflow ();
            }
        }
    }

  MPFR_RET (inexact);
}
