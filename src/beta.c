/* mpfr_beta -- beta function

Copyright 2017 Free Software Foundation, Inc.
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

#define MPFR_NEED_LONGLONG_H /* for MPFR_INT_CEIL_LOG2 */
#include "mpfr-impl.h"

/* use formula (6.2.2) from Abramowitz & Stegun:
   beta(z,w) = gamma(z)*gamma(w)/gamma(z+w) */
int
mpfr_beta (mpfr_ptr r, mpfr_srcptr z, mpfr_srcptr w, mpfr_rnd_t rnd_mode)
{
  mpfr_exp_t emin, emax;
  mpfr_prec_t pmin, prec;
  mpfr_t z_plus_w, tmp, tmp2;
  int inex;
  MPFR_GROUP_DECL (group);
  MPFR_ZIV_DECL (loop);

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (z) || MPFR_IS_SINGULAR (w)))
    {
      /* if z or w is NaN, return NaN */
      if (MPFR_IS_NAN (z) || MPFR_IS_NAN (w))
        {
          MPFR_SET_NAN (r);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (z) || MPFR_IS_INF (w))
        {
          /* by symmetry we can assume z >= w:
             if z = +Inf and w > 0, then r = +0 (including w = +Inf);
             if z = +Inf and w = 0, then r = NaN
               [beta(z,1/log(z)) tends to +Inf whereas
                beta(z,1/log(loz(z))) tends to +0]
             if z = +Inf and w < 0:
                if w is integer or -Inf: r = NaN (pole)
                if -2k-1 < w < -2k:   r = -Inf
                if -2k-2 < w < -2k-1: r = +Inf
             if z = -Inf (then w = -Inf too): r = NaN */
          if (mpfr_cmp (z, w) < 0)
            return mpfr_beta (r, w, z, rnd_mode);
          /* now z >= w */
          printf ("not yet implemented");
        }
      else /* z or w is 0 */
        {
        }
      return 0;
    }

  /* compute the smallest precision such that z + w is exact */
  emax = (MPFR_EXP(z) >= MPFR_EXP(w)) ? MPFR_EXP(z) : MPFR_EXP(w);
  emin = (MPFR_EXP(z) - MPFR_PREC(z) <= MPFR_EXP(w) - MPFR_PREC(w))
    ? MPFR_EXP(z) - MPFR_PREC(z) : MPFR_EXP(w) - MPFR_PREC(w);
  pmin = emax - emin;
  mpfr_init2 (z_plus_w, pmin);
  inex = mpfr_add (z_plus_w, z, w, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  prec = MPFR_PREC(r);
  prec += MPFR_INT_CEIL_LOG2 (prec);
  MPFR_GROUP_INIT_2 (group, prec, tmp, tmp2);
  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      MPFR_GROUP_REPREC_2 (group, prec, tmp, tmp2);
      inex = mpfr_gamma (tmp, z, MPFR_RNDN);
      /* tmp = gamma(z) * (1 + theta) with |theta| <= 2^-prec */
      inex |= mpfr_gamma (tmp2, w, MPFR_RNDN);
      /* tmp2 = gamma(w) * (1 + theta2) with |theta2| <= 2^-prec */
      inex |= mpfr_mul (tmp, tmp, tmp2, MPFR_RNDN);
      /* tmp = gamma(z)*gamma(w) * (1 + theta3)^3 with |theta3| <= 2^-prec */
      inex |= mpfr_gamma (tmp2, z_plus_w, MPFR_RNDN);
      /* tmp2 = gamma(z+w) * (1 + theta4) with |theta4| <= 2^-prec */
      inex |= mpfr_div (tmp, tmp, tmp2, MPFR_RNDN);
      /* tmp = gamma(z)*gamma(w)/gamma(z+w) * (1 + theta5)^5
         with |theta5| <= 2^-prec. For prec >= 3, we have
         |(1 + theta5)^5 - 1| <= 7 * 2^(-prec), thus the error is bounded
         by 7 ulps */
      
      /* if inex=0, then tmp is exactly beta(z,w) */
      if (inex == 0 ||
          MPFR_LIKELY (MPFR_CAN_ROUND (tmp, prec - 3, MPFR_PREC(r), rnd_mode)))
        break;
      MPFR_ZIV_NEXT (loop, prec);
    }
  MPFR_ZIV_FREE (loop);
  inex = mpfr_set (r, tmp, rnd_mode);
  MPFR_GROUP_CLEAR (group);
  mpfr_clear (z_plus_w);
  return inex;
}
