/* mpfr_cbrt -- cube root function.

Copyright 2002-2020 Free Software Foundation, Inc.
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
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

 /* The computation of y = x^(1/3) is done as follows:

    Let x = sign * m * 2^(3*e) where m is an integer

    with 2^(3n-3) <= m < 2^(3n) where n = PREC(y)

    and m = s^3 + r where 0 <= r and m < (s+1)^3

    FIXME: In the code below, the root can be applied to a value of m
    larger than 2^(3n), if PREC(x) > 3*PREC(y) up to some small constant.

    we want that s has n bits i.e. s >= 2^(n-1), or m >= 2^(3n-3)
    i.e. m must have at least 3n-2 bits

    then x^(1/3) = s * 2^e if r=0
         x^(1/3) = (s+1) * 2^e if round up
         x^(1/3) = (s-1) * 2^e if round down
         x^(1/3) = s * 2^e if nearest and r < 3/2*s^2+3/4*s+1/8
                   (s+1) * 2^e otherwise
 */

int
mpfr_cbrt (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  mpz_t m;
  mpfr_exp_t e, sh;
  mpfr_prec_t n, size_m, tmp;
  int inexact, inexact2, negative, r;
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_LOG_FUNC (
    ("x[%Pu]=%.*Rg rnd=%d", mpfr_get_prec (x), mpfr_log_prec, x, rnd_mode),
    ("y[%Pu]=%.*Rg inexact=%d", mpfr_get_prec (y), mpfr_log_prec, y,
     inexact));

  /* special values */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (x))
        {
          MPFR_SET_INF (y);
          MPFR_SET_SAME_SIGN (y, x);
          MPFR_RET (0);
        }
      /* case 0: cbrt(+/- 0) = +/- 0 */
      else /* x is necessarily 0 */
        {
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
          MPFR_SET_ZERO (y);
          MPFR_SET_SAME_SIGN (y, x);
          MPFR_RET (0);
        }
    }

  /* General case */
  MPFR_SAVE_EXPO_MARK (expo);
  mpz_init (m);

  e = mpfr_get_z_2exp (m, x);                /* x = m * 2^e */
  if ((negative = MPFR_IS_NEG(x)))
    mpz_neg (m, m);
  r = e % 3;
  if (r < 0)
    r += 3;
  MPFR_ASSERTD (r >= 0 && r < 3 && (e - r) % 3 == 0);

  /* x = (m*2^r) * 2^(e-r) = (m*2^r) * 2^(3*q) */

  MPFR_LOG_MSG (("e=%" MPFR_EXP_FSPEC "d r=%d\n", (mpfr_eexp_t) e, r));

  MPFR_MPZ_SIZEINBASE2 (size_m, m);
  n = MPFR_PREC (y) + (rnd_mode == MPFR_RNDN);

  /* FIXME: The division by 3 leaves a remainder in [-2,2], which is
     incompatible with the [-2,0] implied by the constraint below. */
  /* we want 3*n-2 <= size_m + 3*sh + r <= 3*n */
  sh = (3 * (mpfr_exp_t) n - (mpfr_exp_t) size_m - r) / 3;
  /* i.e. the following value should be in [-2,0] */
  MPFR_LOG_MSG (("size_m + 3*sh + r - 3*n = %" MPFR_EXP_FSPEC "d\n",
                 3 * sh - (3 * (mpfr_exp_t) n - (mpfr_exp_t) size_m - r)));

  sh = 3 * sh + r;
  if (sh > 0)
    {
      mpz_mul_2exp (m, m, sh);
      e -= sh;
    }
  else if (r > 0)
    {
      mpz_mul_2exp (m, m, r);
      e -= r;
    }

  MPFR_ASSERTD (e % 3 == 0);
  e /= 3;

  /* invariant: x = m*2^(3*e) */

  /* we reuse the variable m to store the cube root, since it is not needed
     any more: we just need to know if the root is exact */
  inexact = mpz_root (m, m, 3) == 0;

  MPFR_MPZ_SIZEINBASE2 (tmp, m);
  sh = tmp - n;
  if (sh > 0) /* we have to flush to 0 the last sh bits from m */
    {
      inexact = inexact || (mpz_scan1 (m, 0) < sh);
      mpz_fdiv_q_2exp (m, m, sh);
      e += sh;
    }

  if (inexact)
    {
      if (negative)
        rnd_mode = MPFR_INVERT_RND (rnd_mode);
      if (rnd_mode == MPFR_RNDU || rnd_mode == MPFR_RNDA
          || (rnd_mode == MPFR_RNDN && mpz_tstbit (m, 0)))
        {
          inexact = 1;
          mpz_add_ui (m, m, 1);
        }
      else
        inexact = -1;
    }

  /* either inexact is not zero, and the conversion is exact, i.e. inexact
     is not changed; or inexact=0, and inexact is set only when
     rnd_mode=MPFR_RNDN and bit (n+1) from m is 1 */
  inexact2 = mpfr_set_z (y, m, MPFR_RNDN);
  MPFR_ASSERTD (inexact == 0 || inexact2 == 0);
  inexact += inexact2;
  MPFR_SET_EXP (y, MPFR_GET_EXP (y) + e);

  if (negative)
    {
      MPFR_CHANGE_SIGN (y);
      inexact = -inexact;
    }

  mpz_clear (m);
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}
