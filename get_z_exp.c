/* mpfr_get_z_exp -- get a multiple-precision integer and an exponent
                     from a floating-point number

Copyright 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

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


#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* puts the mantissa of f into z, and returns 'exp' such that f = z * 2^exp
 *
 * 0 doesn't have an exponent, therefore the returned exponent in this case
 * isn't really important. We choose to return __gmpfr_emin because
 *   1) it is in the exponent range [__gmpfr_emin,__gmpfr_emax],
 *   2) the smaller a number is (in absolute value), the smaller its
 *      exponent is. In other words, the f -> exp function is monotonous
 *      on nonnegative numbers.
 * Note that this is different from the C function frexp().
 */

mp_exp_t
mpfr_get_z_exp (mpz_ptr z, mpfr_srcptr f)
{
  mp_size_t fn;
  int sh;

  MPFR_ASSERTN(MPFR_IS_FP(f));

  if (MPFR_IS_ZERO(f))
    {
      mpz_set_ui (z, 0);
      return __gmpfr_emin;
    }

  fn = MPFR_LIMB_SIZE(f);

  /* check whether allocated space for z is enough */
  if (ALLOC(z) < fn)
    MPZ_REALLOC(z, fn);

  sh = (mp_prec_t) fn * BITS_PER_MP_LIMB - MPFR_PREC(f);
  if (sh)
    mpn_rshift (PTR(z), MPFR_MANT(f), fn, sh);
  else
    MPN_COPY (PTR(z), MPFR_MANT(f), fn);

  SIZ(z) = MPFR_SIGN(f) < 0 ? -fn : fn;

  /* Test if the result is representable. Later, we could choose
     to return MPFR_EXP_MIN if it isn't, or perhaps MPFR_EXP_MAX
     to signal an error. The mantissa would still be meaningful. */
  MPFR_ASSERTN((mp_exp_unsigned_t) MPFR_GET_EXP (f) - MPFR_EXP_MIN
               >= (mp_exp_unsigned_t) MPFR_PREC(f));

  return MPFR_GET_EXP (f) - MPFR_PREC (f);
}
