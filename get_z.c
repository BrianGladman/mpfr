/* mpfr_get_z_exp -- get a multiple-precision integer from
                     a floating-point number

Copyright 2004, 2006, 2007, 2008 Free Software Foundation, Inc.
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

#include "mpfr-impl.h"

void
mpfr_get_z (mpz_ptr z, mpfr_srcptr f, mp_rnd_t rnd)
{
  mpfr_t r;
  mp_exp_t exp = MPFR_EXP (f);

  /* if exp <= 0, then |f|<1, thus |o(f)|<=1 */
  MPFR_ASSERTD (MPFR_PREC_MIN == 2);
  mpfr_init2 (r, (exp <= 0) ? MPFR_PREC_MIN : exp + 1);
  mpfr_rint (r, f, rnd);
  MPFR_ASSERTN (MPFR_IS_FP (r) );
  exp = mpfr_get_z_exp (z, r);
  /* FIXME: It seems that exp < 0 is always true */
  if (exp >= 0)
    mpz_mul_2exp (z, z, exp);
  else
    mpz_div_2exp (z, z, -exp);
  mpfr_clear (r);
}
