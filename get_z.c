/* mpfr_get_z_exp -- get a multiple-precision integer from 
                     a floating-point number

Copyright 2004 Free Software Foundation, Inc.

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

void
mpz_set_fr (mpz_ptr z, mpfr_srcptr f, mp_rnd_t rnd)
{
  mpfr_t r;
  mp_exp_t exp = MPFR_EXP(f);
  
  /* if exp <= 0, then |f|<1, thus |o(f)|<=1 */
  mpfr_init2 (r, (exp <= 0) ? MPFR_PREC_MIN : exp + 1);
  mpfr_rint (r, f, rnd);
  MPFR_ASSERTN (mpfr_number_p (r) );
  exp = mpfr_get_z_exp (z, r);
  if (exp >= 0)
    mpz_mul_2exp (z, z, exp);
  else
    mpz_div_2exp (z, z, -exp);
  mpfr_clear (r);
}
