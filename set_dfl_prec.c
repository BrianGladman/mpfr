/* mpfr_set_default_prec, mpfr_get_default_prec -- set/get default precision

Copyright 1999, 2000, 2001, 2004, 2005 Free Software Foundation, Inc.

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

/* default is IEEE double precision, i.e. 53 bits */
mp_prec_t MPFR_THREAD_ATTR __gmpfr_default_fp_bit_precision \
  = IEEE_DBL_MANT_DIG;

void
mpfr_set_default_prec (mp_prec_t prec)
{
  MPFR_ASSERTN (prec >= MPFR_PREC_MIN && prec <= MPFR_PREC_MAX);
  __gmpfr_default_fp_bit_precision = prec;
}

#undef mpfr_get_default_prec
mp_prec_t
mpfr_get_default_prec (void)
{
  return __gmpfr_default_fp_bit_precision;
}
