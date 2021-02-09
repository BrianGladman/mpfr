/* mpfr_pow_uj -- compute the power of a floating-point by a uintmax_t

Copyright 2021 Free Software Foundation, Inc.
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
#define MPFR_NEED_INTMAX_H
#include "mpfr-impl.h"

#define	POW_U mpfr_pow_uj
#define MPZ_SET_U mpfr_mpz_set_uj
#define	UTYPE uintmax_t
#define LONG_TYPE 0

#define ULONG_BITS (sizeof (unsigned long) * CHAR_BIT)

/* z <- n, assuming uintmax_t is at most twice as wide as unsigned long */
static void
mpfr_mpz_set_uj (mpz_t z, uintmax_t n)
{
  uintmax_t h;

  h = n >> ULONG_BITS;

  mpz_set_ui (z, (unsigned long) h);
  MPFR_ASSERTN((h >> ULONG_BITS) == 0);
  mpz_mul_2exp (z, z, ULONG_BITS);
  mpz_add_ui (z, z, (unsigned long) n);
}

#include "pow_ui.c"
