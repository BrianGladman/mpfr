/* mpfr_div_2exp -- divide a floating-point number by a power of two

Copyright 1999, 2001, 2004, 2006 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include "mpfr-impl.h"

/* Obsolete function, use mpfr_div_2ui or mpfr_div_2si instead. */

#undef mpfr_div_2exp

int
mpfr_div_2exp (mpfr_ptr y, mpfr_srcptr x, unsigned long int n, mp_rnd_t rnd_mode)
{
  return mpfr_div_2ui (y, x, n, rnd_mode);
}
