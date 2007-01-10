/* mpfr_set_exp - set the exponent of a floating-point number

Copyright 2002, 2003, 2004, 2006, 2007 Free Software Foundation, Inc.

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

int
mpfr_set_exp (mpfr_ptr x, mp_exp_t exponent)
{
  if (exponent >= __gmpfr_emin && exponent <= __gmpfr_emax)
    {
      MPFR_EXP(x) = exponent; /* do not use MPFR_SET_EXP of course... */
      return 0;
    }
  else
    {
      return 1;
    }
}
