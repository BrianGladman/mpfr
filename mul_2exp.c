/* mpfr_mul_2exp -- multiply a floating-point number by a power of two

Copyright (C) 1999, 2001 Free Software Foundation, Inc.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

int
mpfr_mul_2exp (mpfr_ptr y, mpfr_srcptr x, unsigned long int n, mp_rnd_t rnd_mode)
{
  int inexact = 0;

  /* Important particular case */ 
  if (y != x)
    inexact = mpfr_set (y, x, rnd_mode);
  return ((MPFR_EXP(y) += n) > __mpfr_emax)
    ? mpfr_set_overflow (y, rnd_mode, MPFR_SIGN(y)) : inexact;
}

