/* mpfr_setmin -- minimum representable floating-point number (raw version)

Copyright 2002 Free Software Foundation.
Contributed by the Spaces project, INRIA Lorraine.

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

#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

/* Note: the flags are not cleared and the current sign is kept. */

void mpfr_setmin (mpfr_ptr x)
{
  mp_size_t xn;
  mp_limb_t *xp;

  MPFR_EXP(x) = __mpfr_emin;
  xn = (MPFR_PREC(x) - 1) / BITS_PER_MP_LIMB;
  xp = MPFR_MANT(x);
  xp[xn] = MPFR_LIMB_HIGHBIT;
  MPN_ZERO(xp, xn);
}
