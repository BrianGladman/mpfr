/* mpfr_reldiff -- compute relative difference of two floating-point numbers.

Copyright (C) 2000 PolKA project, Inria Lorraine and Loria

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

/* reldiff(b, c) = abs(b-c)/b */
void 
#if __STDC__
mpfr_reldiff(mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mp_rnd_t rnd_mode)
#else
mpfr_reldiff(a, b, c, rnd_mode) 
     mpfr_ptr a;
     mpfr_srcptr b;
     mpfr_srcptr c;
     mp_rnd_t rnd_mode;
#endif
{
  if (MPFR_IS_NAN(b) || MPFR_IS_NAN(c)) { MPFR_SET_NAN(a); return; }

  if (!MPFR_NOTZERO(b)) /* reldiff = abs(c)/c = sign(c) */
    mpfr_set_ui(a, MPFR_SIGN(c), rnd_mode);

  else {
    mpfr_sub(a, b, c, rnd_mode);
    mpfr_abs(a, a, rnd_mode); /* for compatibility with MPF */
    mpfr_div(a, a, b, rnd_mode);
  }
}

