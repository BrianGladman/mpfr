/* mpfr_reldiff -- compute relative difference of two floating-point numbers.

Copyright (C) 2000 Free Software Foundation.

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
  if (MPFR_IS_NAN(b) || MPFR_IS_NAN(c))
    { MPFR_CLEAR_FLAGS(a); MPFR_SET_NAN(a); return; }
  if (MPFR_IS_INF(b)) 
    { 
      if (MPFR_IS_INF(c) && (MPFR_SIGN(c) == MPFR_SIGN(b)))
	{ MPFR_CLEAR_FLAGS(a); MPFR_SET_ZERO(a); return; }
      else
	{ MPFR_CLEAR_FLAGS(a); MPFR_SET_NAN(a); return; }
    }

  if (MPFR_IS_INF(c)) 
    {
      if (MPFR_SIGN(a) != MPFR_SIGN(b)) { MPFR_CHANGE_SIGN(a); }
      MPFR_CLEAR_FLAGS(a);
      MPFR_SET_INF(a);
    }

  if (MPFR_IS_ZERO(b)) /* reldiff = abs(c)/c = sign(c) */
    /*    TODO: faire preciser la SEMANTIQUE DE CE FOUTOIR. */
    mpfr_set_ui(a, MPFR_SIGN(c), rnd_mode);
  else {
    mpfr_sub(a, b, c, rnd_mode);
    mpfr_abs(a, a, rnd_mode); /* for compatibility with MPF */
    mpfr_div(a, a, b, rnd_mode);
  }
}

