/* mpfr_set -- copy of a floating-point number

Copyright (C) 1999 Free Software Foundation.

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
#if __STDC__
mpfr_set4 (mpfr_ptr a, mpfr_srcptr b, mp_rnd_t rnd_mode, int signb)
#else
mpfr_set4 (a, b, rnd_mode, signb)
     mpfr_ptr a; 
     mpfr_srcptr b; 
     mp_rnd_t rnd_mode;
     int signb;
#endif
{
  int inex;

  if (MPFR_IS_NAN(b))
  {
    MPFR_CLEAR_FLAGS(a);
    MPFR_SET_NAN(a);
    MPFR_RET(0);
  }

  if (MPFR_IS_INF(b))
  { 
    MPFR_CLEAR_FLAGS(a);
    MPFR_SET_INF(a);
    inex = 0;
  }
  else
  {
    mp_limb_t *ap;
    mp_prec_t aq;
    int carry;

    MPFR_CLEAR_FLAGS(a);

    ap = MPFR_MANT(a);
    aq = MPFR_PREC(a);

    carry = mpfr_round_raw(ap, MPFR_MANT(b), MPFR_PREC(b), (signb < 0),
                           aq, rnd_mode, &inex);
    MPFR_EXP(a) = MPFR_EXP(b);

    if (carry)
    {
      mp_exp_t exp = MPFR_EXP(a);
      if (exp == __mpfr_emax)
      {
        mpfr_set_overflow(a, rnd_mode, signb);
        MPFR_RET(inex);
      }
      else
      {
        MPFR_EXP(a)++;
        ap[(MPFR_PREC(a)-1)/BITS_PER_MP_LIMB] = MP_LIMB_T_HIGHBIT;
      }
    }
  }

  if (MPFR_SIGN(a) * signb < 0) MPFR_CHANGE_SIGN(a);
  MPFR_RET(inex);
}
