/* mpfr_init2 -- initialize a floating-point number with given precision

Copyright 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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

#include <limits.h>
#include <stdlib.h>
#include "mpfr-impl.h"

void
mpfr_init2 (mpfr_ptr x, mp_prec_t p)
{
  mp_size_t xsize;
  mp_ptr tmp;

  /* Check if we can represent the number of limbs
   * associated to the maximum of mpfr_prec_t*/
  MPFR_ASSERTN( MP_SIZE_T_MAX >= (MPFR_PREC_MAX/BYTES_PER_MP_LIMB) );

  /* Check for correct BITS_PER_MP_LIMB and BYTES_PER_MP_LIMB */
  MPFR_ASSERTN( BITS_PER_MP_LIMB == BYTES_PER_MP_LIMB * CHAR_BIT
		&& sizeof(mp_limb_t) == BYTES_PER_MP_LIMB );

  /* Check for correct EXP NAN, ZERO & INF in both mpfr.h and in mpfr-impl.h */
  MPFR_ASSERTN( __MPFR_EXP_NAN  == MPFR_EXP_NAN  );
  MPFR_ASSERTN( __MPFR_EXP_ZERO == MPFR_EXP_ZERO );
  MPFR_ASSERTN( __MPFR_EXP_INF  == MPFR_EXP_INF  );

  MPFR_ASSERTN( MPFR_EMAX_MAX <= (MPFR_EXP_MAX >> 1)  );
  MPFR_ASSERTN( MPFR_EMIN_MIN >= -(MPFR_EXP_MAX >> 1) );

  /* p=1 is not allowed since the rounding to nearest even rule requires at
     least two bits of mantissa: the neighbours of 3/2 are 1*2^0 and 1*2^1,
     which both have an odd mantissa */
  MPFR_ASSERTN(p >= MPFR_PREC_MIN && p <= MPFR_PREC_MAX);

  xsize = (mp_size_t) ((p - 1) / BITS_PER_MP_LIMB) + 1;
  tmp   = (mp_ptr) (*__gmp_allocate_func)(MPFR_MALLOC_SIZE(xsize));

  MPFR_PREC(x) = p;                /* Set prec */
  MPFR_EXP (x) = MPFR_EXP_INVALID; /* make sure that the exp field has a
                                      valid value in the C point of view */
  MPFR_SET_POS(x);                 /* Set a sign */
  MPFR_SET_MANT_PTR(x, tmp);       /* Set Mantissa ptr */
  MPFR_SET_ALLOC_SIZE(x, xsize);   /* Fix alloc size of Mantissa */
  MPFR_SET_NAN(x);                 /* initializes to NaN */
}
