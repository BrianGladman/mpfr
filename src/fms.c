/* mpfr_fms -- Floating multiply-subtract

Copyright 2001-2002, 2004, 2006-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

#include "mpfr-impl.h"

/* The fused-multiply-subtract (fms) of x, y and z is defined by:
   fms(x,y,z)= x*y - z
   Note: this is neither in IEEE 754-2008, nor in LIA-2, but both the
   PowerPC and the Itanium define fms as x*y - z.
*/
int
mpfr_fms (mpfr_ptr s, mpfr_srcptr x, mpfr_srcptr y, mpfr_srcptr z,
          mpfr_rnd_t rnd_mode)
{
  mpfr_t minus_z;

  /* Warning! If s == z (reuse of z for the destination), s and minus_z
     will be different pointers in mpfr_fma, while they actually share
     their significand. So mpfr_fma must either assume a possible reuse
     of z (at least for the significand) in *any* case or compare the
     pointers to the significands for reuse detection, if needed. And
     be careful with sharing that propagates to called functions. */
  MPFR_TMP_INIT_NEG (minus_z, z);
  return mpfr_fma (s, x, y, minus_z, rnd_mode);
}
