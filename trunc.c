/* mpf_trunc, mpf_floor, mpf_ceil -- Assign a float from another float while
   rounding it to an integer.

Copyright (C) 1999, 2001 Free Software Foundation, Inc.

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

#if MPFR_FLOOR
#define _MPFR_FLOOR_OR_CEIL
#define FUNC_NAME mpfr_floor
#undef MPFR_FLOOR
#define MPFR_FLOOR 1
#define MPFR_CEIL 0
#endif

#if MPFR_CEIL
#define _MPFR_FLOOR_OR_CEIL
#define FUNC_NAME mpfr_ceil
#undef MPFR_CEIL
#define MPFR_CEIL 1
#define MPFR_FLOOR 0
#endif

#if MPFR_TRUNC
#undef FUNC_NAME
#define FUNC_NAME mpfr_trunc
#endif

#ifdef _MPFR_FLOOR_OR_CEIL
static int
mpn_zero_p (mp_ptr p, mp_size_t n)
{
  mp_size_t i;

  for (i = 0; i < n; i++)
    {
      if (p[i] != 0)
	return 0;
    }

  return 1;
}
#endif

void
FUNC_NAME (mpfr_ptr r, mpfr_srcptr u)
{
  mp_ptr rp, up;
  mp_size_t rsize, usize;
#ifdef _MPFR_FLOOR_OR_CEIL
  mp_size_t ignored_n;
#endif
  mp_exp_t exp;
  int rw;
  int signu;
  mp_exp_t diff;

  if (MPFR_IS_NAN(u))
    {
      MPFR_SET_NAN(r);
      __mpfr_flags |= MPFR_FLAGS_NAN;
      return;
    }

  MPFR_CLEAR_NAN(r);
  MPFR_SET_SAME_SIGN(r, u);

  if (MPFR_IS_INF(u))
    {
      MPFR_SET_INF(r);
      return;
    }

  MPFR_CLEAR_INF(r);

  if (MPFR_IS_ZERO(u))
    {
      MPFR_SET_ZERO(r);
      return;
    }

  signu = MPFR_SIGN(u);
  rp = MPFR_MANT(r);
  exp = MPFR_EXP(u);
  rsize = MPFR_ESIZE(r);

  /* Single out the case where |u| < 1. */
  if (exp <= 0)
    {
#ifdef _MPFR_FLOOR_OR_CEIL
      if ((MPFR_FLOOR && signu < 0) || (MPFR_CEIL && signu >= 0))
	{
	  rp[rsize-1] = MP_LIMB_T_HIGHBIT;
	  MPN_ZERO(rp, rsize-1);
	  /* sign of result is that of u */
	  MPFR_EXP(r) = 1;
	  return;
	}
#endif
      MPFR_SET_ZERO(r);
      return;
    }

  usize = MPFR_ESIZE(u);

#ifdef _MPFR_FLOOR_OR_CEIL
  ignored_n = 0;
#endif
  up = MPFR_MANT(u);

  if (usize > rsize)
    {
#ifdef _MPFR_FLOOR_OR_CEIL
      ignored_n = usize - rsize;
#endif
      up += usize - rsize;
      usize = rsize;
    }

  diff = (mp_prec_t) usize * BITS_PER_MP_LIMB - exp;
  if (diff > 0)
    {
      diff /= BITS_PER_MP_LIMB;
#ifdef _MPFR_FLOOR_OR_CEIL
      ignored_n += diff;
#endif
      up += diff;
      usize -= diff;
    }

  /* number of non significant bits in low limb of r */
  MPFR_ASSERTN ((mp_prec_t) usize * BITS_PER_MP_LIMB >= exp &&
                (mp_prec_t) usize * BITS_PER_MP_LIMB - exp < BITS_PER_MP_LIMB);
  rw = (mp_prec_t) usize * BITS_PER_MP_LIMB - exp;
  MPN_ZERO(rp, rsize - usize);
  rp += rsize - usize;

#ifdef _MPFR_FLOOR_OR_CEIL
  if (((MPFR_FLOOR && signu < 0) || (MPFR_CEIL && signu >= 0))
      && (!mpn_zero_p (up - ignored_n, ignored_n)
	  || (rw && (up[0] << (BITS_PER_MP_LIMB - rw)))))
    {
      mp_limb_t cy;

      cy = mpn_add_1 (rp, up, usize, MP_LIMB_T_ONE << rw);
      if (cy != 0)
	{
          if (exp == __mpfr_emax)
            {
              mpfr_set_overflow(r, MPFR_FLOOR ? GMP_RNDD : GMP_RNDU,
                                MPFR_SIGN(r));
              return;
            }

	  /* all the bits from "1<<rw" upwards are zero */
	  rp[usize-1] = MP_LIMB_T_HIGHBIT;
	  exp++;
	}
    }
  else
#endif
    {
      if (rp != up)
	MPN_COPY (rp, up, usize);
    }

  /* Put to 0 the remaining bits */
  if (rw != 0)
    rp[0] &= ~((MP_LIMB_T_ONE << rw) - MP_LIMB_T_ONE);

  MPFR_EXP(r) = exp;
}
