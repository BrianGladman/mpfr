/* Exception flags and utilities.

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

#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

unsigned int __mpfr_flags = 0;

mp_exp_t __mpfr_emin = MPFR_EMIN_DEFAULT;
mp_exp_t __mpfr_emax = MPFR_EMAX_DEFAULT;

#undef mpfr_get_emin

mp_exp_t
#if __STDC__
mpfr_get_emin (void)
#else
mpfr_get_emin ()
#endif
{
  return __mpfr_emin;
}

#undef mpfr_set_emin

int
#if __STDC__
mpfr_set_emin (mp_exp_t exponent)
#else
mpfr_set_emin ()
#endif
{
  if (exponent >= MPFR_EMIN_MIN && exponent <= MPFR_EMIN_MAX)
  {
    __mpfr_emin = exponent;
    return 0;
  }
  else
  {
    return 1;
  }
}

#undef mpfr_get_emax

mp_exp_t
#if __STDC__
mpfr_get_emax (void)
#else
mpfr_get_emax ()
#endif
{
  return __mpfr_emax;
}

#undef mpfr_set_emax

int
#if __STDC__
mpfr_set_emax (mp_exp_t exponent)
#else
mpfr_set_emax ()
#endif
{
  if (exponent >= MPFR_EMAX_MIN && exponent <= MPFR_EMAX_MAX)
  {
    __mpfr_emax = exponent;
    return 0;
  }
  else
  {
    return 1;
  }
}

#undef mpfr_clear_flags

void
#if __STDC__
mpfr_clear_flags (void)
#else
mpfr_clear_flags ()
#endif
{
  __mpfr_flags = 0;
}

#undef mpfr_clear_underflow

void
#if __STDC__
mpfr_clear_underflow (void)
#else
mpfr_clear_underflow ()
#endif
{
  __mpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_UNDERFLOW;
}

#undef mpfr_clear_overflow

void
#if __STDC__
mpfr_clear_overflow (void)
#else
mpfr_clear_overflow ()
#endif
{
  __mpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_OVERFLOW;
}

#undef mpfr_check_range

int
#if __STDC__
mpfr_check_range (mpfr_ptr x, mp_rnd_t rnd_mode)
#else
mpfr_check_range ()
#endif
{
  if (MPFR_IS_FP(x) && MPFR_NOTZERO(x))
  { /* x is a non-zero FP */
    mp_exp_t exp = MPFR_EXP(x);
    if (exp < __mpfr_emin)
    {
      mpfr_set_underflow(x, rnd_mode, MPFR_SIGN(x));
      return -1;
    }
    if (exp > __mpfr_emax)
    {
      mpfr_set_overflow(x, rnd_mode, MPFR_SIGN(x));
      return 1;
    }
  }
  return 0;
}

#undef mpfr_underflow_p

int
#if __STDC__
mpfr_underflow_p (void)
#else
mpfr_underflow_p ()
#endif
{
  return __mpfr_flags & MPFR_FLAGS_UNDERFLOW;
}

#undef mpfr_overflow_p

int
#if __STDC__
mpfr_overflow_p (void)
#else
mpfr_overflow_p ()
#endif
{
  return __mpfr_flags & MPFR_FLAGS_OVERFLOW;
}

#undef mpfr_set_underflow

void
#if __STDC__
mpfr_set_underflow (mpfr_ptr x, mp_rnd_t rnd_mode, int sign)
#else
mpfr_set_underflow ()
#endif
{
  MPFR_CLEAR_FLAGS(x);
  if ((rnd_mode == GMP_RNDU && sign > 0)
   || (rnd_mode == GMP_RNDD && sign < 0))
  {
    mp_size_t xn;
    mp_limb_t *xp;
    MPFR_EXP(x) = __mpfr_emin;
    xn = (MPFR_PREC(x)-1)/BITS_PER_MP_LIMB;
    xp = MPFR_MANT(x);
    xp[xn] = MP_LIMB_T_HIGHBIT;
    MPN_ZERO(xp, xn);
  }
  else
  {
    MPFR_SET_ZERO(x);
  }
  if (MPFR_SIGN(x) != sign) { MPFR_CHANGE_SIGN(x); }
  __mpfr_flags |= MPFR_FLAGS_UNDERFLOW;
}

#undef mpfr_set_overflow

void
#if __STDC__
mpfr_set_overflow (mpfr_ptr x, mp_rnd_t rnd_mode, int sign)
#else
mpfr_set_overflow ()
#endif
{
  MPFR_CLEAR_FLAGS(x);
  if ((rnd_mode == GMP_RNDU && sign < 0)
   || (rnd_mode == GMP_RNDD && sign > 0))
  {
    mp_size_t xn, i;
    mp_limb_t *xp;
    MPFR_EXP(x) = __mpfr_emax;
    xn = (MPFR_PREC(x)-1)/BITS_PER_MP_LIMB;
    xp = MPFR_MANT(x);
    for (i = 0; i <= xn; i++)
      xp[i] = MP_LIMB_T_MAX;
  }
  else
  {
    MPFR_SET_INF(x);
  }
  if (MPFR_SIGN(x) != sign) { MPFR_CHANGE_SIGN(x); }
  __mpfr_flags |= MPFR_FLAGS_OVERFLOW;
}
