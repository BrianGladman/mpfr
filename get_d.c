/* mpfr_get_d -- convert a multiple precision floating-point number
                 to a machine double precision float

Copyright 1999, 2000, 2001, 2002 Free Software Foundation, Inc.

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

#include <float.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"
#include "mpfr-impl.h"

static double __mpfr_scale2 _PROTO ((double, int));

#define NaN (0./0.) /* ensures a machine-independent NaN */
#define Infp (1/0.)
#define Infm (-1/0.)

static double
__mpfr_scale2 (double d, int exp)
{
#if _GMP_IEEE_FLOATS
  {
    union ieee_double_extract x;

    if (exp < -2099)
      return 0.0 * d; /* 0 with the correct sign */

    x.d = d;
    if (exp >= 2047 || exp + x.s.exp >= 2047)
      {
        /* Return +-infinity */
        x.s.exp = 2047;
        x.s.manl = x.s.manh = 0;
      }
    else if (exp + x.s.exp < 1)
      {
        exp += x.s.exp;
        if (exp <= -52)
          return 0.0 * d; /* 0 with the correct sign */
        x.s.exp = 1; /* smallest exponent (biased) */
        x.d *= __mpfr_scale2(1.0, exp - 1);
      }
    else
      {
        x.s.exp += exp;
      }
    return x.d;
  }
#else
  {
    double factor;

    if (exp < 0)
      {
        factor = 0.5;
        exp = -exp;
      }
    else
      {
        factor = 2.0;
      }
    while (exp != 0)
      {
        if ((exp & 1) != 0)
          d *= factor;
        exp >>= 1;
        factor *= factor;
      }
    return r;
  }
#endif
}

double
mpfr_get_d2 (mpfr_srcptr src, mp_exp_t e)
{
  double d;
  mpfr_t tmp;
  int s, negative;

  if (MPFR_IS_NAN(src))
    return NaN;

  negative = MPFR_SIGN(src) < 0;

  if (MPFR_IS_INF(src))
    return negative ? Infm : Infp;

  if (MPFR_IS_ZERO(src))
    return negative ? -0.0 : 0.0;

  if (e < -1076)
    { /* Simulate the underflow with the current IEEE rounding mode. */
      d = DBL_MIN;
      d *= negative ? -DBL_MIN : DBL_MIN;
      /* -> double precision forced by the affectation */
    }
  else if (e > 1024)
    { /* Simulate the overflow with the current IEEE rounding mode. */
      d = DBL_MAX;
      d *= negative ? -2 : 2;
    }
  else
    {
      mp_ptr tp;
      mp_size_t np, i;
      double r;

      mpfr_save_emin_emax();

      /* Truncate src to 54 bits
       * --> rounding bit stored to the 54th bit
       * --> sticky bit stored to variable s */
      mpfr_init2(tmp, 54);
      s = mpfr_set(tmp, src, GMP_RNDZ);
      MPFR_ASSERTN(MPFR_IS_FP(tmp) && MPFR_NOTZERO(tmp));

      /* Warning: the rounding may still be incorrect in the rounding
         to the nearest mode when the result is a subnormal because of
         a double rounding (-> 53 bits -> final precision). */
      tp = MPFR_MANT(tmp);
      d = (tp[0] & (MP_LIMB_T_MAX << 11)) / MP_BASE_AS_DOUBLE;
      np = (MPFR_PREC(tmp) - 1) / BITS_PER_MP_LIMB;
      for (i = 1; i <= np; i++)
        d = (d + tp[i]) / MP_BASE_AS_DOUBLE;
      /* d is the mantissa (between 1/2 and 1) of the argument truncated
         to 53 bits */
      r = (((tp[0] >> 9) & 2) + (s != 0)) * (DBL_EPSILON / 8.0);
      d += r; /* double precision forced by the affectation */
      d = __mpfr_scale2(d, e);
      if (negative)
        d = -d;

      mpfr_clear(tmp);
      mpfr_restore_emin_emax();
    }

  return d;
}

double
mpfr_get_d (mpfr_srcptr src)
{
  return mpfr_get_d2 (src, MPFR_EXP(src));
}
