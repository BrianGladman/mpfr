/* mpfr_pow -- power function x^y

Copyright 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* return non zero iff x^y is exact.
   Assumes x and y are ordinary numbers (neither NaN nor Inf),
   and y is not zero.
   Assumes the case x < 0 and y not an integer does not happen.
*/
static int
mpfr_pow_is_exact (mpfr_srcptr x, mpfr_srcptr y)
{
  mp_exp_t d;
  unsigned long i, c;
  mp_limb_t *yp;
  mp_size_t ysize;

  if (mpfr_sgn (y) < 0)
    {
      mp_exp_t b;
      mpfr_t z;
      int res;
      if (mpfr_cmp_si_2exp (x, MPFR_INT_SIGN (x), MPFR_GET_EXP (x) - 1) != 0)
        return 0; /* x is not a power of two */
      /* now x = 2^b, so x^y = 2^(b*y) is exact whenever b*y is an integer */
      b = MPFR_GET_EXP (x) - 1; /* x = 2^b */
      mpfr_init2 (z, MPFR_PREC(y) + BITS_PER_MP_LIMB);
      mpfr_mul_si (z, y, b, GMP_RNDN); /* exact */
      res = mpfr_integer_p (z);
      mpfr_clear (z);
      return res;
    }

  /* compute d such that y = c*2^d with c odd integer */
  ysize = 1 + (MPFR_PREC(y) - 1) / BITS_PER_MP_LIMB;
  d = MPFR_GET_EXP (y) - ysize * BITS_PER_MP_LIMB;
  /* since y is not zero, necessarily one of the mantissa limbs is not zero,
     thus we can simply loop until we find a non zero limb */
  yp = MPFR_MANT(y);
  for (i = 0; yp[i] == 0; i++, d += BITS_PER_MP_LIMB);
  /* now yp[i] is not zero */
  count_trailing_zeros (c, yp[i]);
  d += c;

  if (d < 0)
    {
      mpz_t a;
      mp_exp_t b;

      mpz_init (a);
      b = mpfr_get_z_exp (a, x); /* x = a * 2^b */
      c = mpz_scan1 (a, 0);
      mpz_div_2exp (a, a, c);
      b += c;
      /* now a is odd */
      while (d != 0)
        {
          /* a * 2^b is a square iff (i)  a is a square when b is even
                                     (ii) 2*a is a square when b is odd */
          if (b % 2 != 0)
            {
              mpz_mul_2exp (a, a, 1); /* 2*a */
              b --;
            }
          if (mpz_perfect_square_p (a))
            {
              d ++;
              mpz_sqrt (a, a);
              b = b / 2;
            }
          else
            {
              mpz_clear (a);
              return 0;
            }
        }
      mpz_clear (a);
    }

    return 1;
}

#if 0
/* Return 1 if y is an odd integer, -1 if an even integer,
   0 otherwise (not an integer).
   Assumes y is not zero.
*/
static int
is_odd_even (mpfr_srcptr y)
{
  mp_exp_t expo;
  mp_prec_t prec;
  mp_size_t yn;
  mp_limb_t *yp;
  int res;

  MPFR_ASSERTD(MPFR_IS_FP(y));
  MPFR_ASSERTD(MPFR_NOTZERO(y));

  expo = MPFR_GET_EXP (y);
  if (expo <= 0)
    return 0;  /* |y| < 1 and not 0 */

  prec = MPFR_PREC(y);
  if ((mpfr_prec_t) expo > prec)
    return -1;  /* y is a multiple of 2^(expo-prec), thus an even integer */

  /* 0 < expo <= prec */

  yn = (prec - 1) / BITS_PER_MP_LIMB;  /* index of last limb */
  yn -= (mp_size_t) (expo / BITS_PER_MP_LIMB);
  MPFR_ASSERTN(yn >= 0);
  /* now the index of the last limb containing bits of the fractional part */

  yp = MPFR_MANT(y);
  res = 1; /* default is odd */
  /* test the bit of weight 0 */
  if (expo % BITS_PER_MP_LIMB == 0)
    {
      if (yp[yn])
        return 0; /* not an integer */
      if ((yp[yn+1] & 1) == 0)
        res = -1; /* maybe an even integer */
    }
  else
    {
      if (yp[yn] << (expo % BITS_PER_MP_LIMB))
        return 0; /* not an integer */
      if (yp[yn] << ((expo % BITS_PER_MP_LIMB) - 1) != MPFR_LIMB_HIGHBIT)
        res = -1; /* maybe an even integer */
    }
  while (--yn >= 0)
    if (yp[yn] != 0)
      return 0; /* not an integer */
  return res;
}
#endif

/* Return 1 if y is an odd integer, 0 otherwise. */
static int
is_odd (mpfr_srcptr y)
{
  mp_exp_t expo;
  mp_prec_t prec;
  mp_size_t yn;
  mp_limb_t *yp;

  /* NAN, INF or ZERO are not allowed */
  MPFR_ASSERTD (!MPFR_IS_SINGULAR (y));

  expo = MPFR_GET_EXP (y);
  if (expo <= 0)
    return 0;  /* |y| < 1 and not 0 */

  prec = MPFR_PREC(y);
  if ((mpfr_prec_t) expo > prec)
    return 0;  /* y is a multiple of 2^(expo-prec), thus not odd */

  /* 0 < expo <= prec */

  yn = (prec - 1) / BITS_PER_MP_LIMB;  /* index of last limb */
  yn -= (mp_size_t) (expo / BITS_PER_MP_LIMB);
  MPFR_ASSERTN(yn >= 0);
  /* now the index of the last limb containing bits of the fractional part */

  yp = MPFR_MANT(y);
  if (expo % BITS_PER_MP_LIMB == 0 ? (yp[yn+1] & 1) == 0 || yp[yn] != 0
      : yp[yn] << ((expo % BITS_PER_MP_LIMB) - 1) != MPFR_LIMB_HIGHBIT)
    return 0;
  while (--yn >= 0)
    if (yp[yn] != 0)
      return 0;
  return 1;
}

/* The computation of z = pow(x,y) is done by
   z = exp(y * log(x)) = x^y
   For the special cases, see Section F.9.4.4 of the C standard:
     _ pow(±0, y) = ±inf for y an odd integer < 0.
     _ pow(±0, y) = +inf for y < 0 and not an odd integer.
     _ pow(±0, y) = ±0 for y an odd integer > 0.
     _ pow(±0, y) = +0 for y > 0 and not an odd integer.
     _ pow(-1, ±inf) = 1.
     _ pow(+1, y) = 1 for any y, even a NaN.
     _ pow(x, ±0) = 1 for any x, even a NaN.
     _ pow(x, y) = NaN for finite x < 0 and finite non-integer y.
     _ pow(x, -inf) = +inf for |x| < 1.
     _ pow(x, -inf) = +0 for |x| > 1.
     _ pow(x, +inf) = +0 for |x| < 1.
     _ pow(x, +inf) = +inf for |x| > 1.
     _ pow(-inf, y) = -0 for y an odd integer < 0.
     _ pow(-inf, y) = +0 for y < 0 and not an odd integer.
     _ pow(-inf, y) = -inf for y an odd integer > 0.
     _ pow(-inf, y) = +inf for y > 0 and not an odd integer.
     _ pow(+inf, y) = +0 for y < 0.
     _ pow(+inf, y) = +inf for y > 0. */
int
mpfr_pow (mpfr_ptr z, mpfr_srcptr x, mpfr_srcptr y, mp_rnd_t rnd_mode)
{
  int inexact = 1;
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_LOG_FUNC (("x[%#R]=%R y[%#R]=%R rnd=%d", x, x, y, y, rnd_mode),
		 ("z[%#R]=%R inexact=%d", z, z, inexact));

  if (MPFR_ARE_SINGULAR (x, y))
    {
      /* pow(x, 0) returns 1 for any x, even a NaN. */
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (y)))
	return mpfr_set (z, __gmpfr_one, rnd_mode);
      else if (MPFR_IS_NAN (x))
	{
	  MPFR_SET_NAN (z);
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_NAN (y))
	{
	  /* pow(+1, NaN) returns 1. */
	  if (mpfr_cmp (x, __gmpfr_one) == 0)
	    return mpfr_set (z, __gmpfr_one, rnd_mode);
	  MPFR_SET_NAN (z);
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_INF (y))
	{
	  if (MPFR_IS_INF (x))
	    {
	      if (MPFR_IS_POS (y))
		MPFR_SET_INF (z);
	      else
		MPFR_SET_ZERO (z);
	      MPFR_SET_POS (z);
	      MPFR_RET (0);
	    }
	  else
	    {
	      int cmp;
	      cmp = mpfr_cmpabs (x, __gmpfr_one) * MPFR_INT_SIGN (y);
	      MPFR_SET_POS (z);
	      if (cmp > 0)
		{
		  /* Return +inf. */
		  MPFR_SET_INF (z);
		  MPFR_RET (0);
		}
	      else if (cmp < 0)
		{
		  /* Return +0. */
		  MPFR_SET_ZERO (z);
		  MPFR_RET (0);
		}
	      else
		{
		  /* Return 1. */
		  return mpfr_set (z, __gmpfr_one, rnd_mode);
		}
	    }
	}
      else if (MPFR_IS_INF (x))
	{
	  int negative;
	  /* Determine the sign now, in case y and z are the same object */
	  negative = MPFR_IS_NEG (x) && is_odd (y);
	  if (MPFR_IS_POS (y))
	    MPFR_SET_INF (z);
	  else
	    MPFR_SET_ZERO (z);
	  if (negative)
	    MPFR_SET_NEG (z);
	  else
	    MPFR_SET_POS (z);
	  MPFR_RET (0);
	}
      else
	{
	  int negative;
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
	  /* Determine the sign now, in case y and z are the same object */
	  negative = MPFR_IS_NEG(x) && is_odd (y);
	  if (MPFR_IS_NEG (y))
	    MPFR_SET_INF (z);
	  else
	    MPFR_SET_ZERO (z);
	  if (negative)
	    MPFR_SET_NEG (z);
	  else
	    MPFR_SET_POS (z);
	  MPFR_RET (0);
	}
    }

  if (mpfr_cmp (x, __gmpfr_one) == 0) /* 1^y is always 1 */
    return mpfr_set (z, __gmpfr_one, rnd_mode);

  /* detect overflows: |x^y| >= 2^EMAX when (EXP(x)-1) * y >= EMAX for y > 0,
                                       or   EXP(x) * y     >= EMAX for y < 0 */
  {
    double exy;
    int negative;

    exy = (double) (mpfr_sgn (y) > 0) ? MPFR_EXP(x) - 1 : MPFR_EXP(x);
    exy *= mpfr_get_d (y, GMP_RNDZ);
    if (exy >= (double) __gmpfr_emax)
      {
        negative = MPFR_SIGN(x) < 0 && is_odd (y);
        return mpfr_overflow (z, rnd_mode, negative ? -1 : 1);
      }
  }

  if (mpfr_integer_p (y))
    {
      mpz_t zi;

      mpz_init (zi);
      mpfr_get_z (zi, y, GMP_RNDN);
      inexact = mpfr_pow_z (z, x, zi, rnd_mode);
      mpz_clear (zi);
      return inexact;
    }

  /* x^y for x<0 and y not an integer is not defined */
  if (MPFR_IS_NEG (x))
    {
      MPFR_SET_NAN (z);
      MPFR_RET_NAN;
    }

  MPFR_SAVE_EXPO_MARK (expo);

  /* General case */
  {
    /* Declaration of the intermediary variable */
    mpfr_t t;
    int loop = 0;
    /* Declaration of the size variable */
    mp_prec_t Nz = MPFR_PREC(z);               /* target precision */
    mp_prec_t Nt;                              /* working precision */
    mp_exp_t err, exp_te;                      /* error */
    MPFR_ZIV_DECL (ziv_loop);

    /* compute the precision of intermediary variable */
    /* the optimal number of bits : see algorithms.tex */
    Nt = Nz + 5 + MPFR_INT_CEIL_LOG2 (Nz);

    /* initialise of intermediary variable */
    mpfr_init2 (t, Nt);

    MPFR_ZIV_INIT (ziv_loop, Nt);
    for (;;)
      {
        loop ++;

        /* compute exp(y*ln(x)) */
        mpfr_log (t, x, GMP_RNDU);               /* ln(x) */
        mpfr_mul (t, y, t, GMP_RNDU);            /* y*ln(x) */
	exp_te = MPFR_GET_EXP (t);               /* FIXME: May overflow */
        mpfr_exp (t, t, GMP_RNDN);               /* exp(y*ln(x))*/
                                                 /* FIXME: May overflow */

	/* estimate of the error -- see pow function in algorithms.tex.
           The error on t is at most 1/2 + 3*2^(exp_te+1) ulps, which is
           <= 2^(exp_te+3) for exp_te >= -1, and <= 2 ulps for exp_te <= -2 */
        err = (exp_te >= -1) ? Nt - (exp_te + 3) : Nt - 1;
        if (MPFR_LIKELY (MPFR_CAN_ROUND (t, err, Nz, rnd_mode)))
	  break;

        /* check exact power */
        if (loop == 1)
          {
            if (mpfr_pow_is_exact (x, y))
	      {
		inexact = 0;
		break;
	      }
	  }

        /* reactualisation of the precision */
	MPFR_ZIV_NEXT (ziv_loop, Nt);
        mpfr_set_prec (t, Nt);
      }
    MPFR_ZIV_FREE (ziv_loop);

    if (inexact)
      inexact = mpfr_set (z, t, rnd_mode);
    else /* result is exact: round to nearest and return inexact=0 */
      mpfr_set (z, t, GMP_RNDN);

    mpfr_clear (t);
  }

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (z, inexact, rnd_mode);
}
