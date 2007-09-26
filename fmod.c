/* mpfr_fmod -- compute the floating-point remainder of x/y 

Copyright 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

/* adapted from remquo.c */

# include "mpfr-impl.h"

int
mpfr_fmod (mpfr_ptr rem, mpfr_srcptr x, mpfr_srcptr y, mp_rnd_t rnd)
{
  mp_exp_t ex, ey;
  int compare, inex, sign, signx = MPFR_SIGN(x);
  mpz_t mx, my, r;

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x) || MPFR_IS_SINGULAR (y)))
    {
      if (MPFR_IS_NAN (x) || MPFR_IS_NAN (y) || MPFR_IS_INF (x)
          || MPFR_IS_ZERO (y))
        {
          MPFR_SET_NAN (rem);
          MPFR_RET_NAN;
        }
      else /* either y is Inf and x is 0 or non-special,
              or x is 0 and y is non-special. */
	return mpfr_set (rem, x, rnd);
    }

  /* now neither x nor y is NaN, Inf or zero */

  mpz_init (mx);
  mpz_init (my);
  mpz_init (r);

  ex = mpfr_get_z_exp (mx, x); /* x = mx*2^ex */
  ey = mpfr_get_z_exp (my, y); /* y = my*2^ey */

  /* to get rid of sign problems, we compute it separately:
     quo(-x,-y) = quo(x,y), rem(-x,-y) = -rem(x,y)
     quo(-x,y) = -quo(x,y), rem(-x,y)  = -rem(x,y)
     thus quo = sign(x/y)*quo(|x|,|y|), rem = sign(x)*rem(|x|,|y|) */
  sign = (signx == MPFR_SIGN(y)) ? 1 : -1;
  mpz_abs (mx, mx);
  mpz_abs (my, my);

  /* divide my by 2^k if possible to make operations mod my easier */
  {
    unsigned long k = mpz_scan1 (my, 0);
    ey += k;
    mpz_div_2exp (my, my, k);
  }

  if (ex <= ey)
    {
      /* q = x/y = mx/(my*2^(ey-ex)) */
      mpz_mul_2exp (my, my, ey - ex); /* divide mx by my*2^(ey-ex) */
      mpz_fdiv_qr (mx, r, mx, my); /* 0 <= |r| <= |my|, r has the same
                                      sign as my */
    }
  else /* ex > ey */
    {
      /* Let X = mx*2^(ex-ey) and Y = my. Then both X and Y are integers.
         Assume X = R mod Y, then x = X*2^ey = R*2^ey mod (Y*2^ey=y).
         To be able to perform the rounding, we need the least significant
         bit of the quotient, i.e., one more bit in the remainder, which is
         obtained by dividing by 2Y.
      */
      mpz_mul_2exp (my, my, 1); /* 2Y */
      mpz_set_ui (r, 2);
      mpz_powm_ui (r, r, ex - ey, my); /* 2^(ex-ey) mod my */
      mpz_mul (r, r, mx);
      mpz_mod (r, r, my);

      /* now 0 <= r < 2Y */
      mpz_div_2exp (my, my, 1); /* back to Y */
    }

  if (mpz_cmp_ui (r, 0) == 0)
    inex = mpfr_set_ui (rem, 0, GMP_RNDN);
  else
    {
      /* FIXME: the comparison 2*r < my could be done more efficiently
         at the mpn level */
      mpz_mul_2exp (r, r, 1);
      compare = mpz_cmpabs (r, my);
      mpz_div_2exp (r, r, 1);
      /* if compare != 0, we need to subtract my to r */
      if (compare > 0)
          mpz_sub (r, r, my);
      inex = mpfr_set_z (rem, r, rnd);
      /* if ex > ey, rem should be multiplied by 2^ey, else by 2^ex */
      MPFR_EXP(rem) += (ex > ey) ? ey : ex;
    }

  /* take into account sign of x */
  if (signx < 0)
    {
      mpfr_neg (rem, rem, GMP_RNDN);
      inex = -inex;
    }

  mpz_clear (mx);
  mpz_clear (my);
  mpz_clear (r);

  return inex;
}
