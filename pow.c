/* mpfr_pow_ui, mpfr_ui_pow_ui -- compute the power of a floating-point
                                  number or machine integer

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
#include "mpfr.h"
#include "mpfr-impl.h"

/* sets x to y^n, and return 0 if exact, non-zero otherwise */
int
#if __STDC__
mpfr_pow_ui (mpfr_ptr x, mpfr_srcptr y, unsigned long int n, mp_rnd_t rnd)
#else
mpfr_pow_ui (x, y, n, rnd)
     mpfr_ptr x;
     mpfr_srcptr y;
     unsigned long int n;
     mp_rnd_t rnd;
#endif
{
  long int i, err;
  unsigned long m;
  mpfr_t res;
  mp_prec_t prec;
  int inexact;
  mp_rnd_t rnd1;

  if (MPFR_IS_NAN(y))
    {
      MPFR_SET_NAN(x);
      MPFR_RET_NAN;
    }

  MPFR_CLEAR_NAN(x);

  if (n == 0) /* x^0 = 1 for any x */
    {
      mpfr_set_ui (x, 1, rnd);
      return 0;
    }

  if (MPFR_IS_INF(y)) 
    {
      /* Inf^n = Inf, (-Inf)^n = Inf for n even, -Inf for n odd */
      if ((MPFR_SIGN(y) < 0) && (n % 2 == 1))
	{
	  if (MPFR_SIGN(x) > 0)
	    MPFR_CHANGE_SIGN(x);
	}
      else if (MPFR_SIGN(x) < 0)
	MPFR_CHANGE_SIGN(x);
      MPFR_SET_INF(x);
      return 0;
    }
  
  MPFR_CLEAR_INF(x);

  mpfr_init (res);

  prec = MPFR_PREC(x);

  rnd1 = (MPFR_SIGN(y) > 0) ? GMP_RNDU : GMP_RNDD; /* away */

  do
    {
      prec += 3;
      for (m=n, i=0; m; i++, m>>=1, prec++);
      mpfr_set_prec (res, prec);
      inexact = mpfr_set (res, y, rnd1);
      err = 1;
      /* now 2^(i-1) <= n < 2^i */
      for (i-=2; i>=0; i--)
	{
	  if (mpfr_mul (res, res, res, GMP_RNDU))
	    inexact = 1;
	  err++;
	  if (n & (1 << i))
	    if (mpfr_mul (res, res, y, rnd1))
	      inexact = 1;
	}
      err = prec - err;
      if (err < 0)
	err = 0;
    }
  while (inexact && (mpfr_can_round (res, err,
	  (MPFR_SIGN(res) > 0) ? GMP_RNDU : GMP_RNDD, rnd, MPFR_PREC(x)) == 0));

  if (mpfr_set (x, res, rnd))
    inexact = 1;

  mpfr_clear (res);

  return inexact;
}

/* sets x to y^n, and return 0 if exact, non-zero otherwise */

int
#if __STDC__
mpfr_ui_pow_ui (mpfr_ptr x, unsigned long int y, unsigned long int n,
		     mp_rnd_t rnd)
#else
mpfr_ui_pow_ui (x, y, n, rnd)
     mpfr_ptr x;
     unsigned long int y;
     unsigned long int n;
     mp_rnd_t rnd;
#endif
{
  long int i, err;
  unsigned long m;
  mpfr_t res;
  mp_prec_t prec;
  int inexact;

  MPFR_CLEAR_FLAGS(x);

  if (n == 0) /* x^0 = 1 for any x */
    {
      mpfr_set_ui (x, 1, rnd);
      return 0;
    }

  mpfr_init (res);

  prec = MPFR_PREC(x);

  do
    {
      prec += 3;
      for (m=n, i=0; m; i++, m>>=1, prec++);
      mpfr_set_prec (res, prec);
      inexact = mpfr_set_ui (res, y, GMP_RNDU);
      err = 1;
      /* now 2^(i-1) <= n < 2^i */
      for (i-=2; i>=0; i--)
	{
	  if (mpfr_mul (res, res, res, GMP_RNDU))
	    inexact = 1;
	  err++;
	  if (n & (1 << i))
	    if (mpfr_mul_ui (res, res, y, GMP_RNDU))
	      inexact = 1;
	}
      err = prec - err;
      if (err < 0)
	err = 0;
    }
  while (inexact && (mpfr_can_round (res, err,
	  (MPFR_SIGN(res) > 0) ? GMP_RNDU : GMP_RNDD, rnd, MPFR_PREC(x)) == 0));

  if (mpfr_set (x, res, rnd))
    inexact = 1;

  mpfr_clear (res);

  return inexact;
}
