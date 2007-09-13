/* mpfr_modf -- Integral and fractionnal part.

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

#include "mpfr-impl.h"

/* Set iop to the integral part of u and fop to its fractional part */
int
mpfr_modf (mpfr_ptr iop, mpfr_ptr fop, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mp_exp_t ue;
  mp_prec_t uq;

  MPFR_ASSERTN (iop != fop);

  if ( MPFR_UNLIKELY (MPFR_IS_SINGULAR (u)) )
    {
      if (MPFR_IS_NAN (u))
	{
	  MPFR_SET_NAN (iop);
	  MPFR_SET_NAN (fop);
	  MPFR_RET_NAN;
	}
      MPFR_SET_SAME_SIGN (iop, u);
      MPFR_SET_SAME_SIGN (fop, u);
      if (MPFR_IS_INF (u))
	{
	  MPFR_SET_INF (iop);
	  MPFR_SET_ZERO (fop);
	  MPFR_RET (0);
	}
      else /* u is zero */
	{
	  MPFR_ASSERTD (MPFR_IS_ZERO (u));
	  MPFR_SET_ZERO (iop);
	  MPFR_SET_ZERO (fop);
	  MPFR_RET (0);
	}
    }

  ue = MPFR_GET_EXP (u);
  uq = MPFR_PREC (u);

  if (ue <=0)   /* 0 < |u| < 1 */
    {
      if (fop != u)      
	mpfr_set (fop, u, rnd_mode);
      MPFR_SET_SAME_SIGN (iop, u);
      MPFR_SET_ZERO (iop);
      MPFR_RET (MPFR_INT_SIGN (u) > 0 ? -2 : +2);
    }
  else if (ue >= uq) /* u has no fractional part */
    {
      int inexact;
      inexact = (iop != u)? mpfr_set (iop, u, rnd_mode) : 0;
      MPFR_SET_SAME_SIGN (fop, u);
      MPFR_SET_ZERO (fop);
      MPFR_RET (inexact);
    }
  else /* u has both integral and fractional parts */
    {
      int inexact, inexi, inexf;
      mpfr_t uf, ui;

      /* ui and uf are set with minimal but sufficient precision */    
      mpfr_init2 (ui, ue <= MPFR_PREC_MIN ? MPFR_PREC_MIN : ue);
      inexi = mpfr_trunc (ui, u);
      mpfr_init2 (uf, uq - ue <= MPFR_PREC_MIN ? MPFR_PREC_MIN : uq - ue); 
      inexf = mpfr_frac (uf, u, GMP_RNDZ);
      MPFR_ASSERTD (inexf == 0);

      /* note: while the exponent of the fractional part may be out of range, */
      /*       the int-part exponent can't since it is the same as u */
      inexf = mpfr_set (fop, uf, rnd_mode);
      inexf = mpfr_check_range (fop, inexf, rnd_mode);
      inexi = mpfr_set (iop, ui, rnd_mode);
      mpfr_clear (ui);
      mpfr_clear (uf);

      /* return value like mpfr_trunc():   */
      /* 0 iff iop and fop are exact       */
      /* -1 if u is an integer, u > iop    */
      /* +1 if u is an integer, u < iop    */
      /* -2 if u is not an integer, u > 0  */
      /* +2 if u is not an integer, u < 0  */
      inexact = inexf ? (inexi ? 2 * inexi : -2 * MPFR_INT_SIGN (u)) : (mpfr_zero_p (fop) ? inexi : 2 * inexi);
      MPFR_RET (inexact);
    }
}
