/* mpfr_pow_si -- power function x^y with y a signed int

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

#include "mpfr-impl.h"

 /* The computation of y=pow(x,z) is done by

    y=pow_ui(x,z) if z>0
  else
    y=1/pow_ui(x,z) if z<0
 */

int
mpfr_pow_si (mpfr_ptr y, mpfr_srcptr x, long int n, mp_rnd_t rnd_mode)
{
  if (n >= 0)
    return mpfr_pow_ui (y, x, n, rnd_mode);
  else
    {
      if (MPFR_UNLIKELY( MPFR_IS_SINGULAR(x) ))
	{
	  if (MPFR_IS_NAN(x))
	    {
	      MPFR_SET_NAN(y);
	      MPFR_RET_NAN;
	    }
	  else if (MPFR_IS_INF(x))
	    {
	      MPFR_SET_ZERO(y);
	      if (MPFR_IS_POS(x) || ((unsigned) n & 1) == 0)
		MPFR_SET_POS(y);
	      else
		MPFR_SET_NEG(y);
	      MPFR_RET(0);
	    }
	  else /* x is zero */
	    {
              MPFR_ASSERTD(MPFR_IS_ZERO(x));
	      MPFR_SET_INF(y);
	      if (MPFR_IS_POS(x) || ((unsigned) n & 1) == 0)
		MPFR_SET_POS(y);
	      else
		MPFR_SET_NEG(y);
	      MPFR_RET(0);
	    }
	}
      MPFR_CLEAR_FLAGS(y);

      /* detect exact powers: x^(-n) is exact iff x is a power of 2 */
      if (mpfr_cmp_si_2exp (x, MPFR_SIGN(x), MPFR_EXP(x) - 1) == 0)
        {
          mp_exp_t expx = MPFR_EXP(x); /* warning: x and y may be the same
                                            variable */
          mpfr_set_si (y, (n % 2) ? MPFR_SIGN(x) : 1, rnd_mode);
          MPFR_EXP(y) += n * (expx - 1);
          return 0;
        }

      n = -n;

      /* General case */
      {
        /* Declaration of the intermediary variable */
        mpfr_t t, ti;

        /* Declaration of the size variable */
        mp_prec_t Nx = MPFR_PREC(x);   /* Precision of input variable */
        mp_prec_t Ny = MPFR_PREC(y);   /* Precision of output variable */

        mp_prec_t Nt;   /* Precision of the intermediary variable */
        long int err;   /* Precision of error */
        int inexact;

        /* compute the precision of intermediary variable */
        Nt = MAX(Nx,Ny);
        /* the optimal number of bits : see algorithms.ps */
        Nt = Nt + 3 + __gmpfr_ceil_log2 (Nt);

        mpfr_save_emin_emax ();

        /* initialise of intermediary	variable */
        mpfr_init2 (t, Nt);
        mpfr_init2 (ti, Nt);

        do
          {
            /* reactualisation of the precision */
            mpfr_set_prec (t, Nt);
            mpfr_set_prec (ti, Nt);

            /* compute 1/(x^n) n>0*/
            mpfr_pow_ui (ti, x, (unsigned long int) n, GMP_RNDN);
            mpfr_ui_div (t, 1, ti, GMP_RNDN);

            /* error estimate -- see pow function in algorithms.ps */
            err = Nt - 3;

            /* actualisation of the precision */
            Nt += 10;
          }
        while (err < 0 || (!mpfr_can_round (t, err, GMP_RNDN, GMP_RNDZ,
					    Ny + (rnd_mode == GMP_RNDN))
			   /* An overflow can occurs, producing an underflow */
			   && !MPFR_IS_ZERO(t) ));

        inexact = mpfr_set (y, t, rnd_mode);
        mpfr_clear (t);
        mpfr_clear (ti);
        mpfr_restore_emin_emax ();
        return mpfr_check_range (y, inexact, rnd_mode);
      }
    }
}




