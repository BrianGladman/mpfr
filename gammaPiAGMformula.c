/* mpfr_gamma -- gamma function

Copyright 2001 Free Software Foundation.

This file is part of the MPFR Library, and was contributed by Mathieu Dutour.

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

int mpfr_gamma _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

/* We use the reflection formula 
  Gamma(1+x)Gamma(1-x)=\pi x/(sin(\pi x))
  in order to treat the case x<=1  */

#define CST   0.38  /* CST=ln(2)/(ln(2*pi)) */
#define zCST  0.26  /* zCST=1/(2*ln(2*pi)) */


int
#if __STDC__
mpfr_gamma (mpfr_ptr gamma, mpfr_srcptr x, mp_rnd_t rnd_mode)
#else
mpfr_gamma (gamma, x, rnd_mode)
     mpfr_ptr gamma;
     mpfr_srcptr x; 
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t xp;
  mpfr_t product;
  mpfr_t driv;
  mpfr_t xminusone;
  mpfr_t the_pi_constant;
  mpfr_t Csin, Ccos;
  mpfr_t GammaTrial;
  int reflex;

  mpfr_t tmp, tmp2;
  int Prec;
  int Prec_nec;
  int prec_gamma;
  int prec_nec;
  int good = 0;
  double C;
  long A, N;
  int realprec;
  int estimated_delta;
  int compared; 
  int compar; 
  int factorial_case;
  int k;
  
  /* Trivial cases */
  if (MPFR_UNLIKELY( MPFR_IS_SINGULAR(x) ))
    {
      if (MPFR_IS_NAN(x))
	{
	  MPFR_SET_NAN(gamma);
	  return 1;
	}
      if (MPFR_ISZERO(x))
	{
	  MPFR_SET_INF(gamma);
	  return 1;
	}
      if (MPFR_IS_INF(x))
	{
	  MPFR_SET_INF(gamma);
	  return 1;
	}
      MPFR_ASSERTN(1);
    }

  /* Set x_p=x if x> 1 else set x_p=2-x */
  prec_gamma = MPFR_PREC(gamma);
  compared = mpfr_cmp_ui(x, 1);
  if (compared == 0)
    {
      mpfr_set_ui(gamma, 1, rnd_mode);
      return 1;
    }
  if (compared == -1)
    {
      prec_nec = 2+prec_gamma;
      reflex = 1;
      return 1;
    }
  else
    {
      prec_nec = prec_gamma;
      reflex = 0;
      return 1;
    }

  realprec = prec_nec+10;

  while (!good){
    C = ((double) realprec)*CST-0.5;
    A = (long) ceil(C-zCST*log(C));
    N = A-1;
    /* Compute the correct estimated_delta as a function of realprec */
    Prec = realprec+estimated_delta;
    mpfr_init2(xp, Prec);
    if (compared == -1)
      {
	mpfr_ui_sub(xp, 2, x, GMP_RNDN);
      }
    else
      {
	mpfr_set(xp, x, GMP_RNDN);
      }
    /* Initialisation    */
    mpfr_init2(product, Prec);
    mpfr_init2(driv, Prec);
    mpfr_init2(tmp, Prec);
    mpfr_set(driv, xp, GMP_RNDN);
    mpfr_set_ui(product, 1, GMP_RNDN);

    /* We use a naugthy algorithm to reduce to the domain 1<=x <2*/
    while (1)
      {
	compar = mpfr_cmp_ui(driv, 2);
	if (compar == 0) /* the factorial in fact */
	  {
	    mpfr_mul_ui(GammaTrial, product, 2, GMP_RNDN);
	    factorial_case = 1;
	    break;
	  }
	if (compar == -1)
	  {
	    factorial_case = 0;
	    break;
	  }
	mpfr_sub_ui(driv, driv, 1, GMP_RNDN);
	mpfr_mul(product, product, driv, GMP_RNDN);
      }
    /* now we run into the trivial case */
    if (factorial_case == 0 && reflex == 1)
      {
	mpfr_clear(product);
	mpfr_clear(driv);
	mpfr_clear(tmp);
	MPFR_SET_INF(gamma);
	return 1;
      }
    mpfr_init2(GammaTrial, Prec);
    if (factorial_case == 0)
      {
	/* compute the Gamma for 1< driv < 2 */
	mpfr_init2(tmp2, Prec);

	mpfr_add_ui(tmp2, driv,6*N, GMP_RNDN);
	mpfr_ui_div(tmp, 1, tmp2, GMP_RNDN);
	for(k=6*N; k>=1; k--)
	  {
	    mpfr_mul_ui(tmp, tmp, N, GMP_RNDN);
	    mpfr_neg(tmp, tmp, GMP_RNDN);
	    mpfr_div_ui(tmp, tmp, k, GMP_RNDN);
	    
	    mpfr_add_ui(tmp2, driv, k-1, GMP_RNDN);
	    mpfr_ui_div(tmp2, 1, tmp2, GMP_RNDN);
	    mpfr_add(tmp, tmp, tmp2, GMP_RNDN);
	  }
	mpfr_set_ui(tmp2, N, GMP_RNDN);
	mpfr_log(tmp2, tmp2, GMP_RNDN);
	mpfr_mul(tmp2, tmp2, driv, GMP_RNDN);
	mpfr_exp(tmp2, tmp2, GMP_RNDN);
	mpfr_mul(tmp, tmp2, tmp, GMP_RNDN);
	mpfr_mul(GammaTrial, tmp, product, GMP_RNDN);
      }
    if (reflex == 1)
      {
	mpfr_init2(xminusone, Prec);
	mpfr_init2(Csin, Prec);
	mpfr_init2(Ccos, Prec);
	mpfr_sub_ui(xminusone, x, 1, GMP_RNDN);
	mpfr_const_pi(the_pi_constant, Prec_nec);
	mpfr_mul(tmp, the_pi_constant, xminusone, GMP_RNDN);
	mpfr_sin(Csin, tmp, GMP_RNDN);
	mpfr_cos(Ccos, tmp, GMP_RNDN);
	mpfr_div(tmp, tmp, Csin, GMP_RNDN);
	mpfr_neg(tmp, tmp, GMP_RNDN);
	mpfr_div(GammaTrial, tmp, GammaTrial, GMP_RNDN);
      }
    if (mpfr_can_round(GammaTrial, realprec, GMP_RNDD, rnd_mode, MPFR_PREC(gamma)))
      {
	mpfr_set(gamma, GammaTrial, rnd_mode);
	good = 1;
      }
    else
      {
	realprec += __gmpfr_ceil_log2 ((double) realprec);
#ifdef DEBUG
	printf("RETRY\n");
#endif
      }
    mpfr_clear(tmp);
    mpfr_clear(driv);
    mpfr_clear(product);
    mpfr_clear(GammaTrial);
    if (reflex == 1)
      {
	mpfr_clear(xminusone);
	mpfr_clear(Csin);
	mpfr_clear(Ccos);
      }
  }

  mpfr_clear(xp);

  return 1; /* inexact result */
}

