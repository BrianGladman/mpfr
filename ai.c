/* mpfr_ai -- Airy function Ai

Copyright 2007, 2008, 2009, 2010 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#define MPFR_SMALL_PRECISION MPFR_PREC_MIN

/* Reminder and notations:
   -----------------------

   Ai is the solution of:
        / y'' - x*y = 0
       {    Ai(0)   = 1/ ( 9^(1/3)*Gamma(2/3) )
        \  Ai'(0)   = -1/ ( 3^(1/3)*Gamma(1/3) ) 

   Series development:
       Ai(x) = sum (a_i*x^i)
             = sum (t_i)

   Recurrences:
       a_(i+3) = a_i / ((i+2)*(i+3))
       t_(i+3) = t_i * x^3 / ((i+2)*(i+3))

   Values:
       a_0 = Ai(0)  ~  0.355
       a_1 = Ai'(0) ~ -0.259
*/


/* Airy function Ai evaluated by the most naive algorithm */
int
mpfr_ai(mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd)
{
  MPFR_ZIV_DECL(loop);
  MPFR_SAVE_EXPO_DECL(expo);
  mpfr_prec_t wprec;             /* working precision */
  mpfr_prec_t prec;              /* target precision */
  mpfr_prec_t err;               /* used to estimate the evaluation error */
  mpfr_prec_t correct_bits;      /* estimates the number of correct bits*/
  unsigned long int k;
  unsigned long int cond;        /* condition number of the series */
  unsigned assumed_exponent;     /* used as a lowerbound of |EXP(Ai(x))| */
  int r;
  mpfr_t s;                      /* used to store the partial sum */
  mpfr_t ti, tip1;   /* used to store successive values of t_i */
  mpfr_t x3;                     /* used to store x^3 */
  mpfr_t tmp_sp, tmp2_sp;        /* small precision variables */
  unsigned long int x3u;         /* used to store ceil(x^3) */
  mpfr_t temp;
  int test1, test2;

  /* Logging */
  MPFR_LOG_FUNC( ("x[%#R]=%R rnd=%d", x, x, rnd), ("y[%#R]=%R", y, y) );

  /* Special cases */
  if (MPFR_UNLIKELY(MPFR_IS_NAN(x))) 
    { 
      MPFR_SET_NAN(y); 
      MPFR_RET_NAN; 
    }
  if (MPFR_UNLIKELY(MPFR_IS_INF(x))) 
    return mpfr_set_ui(y, 0, rnd);

  /* Save current exponents range */
  MPFR_SAVE_EXPO_MARK(expo);

  /* FIXME: underflow for large values of |x| ? */


  /* Set initial precision */
  /* If we compute sum(i=0,N-1, t_i), the relative error is bounded by   */
  /*       2*(4N)*2^(1-wprec)*C(|x|)/Ai(x)                               */
  /* where C(|x|) = 1 if 0<=x<=1                                         */
  /*   and C(|x|) = (1/2)*x^(-1/4)*exp(2/3 x^(3/2))  if x >= 1           */

  /* A priori, we do not know N, so we estimate it to ~ prec             */
  /* If 0<=x<=1, we estimate Ai(x) ~ 1/8                                 */
  /* if 1<=x,    we estimate Ai(x) ~ (1/4)*x^(-1/4)*exp(-2/3 * x^(3/2))  */
  /* if x<=0,    ?????                                                   */

  /* We begin with 11 guard bits */
  prec = MPFR_PREC(y)+11;
  MPFR_ZIV_INIT(loop, prec);

  /* This variable is used to store the maximal assumed exponent of     */
  /* Ai(x). More precisely, we assume that |Ai(x)| will be greater than */
  /* 2^{-assumed_exponent}.                                                   */
  assumed_exponent = 10;

  /* The working precision is heuristically chosen in order to obtain  */
  /* approximately prec correct bits in the sum. To sum up: the sum    */
  /* is stopped when the *exact* sum gives ~ prec correct bit. And     */
  /* when it is stopped, the accuracy of the computed sum, with respect*/
  /* to the exact one should be ~prec bits.                            */
  mpfr_init2(tmp_sp, MPFR_SMALL_PRECISION);
  mpfr_init2(tmp2_sp, MPFR_SMALL_PRECISION);
  mpfr_abs(tmp_sp, x, MPFR_RNDU);
  mpfr_pow_ui(tmp_sp, tmp_sp, 3, MPFR_RNDU);
  mpfr_sqrt(tmp_sp, tmp_sp, MPFR_RNDU); /* tmp_sp ~ x^3/2 */

  /* 0.96179669392597567 >~ 2/3 * log2(e). See algorithms.tex */
  mpfr_mul_d(tmp2_sp, tmp_sp, 0.96179669392597567, MPFR_RNDU);
  
  /* cond represents the number of lost bits in the evaluation of the sum */
  if(MPFR_EXP(x)<=0) 
    cond = 0;
  else
    cond = mpfr_get_ui(tmp2_sp, MPFR_RNDU) - (MPFR_EXP(x)-1)/4 - 1;

  wprec = prec + MPFR_INT_CEIL_LOG2(prec) + 5 + cond;
  if(MPFR_IS_POS(x))
    {
      if(MPFR_EXP(x)<=0) wprec += 3;
      else
	wprec += 2 + (MPFR_EXP(x)/4 + 1) + mpfr_get_ui(tmp2_sp, MPFR_RNDU);
    }
  else 
    {
      /* wprec = wprec - EXP(Ai(x)) but we do not know Ai(x) yet */
      wprec += assumed_exponent; /* We cover the case when EXP(Ai(x))>=-10 */
    }

  mpfr_init(ti);  
  mpfr_init(tip1);
  mpfr_init(temp);
  mpfr_init(x3);  
  mpfr_init(s);

  /* ZIV loop */
  for (;;) 
    {
      MPFR_LOG_MSG(("Working precision: %d\n", wprec, 0));
      mpfr_set_prec(ti, wprec);
      mpfr_set_prec(tip1, wprec);
      mpfr_set_prec(temp, wprec+4);
      mpfr_set_prec(x3, wprec);
      mpfr_set_prec(s, wprec);
      
      mpfr_sqr(x3, x, MPFR_RNDU);
      mpfr_mul(x3, x3, x, (MPFR_IS_POS(x)?MPFR_RNDU:MPFR_RNDD));  /* x3 = x^3 */
      if (MPFR_IS_NEG(x)) MPFR_CHANGE_SIGN(x3);
      x3u = mpfr_get_ui(x3, MPFR_RNDU);   /* x3u >= ceil(x^3) */
      if (MPFR_IS_NEG(x)) MPFR_CHANGE_SIGN(x3);

      
      mpfr_set_ui(ti, 9, MPFR_RNDN);
      mpfr_cbrt(ti, ti, MPFR_RNDN);
      mpfr_set_ui(temp, 2, MPFR_RNDN);
      mpfr_div_ui(temp, temp, 3, MPFR_RNDN);
      mpfr_gamma(temp, temp, MPFR_RNDN);
      mpfr_mul(ti, ti, temp, MPFR_RNDN);
      mpfr_ui_div(ti, 1, ti , MPFR_RNDN); /* ti = 1 / ( Gamma(2/3)*9^(1/3) ) */
      
      mpfr_set_ui(tip1, 3, MPFR_RNDN);
      mpfr_cbrt(tip1, tip1, MPFR_RNDN);
      mpfr_set_ui(temp, 1, MPFR_RNDN);
      mpfr_div_ui(temp, temp, 3, MPFR_RNDN);
      mpfr_gamma(temp, temp, MPFR_RNDN);
      mpfr_mul(tip1, tip1, temp, MPFR_RNDN);
      mpfr_neg(tip1, tip1, MPFR_RNDN);
      mpfr_div(tip1, x, tip1, MPFR_RNDN); /* tip1 = -x/(Gamma(1/3)*3^(1/3)) */
      
      mpfr_add(s, ti, tip1, MPFR_RNDN);
      
      
      /* Evaluation of the series */
      k = 2;
      for(;;)
	{
	  mpfr_mul(ti, ti, x3, MPFR_RNDN);
	  mpfr_mul(tip1, tip1, x3, MPFR_RNDN);
	  
	  if (k+2 <= ULONG_MAX/(k+1)) 
	    {
	      mpfr_div_ui(ti, ti, k*(k+1), MPFR_RNDN);
	      mpfr_div_ui(tip1, tip1, (k+1)*(k+2), MPFR_RNDN);
	    }
	  else 
	    {
	      mpfr_div_ui(ti, ti, k, MPFR_RNDN);
	      mpfr_div_ui(ti, ti, k+1, MPFR_RNDN);
	      mpfr_div_ui(tip1, tip1, k+1, MPFR_RNDN);
	      mpfr_div_ui(tip1, tip1, k+2, MPFR_RNDN);
	    }
	  
	  k += 3;
	  mpfr_add(s, s, ti, MPFR_RNDN);
	  mpfr_add(s, s, tip1, MPFR_RNDN);
	  
	  /* FIXME: if s==0 */
	  test1 = MPFR_IS_ZERO(ti)
	    || ( MPFR_EXP(ti) + (mp_exp_t)prec + 3 <= MPFR_EXP(s) );
	  test2 = MPFR_IS_ZERO(tip1) 
	    || ( MPFR_EXP(tip1) + (mp_exp_t)prec + 3 <= MPFR_EXP(s) );

	  if ( test1 && test2 && (x3u <= k*(k+1)/2) ) 
	    break; /* FIXME: if k*(k+1) overflows */
	}

      MPFR_LOG_MSG(("Truncation rank: %d\n", k, 0));
    
      err = 4 + MPFR_INT_CEIL_LOG2(k) + cond - MPFR_EXP(s);

      /* err is the number of bits lost due to the evaluation error */
      /* wprec-(prec+1): number of bits lost due to the approximation error */
      MPFR_LOG_MSG(("Roundoff error: %d\n", err, 0));
      MPFR_LOG_MSG(("Approxim error: %d\n", wprec-prec-1, 0));

      if (wprec < err+1) 
	correct_bits=0;
      else
	{
	  if (wprec<err+prec+1)
	    correct_bits=(unsigned)((signed)wprec-(signed)err-1);
	  else 
	    correct_bits = prec;
	}
      
      if (MPFR_LIKELY(MPFR_CAN_ROUND(s, correct_bits, MPFR_PREC(y), rnd))) 
	break;
    
      if (correct_bits==0)
	{
	  assumed_exponent *= 2;
	  MPFR_LOG_MSG(("Not a single bit correct (assumed_exponent=%d)\n",0));
	  wprec = prec + 5 + MPFR_INT_CEIL_LOG2(k) + cond + assumed_exponent;
	}
      else 
	{
	  if (correct_bits < prec)
	    { /* The precision was badly chosen */
	      MPFR_LOG_MSG(("Bad assumption on the exponent of Ai(x) (E=%d)\n",MPFR_EXP(s),0));
	      wprec = prec + err + 1;
	    }
	  else 
	    { /* We are really in a bad case of the TMD */
	      MPFR_ZIV_NEXT(loop, prec);
	      
	      /* We update wprec */
	      /* We assume that K will not be multiplied by more than 4 */
	      wprec = prec + (MPFR_INT_CEIL_LOG2(k)+2) + 5 + cond - MPFR_EXP(s);
	    }
	}
      
    } /* End of ZIV loop */

  MPFR_ZIV_FREE(loop);
  MPFR_SAVE_EXPO_FREE(expo);

  r = mpfr_set(y,s,rnd);
  
  mpfr_clear(ti);
  mpfr_clear(tip1);
  mpfr_clear(temp);
  mpfr_clear(x3);
  mpfr_clear(s);
  mpfr_clear(tmp_sp);
  mpfr_clear(tmp2_sp);
  return r;
}
