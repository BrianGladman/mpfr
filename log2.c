/* mpfr_const_log2 -- compute natural logarithm of 2

Copyright (C) 1999 PolKA project, Inria Lorraine and Loria

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
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

mpfr_t __mpfr_const_log2; /* stored value of log(2) with rnd_mode=GMP_RNDZ */
int __mpfr_const_log2_prec=0; /* precision of stored value */


#define A
#define A1 1
#define A2 1
#undef B
#define C
#define C1 2
#define C2 1
#define NO_FACTORIAL
#undef R_IS_RATIONAL
#define GENERIC mpfr_aux_log2
#include "generic.c" 
#undef A
#undef A1
#undef A2
#undef NO_FACTORIAL
#undef GENERIC
#undef C
#undef C1
#undef C2

int
#if __STDC__
mpfr_const_aux_log2(mpfr_ptr mylog, mp_rnd_t rnd_mode)
#else
mpfr_const_aux_log2(mylog, rnd_mode) mpfr_ptr mylog; mp_rnd_t rnd_mode;
#endif
{
  int prec;
  mpfr_t tmp1, tmp2, result,tmp3; 
  mpz_t cst;
  int good = 0;
  int logn;
  int prec_i_want = PREC(mylog);
  int prec_x;
  mpz_init(cst);
  logn =  (int) ceil(log
                      ((double) PREC(mylog))
                      /log(2.0)); 
  prec_x = prec_i_want + logn;
  while (!good){
    prec = (int) ceil(log
		      ((double) (prec_x))
		      /log(2.0));  
    mpfr_init2(tmp1, prec_x);
    mpfr_init2(result, prec_x);
    mpfr_init2(tmp2, prec_x);
    mpfr_init2(tmp3, prec_x);
    mpz_set_ui(cst, 1);
    mpfr_aux_log2(tmp1, cst, 4, prec-2);
    mpfr_div_2exp(tmp1, tmp1, 4,GMP_RNDD);
    mpfr_mul_ui(tmp1, tmp1, 15,GMP_RNDD);

    mpz_set_ui(cst, 3);
    mpfr_aux_log2(tmp2, cst, 7, prec-2);
    mpfr_div_2exp(tmp2, tmp2, 7,GMP_RNDD);
    mpfr_mul_ui(tmp2, tmp2, 5*3,GMP_RNDD);
    mpfr_sub(result, tmp1, tmp2, GMP_RNDD);

    mpz_set_ui(cst, 13);
    mpfr_aux_log2(tmp3, cst, 8, prec-2);
    mpfr_div_2exp(tmp3, tmp3, 8,GMP_RNDD);
    mpfr_mul_ui(tmp3, tmp3, 3*13,GMP_RNDD);
    mpfr_sub(result, result, tmp3, GMP_RNDD);

    mpfr_clear(tmp1);
    mpfr_clear(tmp2);
    mpfr_clear(tmp3);
    if (mpfr_can_round(result, prec_x, GMP_RNDD, rnd_mode, prec_i_want)){
      mpfr_set(mylog, result, rnd_mode);
      good = 1;
    } else
      {
	prec_x += logn;
      }
    mpfr_clear(result);
  }
  mpz_clear(cst);
  return 0;
}

/* set x to log(2) rounded to precision PREC(x) with direction rnd_mode 

   use formula log(2) = sum(1/k/2^k, k=1..infinity)

   whence 2^N*log(2) = S(N) + R(N)

   where S(N) = sum(2^(N-k)/k, k=1..N-1)
   and   R(N) = sum(1/k/2^(k-N), k=N..infinity) < 2/N

   Let S'(N) = sum(floor(2^(N-k)/k), k=1..N-1)

   Then 2^N*log(2)-S'(N) <= N-1+2/N <= N for N>=2.
*/
void 
#if __STDC__
mpfr_const_log2(mpfr_ptr x, mp_rnd_t rnd_mode)
#else
mpfr_const_log2(x, rnd_mode) mpfr_ptr x; mp_rnd_t rnd_mode;
#endif
{
  int N, oldN, k, precx; mpz_t s, t, u;

  precx = PREC(x);

  /* has stored value enough precision ? */
  if (precx <= __mpfr_const_log2_prec) {
    if (rnd_mode==GMP_RNDZ || rnd_mode==GMP_RNDD ||
	mpfr_can_round(__mpfr_const_log2, __mpfr_const_log2_prec, GMP_RNDZ, 
		       rnd_mode, precx))
      {
	mpfr_set(x, __mpfr_const_log2, rnd_mode); return; 
      }
  }

  /* need to recompute */
  if (precx < 30000){ /* use nai"ve Taylor series evaluation */
     N=2;
     do {
       oldN = N;
       N = precx + (int)ceil(log((double)N)/log(2.0));
     } while (N != oldN);
     mpz_init_set_ui(s,0);
     mpz_init(u);
     mpz_init_set_ui(t,1); 
   #if 0
     /* use log(2) = sum(1/k/2^k, k=1..infinity) */
     mpz_mul_2exp(t, t, N);
     for (k=1;k<N;k++) {
       mpz_div_2exp(t, t, 1);
       mpz_fdiv_q_ui(u, t, k);
       mpz_add(s, s, u);
     }
   #else
     /* use log(2) = sum((6*k-1)/(2*k^2-k)/2^(2*k+1), k=1..infinity) */
     mpz_mul_2exp(t, t, N-1);
     for (k=1;k<N/2;k++) {
       mpz_div_2exp(t, t, 2);
       mpz_mul_ui(u, t, 6*k-1);
       mpz_fdiv_q_ui(u, u, k*(2*k-1));
       mpz_add(s, s, u);
     }
   #endif
     mpfr_set_z(x, s, rnd_mode);
     EXP(x) -= N;
     mpz_clear(s); mpz_clear(t); mpz_clear(u);
  } else
    {
      /* use binary splitting method */
      mpfr_const_aux_log2(x, rnd_mode);
    }

  /* store computed value */
  if (__mpfr_const_log2_prec==0) mpfr_init2(__mpfr_const_log2, precx);
  else mpfr_set_prec(__mpfr_const_log2, precx);
  mpfr_set(__mpfr_const_log2, x, GMP_RNDZ);
  __mpfr_const_log2_prec=precx;

}
