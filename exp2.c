/* mpfr_exp -- exponential of a floating-point number

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
#include "mpfr.h"

int mpfr_exp2_aux  (mpfr_t, mpfr_t, mpfr_t, int);
int mpfr_exp2_aux2 (mpfr_t, mpfr_t, mpfr_t, int);

/* #define DEBUG */

#define LOG2 0.69314718055994528622 /* log(2) rounded to zero on 53 bits */

/* use Brent's formula exp(x) = (1+r+r^2/2!+r^3/3!+...)^(2^K)*2^n
   where x = n*log(2)+(2^K)*r
   together with Brent-Kung O(t^(1/2)) algorithm for the evaluation of
   power series. The resulting complexity is O(n^(1/3)*M(n)).
*/
int 
#if __STDC__
mpfr_exp2(mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode) 
#else
mpfr_exp2(y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mp_rnd_t rnd_mode;
#endif
{
  int n, expx, K, precy, q, k, l, err; 
  mpfr_t r, s, t;

  if (FLAG_NAN(x)) { SET_NAN(y); return 1; }
  if (!NOTZERO(x)) { mpfr_set_ui(y, 1, GMP_RNDN); return 0; }

  expx = EXP(x);
  precy = PREC(y);
#ifdef DEBUG
  printf("EXP(x)=%d\n",expx);
#endif

  /* if x > (2^31-1)*ln(2), then exp(x) > 2^(2^31-1) i.e. gives +infinity */
  if (expx > 30) {
    if (MPFR_SIGN(x)>0) { printf("+infinity"); return 1; }
    else { SET_ZERO(y); return 1; }
  }

  /* if x < 2^(-precy), then exp(x) i.e. gives 1 +/- 1 ulp(1) */
  if (expx < -precy) { int signx = MPFR_SIGN(x);
    mpfr_set_ui(y, 1, rnd_mode);
    if (signx>0 && rnd_mode==GMP_RNDU) mpfr_add_one_ulp(y);
    else if (signx<0 && (rnd_mode==GMP_RNDD || rnd_mode==GMP_RNDZ)) 
      mpfr_sub_one_ulp(y);
    return 1; }

  n = (int) floor(mpfr_get_d(x)/LOG2);

  K = (int) (precy<500) ? sqrt( (double) precy )
    : pow( 8.0 * (double) precy, 1.0/3.0);
  l = (precy-1)/K + 1;
  err = K + (int) ceil(log(2.0*(double)l+18.0)/LOG2);
  /* add K extra bits, i.e. failure probability <= 1/2^K = O(1/precy) */
  q = precy + err + K + 3;
  mpfr_init2(r, q); mpfr_init2(s, q); mpfr_init2(t, q);
  /* the algorithm consists in computing an upper bound of exp(x) using
     a precision of q bits, and see if we can round to PREC(y) taking
     into account the maximal error. Otherwise we increase q. */
  do {
#ifdef DEBUG
  printf("n=%d K=%d l=%d q=%d\n",n,K,l,q);
#endif

  /* if n<0, we have to get an upper bound of log(2)
     in order to get an upper bound of r = x-n*log(2) */
  mpfr_const_log2(s, (n>=0) ? GMP_RNDZ : GMP_RNDU);
#ifdef DEBUG
  printf("n=%d log(2)=",n); mpfr_print_raw(s); putchar('\n');
#endif
  mpfr_mul_ui(r, s, (n<0) ? -n : n, (n>=0) ? GMP_RNDZ : GMP_RNDU); 
  if (n<0) mpfr_neg(r, r, GMP_RNDD);
  /* r = floor(n*log(2)) */

#ifdef DEBUG
  printf("x=%1.20e\n",mpfr_get_d(x));
  printf(" ="); mpfr_print_raw(x); putchar('\n');
  printf("r=%1.20e\n",mpfr_get_d(r));
  printf(" ="); mpfr_print_raw(r); putchar('\n');
#endif
  mpfr_sub(r, x, r, GMP_RNDU);
  if (MPFR_SIGN(r)<0) { /* initial approximation n was too large */
    n--;
    mpfr_mul_ui(r, s, (n<0) ? -n : n, GMP_RNDZ);
    if (n<0) mpfr_neg(r, r, GMP_RNDD);
    mpfr_sub(r, x, r, GMP_RNDU);
  }
#ifdef DEBUG
  printf("x-r=%1.20e\n",mpfr_get_d(r));
  printf(" ="); mpfr_print_raw(r); putchar('\n');
  if (MPFR_SIGN(r)<0) { fprintf(stderr,"Error in mpfr_exp: r<0\n"); exit(1); }
#endif
  mpfr_div_2exp(r, r, K, GMP_RNDU); /* r = (x-n*log(2))/2^K */

  /* s <- 1 + r/1! + r^2/2! + ... + r^l/l! */
  l = (precy<500) ? mpfr_exp2_aux(s, r, t, q) /* naive method */
    : mpfr_exp2_aux2(s, r, t, q); /* Brent/Kung method */

#ifdef DEBUG
  printf("l=%d q=%d (K+l)*q^2=%1.3e\n", l, q, (K+l)*(double)q*q);
#endif

  /* add 2 ulp to take into account rest of summation */
  mpfr_add_one_ulp(s);
  mpfr_add_one_ulp(s);
  
  for (k=0;k<K;k++) {
    mpfr_mul(s, s, s, GMP_RNDU);
#ifdef DEBUG
    printf("k=%d s=%1.20e\n",k,mpfr_get_d(s)); 
#endif
  }

  if (n>0) mpfr_mul_2exp(s, s, n, GMP_RNDU);
  else mpfr_div_2exp(s, s, -n, GMP_RNDU);

  /* error is at most 2^K*(2l+18) ulp */
  l = 2*l+17; k=0; while (l) { k++; l >>= 1; }
  /* now k = ceil(log(2l+18)/log(2)) */
  K += k;
#ifdef DEBUG
    printf("after mult. by 2^n:\n");
    if (EXP(s)>-1024) printf("s=%1.20e\n",mpfr_get_d(s)); 
    printf(" ="); mpfr_print_raw(s); putchar('\n');
    printf("err=%d bits\n", K);
#endif

  l = mpfr_can_round(s, q-K, GMP_RNDU, rnd_mode, precy);
  if (l==0) {
#ifdef DEBUG
     printf("not enough precision, use %d\n", q+BITS_PER_MP_LIMB);
     printf("q=%d q-K=%d precy=%d\n",q,q-K,precy);
#endif
     q += BITS_PER_MP_LIMB;
     mpfr_set_prec(r, q); mpfr_set_prec(s, q); mpfr_set_prec(t, q);
  }
  } while (l==0);

  mpfr_set(y, s, rnd_mode);

  mpfr_clear(r); mpfr_clear(s); mpfr_clear(t);
  return 1;
}

/* s <- 1 + r/1! + r^2/2! + ... + r^l/l! while EXP(r^l/l!)+EXPR(r)>-q
   using naive method with O(l) multiplications
   and auxiliary variable t.
   Return l
*/
int mpfr_exp2_aux(mpfr_t s, mpfr_t r, mpfr_t t, int q)
{
  int expr, l;

  mpfr_set_ui(s, 1, GMP_RNDU);
  mpfr_set_ui(t, 1, GMP_RNDU);

  l = 1; expr = EXP(r);
  do {
    mpfr_mul(t, t, r, GMP_RNDU);
    mpfr_div_ui(t, t, l, GMP_RNDU);
    mpfr_add(s, s, t, GMP_RNDU);
    l++;
  } while (EXP(t)+expr > -q);

  return l;
}

#define MY_INIT(x, p, s) { \
   PREC(x) = p; \
   MANT(x) = (mp_ptr) TMP_ALLOC(s*BYTES_PER_MP_LIMB); \
   SIZE(x) = s; \
   EXP(x) = 0; }

/* s <- 1 + r/1! + r^2/2! + ... + r^l/l! while EXP(r^l/l!)+EXPR(r)>-q
   using Brent/Kung method with O(sqrt(l)) multiplications
   and auxiliary variable t.
   Return l.
*/
int mpfr_exp2_aux2(mpfr_t s, mpfr_t r, mpfr_t t, int q)
{
  int expr, l, m, i, sizer; mpfr_t *R; mp_limb_t c;
  TMP_DECL(marker);

  mpfr_set_ui(s, 0, GMP_RNDU);
  expr = EXP(r);

  /* estimate value of l */
  l = q / (-expr);
  m = (int) sqrt(2.0 * (double) l);
  TMP_MARK(marker);
  R = (mpfr_t*) TMP_ALLOC((m+1)*sizeof(mpfr_t)); /* R[i] stands for r^i */
  sizer = (PREC(r)+mp_bits_per_limb-1)/mp_bits_per_limb;
  for (i=0;i<=m;i++) MY_INIT(R[i], PREC(r), sizer);
  mpfr_mul(R[2], r, r, GMP_RNDU);
  for (i=3;i<=m;i++) mpfr_mul(R[i], R[i-1], r, GMP_RNDU);
  mpfr_set_ui(R[0], 1, GMP_RNDU); /* R[0]=1 */
  mpfr_set_ui(R[1], 1, GMP_RNDU); /* R[1] contains r^l/l! */

  l = 0; 
  do {
    mpfr_set(t, R[m-1], GMP_RNDU);
    /* computes t = 1 + r/(l+1) + ... + r^(m-1)*l!/(l+m-1)! */
    for (i=m-2;i>1;i--) {
      mpfr_div_ui(t, t, l+i+1, GMP_RNDU);
      mpfr_add(t, t, R[i], GMP_RNDU);
    }
    mpfr_div_ui(t, t, l+2, GMP_RNDU);
    mpfr_add(t, t, r, GMP_RNDU);
    mpfr_div_ui(t, t, l+1, GMP_RNDU);
    mpfr_add(t, t, R[0], GMP_RNDU);
    
    /* now multiplies t by r^l/l! and adds to s */
    mpfr_mul(t, t, R[1], GMP_RNDU);
    mpfr_add(s, s, t, GMP_RNDU);

    /* updates R[1], should be done using binary splitting too */
    mpfr_mul(R[1], R[1], R[m], GMP_RNDU);
    for (i=1,c=1;i<=m;i++) {
      if (l+i > ~((mp_limb_t)0)/c) {
	mpfr_div_ui(R[1], R[1], c, GMP_RNDU);
	c = l+i;
      }
      else c *= l+i;
    }
    if (c != 1) mpfr_div_ui(R[1], R[1], c, GMP_RNDU);
    l+=m;
  } while (EXP(R[1]) > -q);

  TMP_FREE(marker);
  return l;
}
