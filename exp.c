#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* #define DEBUG */

#define LOG2 0.69314718055994528622 /* log(2) rounded to zero on 53 bits */

/* use Brent's formula exp(x) = (1+r+r^2/2!+r^3/3!+...)^(2^K)*2^n
   where x = n*log(2)+(2^K)*r
   number of operations = O(K+prec(r)/K)
*/
int mpfr_exp(mpfr_ptr y, mpfr_srcptr x, unsigned char rnd_mode) 
{
  int n, expx, K, precy, q, k, l, expr, err; mp_limb_t xp;
  mpfr_t r, s, t;

  if (FLAG_NAN(x)) { SET_NAN(y); return 1; }

  expx = EXP(x);
  precy = PREC(y);
#ifdef DEBUG
  printf("EXP(x)=%d\n",expx);
#endif

  /* if x > (2^31-1)*ln(2), then exp(x) > 2^(2^31-1) i.e. gives +infinity */
  if (expx > 30) {
    if (SIGN(x)>0) { printf("+infinity"); return 1; }
    else { SET_ZERO(y); return 1; }
  }

  /* if x < 2^(-precy), then exp(x) i.e. gives 1 +/- 1 ulp(1) */
  if (expx < -precy) { int signx = SIGN(x);
    mpfr_set_ui(y, 1, rnd_mode);
    if (signx>0 && rnd_mode==GMP_RNDU) mpfr_add_one_ulp(y);
    else if (signx<0 && (rnd_mode==GMP_RNDD || rnd_mode==GMP_RNDZ)) 
      mpfr_sub_one_ulp(y);
    return 1; }

  n = (int) floor(mpfr_get_d(x)/LOG2);

  /* K = (int) log( (double) precy + 2 ); */
  K = (int) sqrt( (double) precy );
  l = (precy-1)/K + 1;
  err = K + (int) ceil(log(2.0*(double)l+18.0)/LOG2);
  /* add K extra bits, i.e. failure probability <= 1/2^K = O(1/precy) */
  q = precy + err + K + 3;
  mpfr_init2(r, q); mpfr_init2(s, q); mpfr_init2(t, q);
  do {
#ifdef DEBUG
  printf("n=%d K=%d l=%d q=%d\n",n,K,l,q);
#endif

  mpfr_log2(r, GMP_RNDZ);
#ifdef DEBUG
  printf("n=%d r=",n); mpfr_print_raw(r); putchar('\n');
#endif
  mpfr_mul_ui(r, r, (n<0) ? -n : n, GMP_RNDZ); /* r = n*log(2) */
  if (n<0) mpfr_neg(r, r, GMP_RNDD);

#ifdef DEBUG
  printf("x=%1.20e\n",mpfr_get_d(x));
  printf("r=%1.20e\n",mpfr_get_d(r));
  printf(" ="); mpfr_print_raw(r); putchar('\n');
#endif
  mpfr_sub(r, x, r, GMP_RNDU);
#ifdef DEBUG
  printf("x-r=%1.20e\n",mpfr_get_d(r));
  if (SIGN(r)<0) { fprintf(stderr,"Error in mpfr_exp: r<0\n"); exit(1); }
#endif
  mpfr_div_2exp(r, r, K, GMP_RNDU); /* r = (x-n*log(2))/2^K */
  mpfr_set_ui(s, 1, GMP_RNDU);
  mpfr_set_ui(t, 1, GMP_RNDU);

  l = 1; expr = EXP(r);
  do {
    mpfr_mul(t, t, r, GMP_RNDU);
    mpfr_div_ui(t, t, l, GMP_RNDU);
    mpfr_add(s, s, t, GMP_RNDU);
#ifdef DEBUG
    printf("l=%d t=%1.20e\n",l,mpfr_get_d(t)); 
    printf("s=%1.20e\n",mpfr_get_d(s)); 
#endif
    l++;
  } while (EXP(t)+expr > -q);
#ifdef DEBUG
  fprintf(stderr,"l=%d q=%d (K+l)*q^2=%1.3e\n", l, q, (K+l)*(double)q*q);
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
     fprintf(stderr,"not enough precision, use %d\n", q+BITS_PER_MP_LIMB);
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

