#include <math.h>
#include "gmp.h"
#include "mpfr.h"

/* #define DEBUG
#define DEBUG2 */

/* q <- n/d using Goldschmidt's iteration */
void mpfr_div(q, n, d, rnd_mode) 
mpfr_ptr q; mpfr_srcptr n, d; unsigned char rnd_mode;
{
  mpfr_t eps, tmp, one; int expd, i, prec, precq, sh, guard, err;
  mp_limb_t cc;

  if (q==n || q==d) {
    fprintf(stderr, "destination equals source in mpfr_div\n");
    exit(1);
  }
#ifdef DEBUG
  printf("enter mpfr_div, prec(q)=%d n=%1.20e prec(n)=%d d=%1.20e prec(d)=%d rnd=%d\n",PREC(q),mpfr_get_d(n),PREC(n),mpfr_get_d(d),PREC(d),rnd_mode); 
  printf("n="); mpfr_print_raw(n); putchar('\n');
  printf("d="); mpfr_print_raw(d); putchar('\n');
#endif
  if ((MANT(n)[(PREC(n)-1)/mp_bits_per_limb] & 
      ((mp_limb_t)1<<(mp_bits_per_limb-1)))==0) {
    printf("Error in mpfr_div: n is not msb-normalized\n"); exit(1);
  }
  if (FLAG_NAN(n) || FLAG_NAN(d)) {
#ifdef DEBUG
    printf("dividend or divisor is NaN\n");
#endif
    SET_NAN(q); return;
  }
  prec = precq = PREC(q);
  for (i=0;i<2;i++)
    prec = precq + (int) ceil(log(2.0*ceil(log((double)prec)/log(2.0))+7.0)/
			   log(2.0));
  err = prec-precq; /* the error is at most 2^err ulp */
  prec++; /* add one bit otherwise mfpr_can_round will always fail */
  prec = prec-mp_bits_per_limb;
  do {
    prec += mp_bits_per_limb;
#ifdef DEBUG2
  printf("PREC(q)=%d prec=%d\n",precq,prec);
  printf("n=%1.20e d=%1.20e rnd=%d\n",mpfr_get_d(n),mpfr_get_d(d),rnd_mode); 
#endif
  mpfr_set_prec(q, prec, GMP_RNDZ);
  mpfr_set(q, n, GMP_RNDZ);
  mpfr_init2(eps, prec); mpfr_init2(tmp, prec); mpfr_init2(one, prec);
  expd = EXP(d);
  if (mpfr_cmp_si_2exp(d, SIGN(d), expd-1)==0) {
    /* d is an exact power of two */
    if (--expd>=0) mpfr_div_2exp(q, n, expd, rnd_mode);
    else mpfr_mul_2exp(q, n, -expd, rnd_mode);
    if (SIGN(d)<0) mpfr_neg(q, q, rnd_mode);
    return;
  }
  mpfr_set_ui(one, 1, GMP_RNDZ); 
  mpfr_mul_2exp(eps, one, expd, GMP_RNDZ); /* eps = 2^expd */
  if (SIGN(d)<0) mpfr_add(tmp, eps, d, GMP_RNDZ);
  else mpfr_sub(tmp, eps, d, GMP_RNDZ);
  i=0; while (i<ABSSIZE(tmp) && MANT(tmp)[i]==MANT(d)[i]) i++;
  if (i==ABSSIZE(tmp)) {
    /* if tmp=abs(d), then eps=2*d whence d is an exact power of two */
    PREC(q) = precq;
    mpfr_set(q, n, rnd_mode);
    EXP(q) -= expd-1;
    if (SIGN(d)<0) CHANGE_SIGN(q);
    return;
  }
  mpfr_mul_2exp(eps, tmp, -expd, GMP_RNDZ);
  for (;;) {
#ifdef DEBUG2
    printf("q=%1.20e eps=%1.20e\n",mpfr_get_d(q),mpfr_get_d(eps));
#endif
    /* q[n+1] = q[n] * (1 + e[n])
       Numerically, it's better to compute q + (q*e) than q * (1+e).
       Rounding towards zero guarantees we get a lower bound of n/d.
     */
    mpfr_mul(tmp, q, eps, GMP_RNDZ);
    mpfr_add(one, q, tmp, GMP_RNDZ);
    mpfr_set(q, one, GMP_RNDZ);
    /* e[n+1] = e[n]^2 */
    mpfr_mul(tmp, eps, eps, GMP_RNDZ);
    mpfr_set(eps, tmp, GMP_RNDZ);
    if (-EXP(eps)>prec) break; /* further iterations won't change
					   q since we round towards zero. */
  }
     if (SIGN(d)<0) CHANGE_SIGN(q);
     cc = mpfr_can_round(q, prec-err, GMP_RNDZ, rnd_mode, precq);
     if (cc==0) {
#ifdef DEBUG
       printf("not enough precision\n");
       printf("q="); mpfr_print_raw(q); putchar('\n');
#endif
       if (prec>2*precq) { printf("does not converge\n"); exit(1); }
     }
  } while (cc==0);
  mpfr_round(q, rnd_mode, precq);
  mpfr_mul_2exp(q, q, -expd, GMP_RNDZ);
  mpfr_set_prec(q, precq, rnd_mode);
  mpfr_clear(eps); mpfr_clear(tmp); mpfr_clear(one);
#ifdef DEBUG
  printf("q = %1.20e\n", mpfr_get_d(q));
  printf("n/d=%1.20e\n", mpfr_get_d(n)/mpfr_get_d(d));
#endif
}
