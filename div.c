#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* q <- n/d using Goldschmidt's iteration
   returns 0 iff result is exact
 */
int mpfr_div(Q, n, d, rnd_mode) 
mpfr_ptr Q; mpfr_srcptr n, d; unsigned char rnd_mode;
{
  mpfr_t eps, tmp, one; mp_limb_t cc; mpfr_ptr q;
  int expd, i, prec, precq, sh, guard, err, maxprec, exact=0;
  TMP_DECL(marker);

  if (FLAG_NAN(n) || FLAG_NAN(d)) { SET_NAN(Q); return 1; }

  if (!NOTZERO(n)) { SET_ZERO(Q); return 0; }
  if (!NOTZERO(d)) { fprintf(stderr, "division by zero\n"); exit(1); }

  if (Q==n || Q==d) {
    TMP_MARK(marker);
    q = (mpfr_ptr) TMP_ALLOC(sizeof(__mpfr_struct));
    mpfr_init2(q, PREC(Q));
  }
  else q=Q;

  prec = precq = PREC(q);
  /* maxprec is the maximum number of consecutive 0's or 1's in the quotient */
  maxprec = PREC(n); if (PREC(d)>maxprec) maxprec=PREC(d);
  for (i=0;i<2;i++)
    prec = precq + (int) ceil(log(2.0*ceil(log((double)prec)/log(2.0))+7.0)/
			   log(2.0));
  err = prec-precq; /* the error is at most 2^err ulp */
  /* adjust to use complete limbs: the following formula guarantees we get
     at least one guard bit */
  prec = (prec/mp_bits_per_limb)*mp_bits_per_limb;
  do {
    prec += mp_bits_per_limb;

  mpfr_set_prec(q, prec, GMP_RNDZ);
  mpfr_set(q, n, GMP_RNDZ);
  mpfr_init2(eps, prec); mpfr_init2(tmp, prec); mpfr_init2(one, prec);
  expd = EXP(d);
  if (mpfr_cmp_si_2exp(d, SIGN(d), expd-1)==0) {
    /* d is an exact power of two */
    if (--expd>=0) mpfr_div_2exp(q, n, expd, rnd_mode);
    else mpfr_mul_2exp(q, n, -expd, rnd_mode);
    if (SIGN(d)<0) mpfr_neg(q, q, rnd_mode);
    return exact;
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
    return exact;
  }
  mpfr_mul_2exp(eps, tmp, -expd, GMP_RNDZ);
  for (;;) {
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
       if (prec>precq+maxprec+err) {
	 rnd_mode=GMP_RNDN; cc=1; exact=0; /* result is exact */
       }
     }
  } while (cc==0);
  mpfr_round(q, rnd_mode, precq);
  mpfr_mul_2exp(q, q, -expd, GMP_RNDZ);
  mpfr_set_prec(q, precq, rnd_mode);
  mpfr_clear(eps); mpfr_clear(tmp); mpfr_clear(one);
  if (Q==n || Q==d) {
    mpfr_set(Q, q, rnd_mode);
    mpfr_clear(q);
    TMP_FREE(marker);
  }
  return exact;
}
