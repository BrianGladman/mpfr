#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* #define DEBUG */

void mpfr_sqrt(mpfr_ptr X, mpfr_srcptr a, unsigned char rnd_mode)
{
  int p, q, err, i, e, n; mpfr_t t, u; mpfr_ptr x;

#ifdef DEBUG
  printf("enter mpfr_sqrt, a=%1.20e, rnd=%d\n",mpfr_get_d(a), rnd_mode);
  printf("a="); mpfr_print_raw(a); putchar('\n');
#endif
  /* use Newton's iteration x[n+1] = 1/2*(x[n]+a/x[n]),
     the error e[n] = x[n]-sqrt(a) satisfies e[n+1] <= e[n]/2/sqrt(a) */
  if (FLAG_NAN(a) || SIGN(a)<0) { SET_NAN(X); return; }
  if (X==a) {
    x = (mpfr_ptr) (*_mp_allocate_func) (sizeof(__mpfr_struct));
    mpfr_init2(x, PREC(X));
  }
  else x=X;
  e = EXP(a)/2; if (2*e<EXP(a)) e++;
#ifdef DEBUG
  printf("e=%d\n",e);
#endif
  /* now 2^(2*e-2) <= a <= 2^(2*e) i.e. 1/4 <= a/2^(2e) <= 1 */
  q = p = PREC(x);
  for (i=0; i<3; i++)
    q = p + (int) ceil(log(4.0*ceil(log((double)q)/log(2.0))+2.0)/log(2.0));
  err = q-p; /* the error is at most 2^err ulp */
  q = (q/mp_bits_per_limb)*mp_bits_per_limb; /* adjust to entire limb */
  mpfr_init(t); mpfr_init(u);
  do {
    q += mp_bits_per_limb;
    if (q>3*p+mp_bits_per_limb) {
      /* try to detect exact roots */
      mpfr_mul(t, x, x, GMP_RNDU);
      if (mpfr_cmp(t, a)==0) break; /* exact root since x>=sqrt(a) */
      fprintf(stderr, "no convergence in mpfr_sqrt for a=%1.20e, rnd=%d\n",
	      mpfr_get_d(a), rnd_mode); exit(1);
    }
#ifdef DEBUG
    printf("prec=%d q=%d err=%d\n",p,q,err);
#endif
    mpfr_set_prec(t, q, GMP_RNDU); 
    mpfr_set_prec(x, q, GMP_RNDU); 
    mpfr_set_prec(u, q, GMP_RNDU);
    mpfr_set_ui(x, 1, GMP_RNDU);
    EXP(x) += e;
    n = (int) ceil(log((double) q)/log(2.0));
    for (i=0; i<n; i++) {
      mpfr_div(t, a, x, GMP_RNDU);
      mpfr_add(u, x, t, GMP_RNDU);
      mpfr_div_2exp(x, u, 1, GMP_RNDU);
#ifdef DEBUG
      printf("i=%d t=%1.20e u=%1.20e x=%1.20e\n",i,mpfr_get_d(t),mpfr_get_d(u),
	     mpfr_get_d(x));
      printf("t="); mpfr_print_raw(t); putchar('\n');
      printf("u="); mpfr_print_raw(u); putchar('\n');
      printf("x="); mpfr_print_raw(x); putchar('\n');
#endif
    }
  } while (mpfr_can_round(x, q-err, GMP_RNDU, rnd_mode, p)==0);
  mpfr_round(x, rnd_mode, p);
  mpfr_clear(t); mpfr_clear(u);
  if (X==a) {
    mpfr_set(X, x, rnd_mode);
    mpfr_clear(x);
  }
}


