#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void 
#if __STDC__
mpfr_sqrt(mpfr_ptr X, mpfr_srcptr a, unsigned char rnd_mode)
#else
mpfr_sqrt(X, a, rnd_mode)
     mpfr_ptr X;
     mpfr_srcptr a;
     unsigned char rnd_mode;
#endif
{
  int p, q, err, i, e, n; mpfr_t t, u; mpfr_ptr x;

  /* use Newton's iteration x[n+1] = 1/2*(x[n]+a/x[n]),
     the error e[n] = x[n]-sqrt(a) satisfies e[n+1] <= e[n]/2/sqrt(a) */
  if (FLAG_NAN(a) || SIGN(a)<0) { SET_NAN(X); return; }

  x = (mpfr_ptr) (*_mp_allocate_func) (sizeof(__mpfr_struct));
  mpfr_init2(x, PREC(X));
  mpfr_set(x, X, rnd_mode);

  e = EXP(a)/2; if (2*e<EXP(a)) e++;
  /* now 2^(2*e-2) <= a < 2^(2*e) i.e. 1/4 <= a/2^(2e) < 1 */
  q = p = PREC(x);
  for (i=0; i<3; i++)
    q = p + (int) ceil(log(4.0*ceil(log((double)q)/log(2.0))+2.0)/log(2.0));
  err = q-p; /* the error is at most 2^err ulp */
  q = (q/mp_bits_per_limb)*mp_bits_per_limb; /* adjust to entire limb */
  mpfr_init2(t, q+mp_bits_per_limb); mpfr_init2(u, q+mp_bits_per_limb);
  do {
    q += mp_bits_per_limb;
    if (q>3*p+mp_bits_per_limb) {
      /* try to detect exact roots */
      mpfr_round(x, rnd_mode, p);
      mpfr_mul(t, x, x, rnd_mode);
      if (mpfr_cmp(t, a)==0) goto youpi; /* exact root */
      fprintf(stderr, "no convergence in mpfr_sqrt for a=%1.20e, rnd=%d\n",
	      mpfr_get_d(a), rnd_mode); exit(1);
    }
    mpfr_set_prec(t, q); 
    mpfr_set_prec(x, q); 
    mpfr_set_prec(u, q);
    mpfr_set_ui(x, 1, GMP_RNDU);
    EXP(x) += e;
    n = (int) ceil(log((double) q)/log(2.0));
    for (i=0; i<n; i++) {
      mpfr_div(t, a, x, GMP_RNDU);
      mpfr_add(u, x, t, GMP_RNDU);
      mpfr_div_2exp(x, u, 1, GMP_RNDU);
    }
  } while (mpfr_can_round(x, q-err, GMP_RNDU, rnd_mode, p)==0);
  mpfr_round(x, rnd_mode, p);
 youpi:
  mpfr_clear(t); mpfr_clear(u);
  mpfr_set(X, x, rnd_mode);
  mpfr_clear(x); (*_mp_free_func)(x, sizeof(__mpfr_struct)); 
}


