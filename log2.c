#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

mpfr_t _mpfr_log2; /* stored value of log(2) with rnd_mode=GMP_RNDZ */
int _mpfr_log2_prec=0; /* precision of stored value */

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
mpfr_log2(mpfr_ptr x, unsigned char rnd_mode)
#else
mpfr_log2(x, rnd_mode) mpfr_ptr x; unsigned char rnd_mode;
#endif
{
  int N, oldN, k, precx; mpz_t s, t, u;

  precx = PREC(x);

  /* has stored value enough precision ? */
  if (precx <= _mpfr_log2_prec) {
    if (rnd_mode==GMP_RNDZ || rnd_mode==GMP_RNDD ||
	mpfr_can_round(_mpfr_log2, _mpfr_log2_prec, GMP_RNDZ, rnd_mode, precx))
      {
	mpfr_set(x, _mpfr_log2, rnd_mode); return; 
      }
  }

  /* need to recompute */
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

  /* stored computed value */
  if (_mpfr_log2_prec==0) mpfr_init2(_mpfr_log2, precx);
  else mpfr_set_prec(_mpfr_log2, precx);
  mpfr_set(_mpfr_log2, x, GMP_RNDZ);
  _mpfr_log2_prec=precx;

  mpz_clear(s); mpz_clear(t); mpz_clear(u);
}
