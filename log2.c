#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

/* set x to log(2) rounded to precision PREC(x) with direction rnd_mode 

   use formula log(2) = sum(1/k/2^k, k=1..infinity)

   whence 2^N*log(2) = S(N) + R(N)

   where S(N) = sum(2^(N-k)/k, k=1..N-1)
   and   R(N) = sum(1/k/2^(k-N), k=N..infinity) < 2/N

   Let S'(N) = sum(floor(2^(N-k)/k), k=1..N-1)

   Then 2^N*log(2)-S'(N) <= N-1+2/N <= N for N>=2.
*/
mpfr_log2(x, rnd_mode) mpfr_ptr x; unsigned char rnd_mode;
{
  int N, oldN, k; mpz_t s, t, u;

  N=2;
  do {
    oldN = N;
    N = PREC(x) + (int)ceil(log((double)N)/log(2.0));
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
}
