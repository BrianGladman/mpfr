/* bernoulli -- internal function to compute Bernoulli numbers.

Copyright 2005-2014 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

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

#include "mpfr-impl.h"

/* assuming b[0]...b[2(n-1)] are computed, computes and stores B[2n]*(2n+1)!

   t/(exp(t)-1) = sum(B[j]*t^j/j!, j=0..infinity)
   thus t = (exp(t)-1) * sum(B[j]*t^j/j!, n=0..infinity).
   Taking the coefficient of degree n+1 > 1, we get:
   0 = sum(1/(n+1-k)!*B[k]/k!, k=0..n)
   which gives:
   B[n] = -sum(binomial(n+1,k)*B[k], k=0..n-1)/(n+1).

   Let C[n] = B[n]*(n+1)!.
   Then C[n] = -sum(binomial(n+1,k)*C[k]*n!/(k+1)!,  k=0..n-1),
   which proves that the C[n] are integers.
*/
void
mpfr_bernoulli_internal (mpz_t *b, unsigned long n)
{
  if (n == 0)
    {
      mpz_init_set_ui (b[0], 1);
    }
  else
    {
      mpz_t t;
      unsigned long k;

      mpz_init (b[n]);
      /* b[n] = -sum(binomial(2n+1,2k)*C[k]*(2n)!/(2k+1)!,  k=0..n-1) */
      mpz_init_set_ui (t, 2 * n + 1);
      mpz_mul_ui (t, t, 2 * n - 1);
      mpz_mul_ui (t, t, 2 * n);
      mpz_mul_ui (t, t, n);
      mpz_fdiv_q_ui (t, t, 3); /* exact: t=binomial(2*n+1,2*k)*(2*n)!/(2*k+1)!
                               for k=n-1 */
      mpz_mul (b[n], t, b[n-1]);
      for (k = n - 1; k-- > 0;)
        {
          mpz_mul_ui (t, t, 2 * k + 1);
          mpz_mul_ui (t, t, 2 * k + 2);
          mpz_mul_ui (t, t, 2 * k + 2);
          mpz_mul_ui (t, t, 2 * k + 3);
          mpz_fdiv_q_ui (t, t, 2 * (n - k) + 1);
          mpz_fdiv_q_ui (t, t, 2 * (n - k));
          mpz_addmul (b[n], t, b[k]);
        }
      /* take into account C[1] */
      mpz_mul_ui (t, t, 2 * n + 1);
      mpz_fdiv_q_2exp (t, t, 1);
      mpz_sub (b[n], b[n], t);
      mpz_neg (b[n], b[n]);
      mpz_clear (t);
    }
  return;
}

static MPFR_THREAD_ATTR mpz_t *bernoulli_table = NULL;
static MPFR_THREAD_ATTR unsigned long bernoulli_size = 0;
static MPFR_THREAD_ATTR unsigned long bernoulli_alloc = 0;

mpz_srcptr
mpfr_bernoulli_cache (unsigned long n)
{
  unsigned long i;

  if (n >= bernoulli_size)
    {
      if (bernoulli_alloc == 0)
        {
          bernoulli_alloc = MAX(16, n + n/4);
          bernoulli_table = (mpz_t *)
            (*__gmp_allocate_func) (bernoulli_alloc * sizeof (mpz_t));
          bernoulli_size  = 0;
        }
      else if (n >= bernoulli_alloc)
        {
          bernoulli_table = (mpz_t *) (*__gmp_reallocate_func)
            (bernoulli_table, bernoulli_alloc * sizeof (mpz_t),
             (n + n/4) * sizeof (mpz_t));
          bernoulli_alloc = n + n/4;
        }
      MPFR_ASSERTD (bernoulli_alloc > n);
      MPFR_ASSERTD (bernoulli_size >= 0);
      for(i = bernoulli_size; i <= n; i++)
        mpfr_bernoulli_internal (bernoulli_table, i);
      bernoulli_size = n+1;
    }
  MPFR_ASSERTD (bernoulli_size > n);
  return bernoulli_table[n];
}

void
mpfr_bernoulli_freecache (void)
{
  unsigned long i;

  if (bernoulli_table != NULL)
    {
      for (i = 0; i < bernoulli_size; i++)
        {
          mpz_clear (bernoulli_table[i]);
        }
      (*__gmp_free_func) (bernoulli_table, bernoulli_alloc * sizeof (mpz_t));
      bernoulli_table = NULL;
      bernoulli_alloc = 0;
      bernoulli_size = 0;
    }
}
