/* mpfr_free_cache - Free the cache used by MPFR for internal consts.

Copyright 2004-2017 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

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

/* Default value for the cache of mpz_t */
#ifndef MPFR_MY_MPZ_INIT
#  define MPFR_MY_MPZ_INIT 32
#endif

/* If the number of value to cache is not zero */
#if MPFR_MY_MPZ_INIT

/* Index in the stack table of mpz_t and stack table of mpz_t */
static MPFR_THREAD_ATTR int n_alloc = 0;
static MPFR_THREAD_ATTR __mpz_struct mpz_tab[MPFR_MY_MPZ_INIT];

MPFR_HOT_FUNCTION_ATTR void
mpfr_mpz_init (mpz_t z)
{
  if (MPFR_LIKELY (n_alloc > 0))
    {
      /* Get a mpz_t from the MPFR stack of previously used mpz_t.
         It reduces memory pressure, and it allows to reuse
         a mpz_t which should be sufficiently big. */
      MPFR_ASSERTD (n_alloc <= numberof (mpz_tab));
      memcpy (z, &mpz_tab[--n_alloc], sizeof (mpz_t));
      SIZ(z) = 0;
    }
  else
    {
      /* Call real GMP function */
      (__gmpz_init)(z);
    }
}

MPFR_HOT_FUNCTION_ATTR void
mpfr_mpz_init2 (mpz_t z, mp_bitcnt_t n)
{
  if (MPFR_LIKELY (n_alloc > 0))
    {
      /* Get a mpz_t from the MPFR stack of previously used mpz_t.
         It reduces memory pressure, and it allows to reuse
         a mpz_t which should be sufficiently big. */
      MPFR_ASSERTD (n_alloc <= numberof (mpz_tab));
      memcpy (z, &mpz_tab[--n_alloc], sizeof (mpz_t));
      SIZ(z) = 0;
    }
  else
    {
      /* Call real GMP function */
      (__gmpz_init2)(z, n);
    }
}

MPFR_HOT_FUNCTION_ATTR void
mpfr_mpz_clear (mpz_t z)
{
  if (MPFR_LIKELY (n_alloc < numberof (mpz_tab)))
    {
      /* Push back the mpz_t inside the stack of the used mpz_t */
      MPFR_ASSERTD (n_alloc >= 0);
      memcpy (&mpz_tab[n_alloc++], z, sizeof (mpz_t));
    }
  else
    {
      /* Call real GMP function */
      (__gmpz_clear)(z);
    }
}

#endif

/* Theses caches may be global to all threads or local to the current */
static void
mpfr_free_const_caches (void)
{
#ifndef MPFR_USE_LOGGING
  mpfr_clear_cache (__gmpfr_cache_const_pi);
  mpfr_clear_cache (__gmpfr_cache_const_log2);
#else
  mpfr_clear_cache (__gmpfr_normal_pi);
  mpfr_clear_cache (__gmpfr_normal_log2);
  mpfr_clear_cache (__gmpfr_logging_pi);
  mpfr_clear_cache (__gmpfr_logging_log2);
#endif
  mpfr_clear_cache (__gmpfr_cache_const_euler);
  mpfr_clear_cache (__gmpfr_cache_const_catalan);
}

/* Theses caches are always local to a thread */
static void
mpfr_free_local_cache (void)
{
  /* Before mpz caching */
  mpfr_bernoulli_freecache();

#if MPFR_MY_MPZ_INIT
  { /* Avoid mixed declarations and code for ISO C90 support. */
    int i;
    MPFR_ASSERTD (n_alloc >= 0 && n_alloc <= numberof (mpz_tab));
    for (i = 0; i < n_alloc; i++)
      (__gmpz_clear)(&mpz_tab[i]);
    n_alloc = 0;
  }
#endif
}

void
mpfr_free_cache (void)
{
  mpfr_free_local_cache();
  mpfr_free_const_caches ();
}

void
mpfr_free_cache2 (mpfr_free_cache_t way)
{
  if (way & MPFR_FREE_LOCAL_CACHE)
    {
      mpfr_free_local_cache();
#if !defined (WANT_SHARED_CACHE)
      mpfr_free_const_caches ();
#endif
    }
  if (way & MPFR_FREE_GLOBAL_CACHE)
    {
#if defined (WANT_SHARED_CACHE)
      mpfr_free_const_caches ();
#endif
    }
}
