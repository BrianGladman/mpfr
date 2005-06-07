/* tsum -- test file for the list summation function

Copyright 2004, 2005 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdlib.h>
#include <stdio.h>
#include "mpfr-test.h"


static int
check_is_sorted (unsigned long n, mpfr_srcptr *perm)
{
  unsigned long i;

  for (i = 0; i < n - 1; i++)
    if (MPFR_GET_EXP(perm[i]) < MPFR_GET_EXP(perm[i+1]))
      return 0;
  return 1;
}

static int
mpfr_list_sum (mpfr_ptr ret, mpfr_t *tab, unsigned long n, mp_rnd_t rnd)
{
    mpfr_ptr *tabtmp;
    unsigned long i;
    int inexact;
    MPFR_TMP_DECL(marker);
    
    MPFR_TMP_MARK(marker);
    tabtmp = (mpfr_ptr *) MPFR_TMP_ALLOC(n * sizeof(mpfr_srcptr));
    for (i = 0; i < n; i++)
        tabtmp[i] = tab[i];
    
    inexact = mpfr_sum (ret, tabtmp, n, rnd);
    MPFR_TMP_FREE(marker);
    return inexact;
}


static mp_prec_t
get_prec_max (mpfr_t *tab, unsigned long n, mp_prec_t f)
{
  mp_prec_t res;
  mp_exp_t min, max;
  unsigned long i;

  for (i = 0; MPFR_IS_ZERO (tab[i]); i++)
    MPFR_ASSERTD (i < n);
  min = max = MPFR_GET_EXP(tab[i]);
  for (i++; i < n; i++)
    {
      if (!MPFR_IS_ZERO (tab[i])) {
        if (MPFR_GET_EXP(tab[i]) > max)
          max = MPFR_GET_EXP(tab[i]);
        if (MPFR_GET_EXP(tab[i]) < min)
          min = MPFR_GET_EXP(tab[i]);
      }
    }
  res = max - min;
  res += f;
  res += __gmpfr_ceil_log2 (n) + 1;
  return res;
}


static void
algo_exact (mpfr_t somme, mpfr_t *tab, unsigned long n, mp_prec_t f)
{
  unsigned long i;
  mp_prec_t prec_max;

  prec_max = get_prec_max(tab, n, f);
  mpfr_set_prec (somme, prec_max);
  mpfr_set_ui (somme, 0, GMP_RNDN);
  for (i = 0; i < n; i++)
    {
      if (mpfr_add(somme, somme, tab[i], GMP_RNDN))
	{
          printf ("FIXME: algo_exact is buggy.\n");
          exit (1);
	}
    }
}

/* Test the sorting function */
static void
test_sort (mp_prec_t f, unsigned long n)
{
  mpfr_t *tab;
  mpfr_ptr *tabtmp;
  mpfr_srcptr *perm;
  unsigned long i;

  /* Init stuff */
  tab = (mpfr_t *) malloc (n * sizeof(mpfr_t));
  for (i = 0; i < n; i++)
    mpfr_init2 (tab[i], f);
  tabtmp = (mpfr_ptr *) malloc (n * sizeof(mpfr_ptr));
  perm = (mpfr_srcptr *) malloc (n * sizeof(mpfr_srcptr));

  for (i = 0; i < n; i++)
    {
      mpfr_random (tab[i]);
      tabtmp[i] = tab[i];
    }

  mpfr_count_sort (tabtmp, n, perm);

  if (check_is_sorted (n, perm) == 0)
    {
      printf ("mpfr_count_sort incorrect.\n");
      for (i = 0; i < n; i++)
        mpfr_dump (perm[i]);
      exit (1);
    }

  /* Clear stuff */
  for (i = 0; i < n; i++)
    mpfr_clear (tab[i]);
  free (tabtmp);
  free (perm);
  free (tab);
}

static void
test_sum (mp_prec_t f, unsigned long n)
{
  mpfr_t sum, real_sum, real_non_rounded;
  mpfr_t *tab;
  unsigned long i;
  int rnd_mode;

  /* Init */
  tab = (mpfr_t *) malloc (n * sizeof(mpfr_t));
  for (i = 0; i < n; i++)
    mpfr_init2 (tab[i], f);
  mpfr_inits2 (f, sum, real_sum, real_non_rounded, NULL);

  /* First Uniform */
  for (i = 0; i < n; i++)
    mpfr_random (tab[i]);
  algo_exact (real_non_rounded, tab, n, f);
  for (rnd_mode = 0; rnd_mode < GMP_RND_MAX; rnd_mode++)
    {
      mpfr_list_sum (sum, tab, n, (mp_rnd_t) rnd_mode);
      mpfr_set (real_sum, real_non_rounded, (mp_rnd_t) rnd_mode);
      if (mpfr_cmp (real_sum, sum) != 0)
        {
          printf ("mpfr_list_sum incorrect.\n");
          mpfr_dump (real_sum);
          mpfr_dump (sum);
          exit (1);
        }
    }

  /* Then non uniform */
  for (i = 0; i < n; i++)
    {
      mpfr_random (tab[i]);
      mpfr_set_exp (tab[i], randlimb () %1000);
    }
  algo_exact (real_non_rounded, tab, n, f);
  for (rnd_mode = 0; rnd_mode < GMP_RND_MAX; rnd_mode++)
    {
      mpfr_list_sum (sum, tab, n, (mp_rnd_t) rnd_mode);
      mpfr_set (real_sum, real_non_rounded, (mp_rnd_t) rnd_mode);
      if (mpfr_cmp (real_sum, sum) != 0)
        {
          printf ("mpfr_list_sum incorrect.\n");
          mpfr_dump (real_sum);
          mpfr_dump (sum);
          exit (1);
        }
    }

  /* Clear stuff */
  for (i = 0; i < n; i++)
    mpfr_clear (tab[i]);
  mpfr_clears (sum, real_sum, real_non_rounded, NULL);
  free (tab);
}

int
main (void)
{
  mp_prec_t p;
  unsigned long n;

  tests_start_mpfr ();

  test_sort (1764, 1026);
  for (p = 2 ; p < 1764 ; p+=17)
    for (n = 2 ; n < 1026 ; n+=42+p)
      test_sum (p, n);

  tests_end_mpfr ();
  return 0;
}
