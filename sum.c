/* Sum -- efficiently sum a list of floating-point numbers

Copyright 2004, 2005 Free Software Foundation, Inc.

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* Performs a counting sort of the entries */

static void heap_sort_exp_clean (mpfr_ptr const tab[], unsigned long n, 
                                 mpfr_srcptr *perm);

void mpfr_count_sort (mpfr_ptr const tab[], unsigned long n, 
                 mpfr_srcptr *perm)
{
    mp_exp_t min, max;
    unsigned long i;
    unsigned long *account;
    unsigned long exp_num;
    unsigned long target_rank;
    MPFR_TMP_DECL(marker);
    
    MPFR_TMP_MARK(marker);
    min = max = MPFR_GET_EXP(tab[0]);

    for (i = 1; i < n; i++)
    {
        if (MPFR_GET_EXP(tab[i]) < min)
            min = MPFR_GET_EXP(tab[i]);
        if (MPFR_GET_EXP(tab[i]) > max)
            max = MPFR_GET_EXP(tab[i]);
    }

    exp_num = max - min + 1;
    if (exp_num > (unsigned long) 42 * MPFR_INT_CEIL_LOG2 (n))
      /* FIXME : better test */
    {
        heap_sort_exp_clean (tab, n, perm);
        return;
    }
    account = (unsigned long *) MPFR_TMP_ALLOC(exp_num * sizeof(*account));
    for (i = 0; i < exp_num; i++)
        account[i] = 0;
    for (i = 0; i < n; i++)
        account[MPFR_GET_EXP(tab[i]) - min]++;
    for (i = exp_num - 1; i >= 1; i--)
        account[i - 1] += account[i];
    for (i = 0; i < n; i++)
    {
        target_rank = --account[MPFR_GET_EXP(tab[i]) - min];
        perm[target_rank] = tab[i];
    }
    
    MPFR_TMP_FREE(marker);
}

/* Performs a heap sort of the entries */

static void heap_sort_exp_clean (mpfr_ptr const tab[], unsigned long n, 
                                 mpfr_srcptr *perm)
{
  unsigned long dernier_traite;
  unsigned long i, pere;
  mpfr_srcptr tmp;
  unsigned long fils_gauche, fils_droit, fils_indigne;
  /* Reminder of a heap structure :
     node(i) has for left son node(2i +1) and right son node(2i)
     and father(node(i)) = node((i - 1) / 2)
  */
  
  /* initialize the permutation to identity */

  for (i = 0; i < n; i++)
    perm[i] = tab[i];

  /* insertion phase */

  for (dernier_traite = 1; dernier_traite < n; dernier_traite++)
    {
      i = dernier_traite;
      while (i > 0)
	{
	  pere = (i - 1) / 2;
	  if (MPFR_GET_EXP(perm[pere]) > MPFR_GET_EXP(perm[i]))
	    {
	      tmp = perm[pere];
	      perm[pere] = perm[i];
	      perm[i] = tmp;
	      i = pere;
	    }
	  else
	    break;
	}
    }

  /* extraction phase */
  
  for (dernier_traite = n - 1; dernier_traite > 0; dernier_traite--)
    {
      tmp = perm[0];
      perm[0] = perm[dernier_traite];
      perm[dernier_traite] = tmp;

      i = 0;
      while (1)
	{
	  fils_gauche = 2 * i + 1;
	  fils_droit = fils_gauche + 1;
	  if (fils_gauche < dernier_traite)
	    {
	      if (fils_droit < dernier_traite)
		{
		  if (MPFR_GET_EXP(perm[fils_droit]) < MPFR_GET_EXP(perm[fils_gauche]))
		    fils_indigne = fils_droit;
		  else
		    fils_indigne = fils_gauche;

		  if (MPFR_GET_EXP(perm[i]) > MPFR_GET_EXP(perm[fils_indigne]))
		    {
		      tmp = perm[i];
		      perm[i] = perm[fils_indigne];
		      perm[fils_indigne] = tmp;
		      i = fils_indigne;
		    }
		  else
		    break;
		}
	      else /* on a un fils gauche, pas de fils droit */
		{
		  if (MPFR_GET_EXP(perm[i]) > MPFR_GET_EXP(perm[fils_gauche]))
		    {
		      tmp = perm[i];
		      perm[i] = perm[fils_gauche];
		      perm[fils_gauche] = tmp;
		    }
		  break;
		}
	    }
	  else /* on n'a pas de fils */
	    break;
	}
    }
}


/* Sum a list of float with order given by permutation perm,
 * intermediate size set to F.
 * Internal use function.
 */

static int mpfr_list_sum_once (mpfr_ptr ret, mpfr_srcptr const tab[], 
                               unsigned long n, mp_prec_t F)
{
  unsigned long i;
  mpfr_t sum;
  int error_trap;

  if (MPFR_UNLIKELY (n == 1))
    return mpfr_set (ret, tab[0], GMP_RNDN);

  mpfr_init2 (sum, F);

  error_trap = mpfr_set (sum, tab[0], GMP_RNDN);
  for (i = 1; i < n - 1; i++)
    error_trap |= mpfr_add (sum, sum, tab[i], GMP_RNDN);
  error_trap |= mpfr_add (ret, sum, tab[n - 1], GMP_RNDN);

  mpfr_clear (sum);
  return error_trap;
}

/* Sum a list of floating-point numbers.
 * FIXME : add reference to Demmel-Hida's paper.
*/

int mpfr_sum (mpfr_ptr ret, mpfr_ptr const tab[], unsigned long n, 
              mp_rnd_t rnd)
{
  mp_prec_t prec;
  unsigned long k;
  mpfr_srcptr *perm;
  int error_trap;
  mpfr_t cur_sum;
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_TMP_DECL(marker);
    
  MPFR_TMP_MARK(marker);
  if (MPFR_UNLIKELY (n == 0)) {
    MPFR_SET_ZERO (ret);
    MPFR_SET_POS (ret);
    return 0;
  }

  /* Sort */
  perm = (mpfr_srcptr *) MPFR_TMP_ALLOC(n * sizeof(mpfr_srcptr)); 
  mpfr_count_sort (tab, n, perm);

  /* Initial precision */
  prec = MAX (MPFR_PREC (tab[0]), MPFR_PREC (ret));
  k = MPFR_INT_CEIL_LOG2 (n) + 1;
  prec += k + 2;
  mpfr_init2 (cur_sum, prec);

  MPFR_SAVE_EXPO_MARK (expo);

  /* Ziv Loop */
  MPFR_ZIV_INIT (loop, prec);
  for (;;) {
    error_trap = mpfr_list_sum_once (cur_sum, perm, n, prec + k);
    if (MPFR_LIKELY (error_trap == 0
		     || mpfr_can_round (cur_sum, 
					MPFR_GET_EXP (cur_sum) - prec + 2,
					GMP_RNDN, rnd, MPFR_PREC (ret))))
      break;
    MPFR_ZIV_NEXT (loop, prec);
    mpfr_set_prec (cur_sum, prec);
  }
  MPFR_ZIV_FREE (loop);

  error_trap |= mpfr_set (ret, cur_sum, rnd);
  mpfr_clear (cur_sum);
  MPFR_TMP_FREE(marker);

  MPFR_SAVE_EXPO_FREE (expo);
  error_trap |= mpfr_check_range (ret, 0, rnd);
  return error_trap; /* It doesn't return the ternary value */
}


/* __END__ */
