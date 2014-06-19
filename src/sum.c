/* Sum -- efficiently sum a list of floating-point numbers

Copyright 2004-2014 Free Software Foundation, Inc.
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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* FIXME[VL].

Preliminary note: The previous mpfr_sum algorithm was a high-level one
and could take too much memory and too much time in cases with numbers
of different magnitudes and cancellation. Here we will use functions
from the mpn layer of GMP only (no mpfr_add), as the algorithm is based
on additions by blocks, involving virtual splitting of MPFR numbers.

The general ideas of the algorithm:

0. Handle the case n <= 1 separately (+0 for n = 0, mpfr_set for n = 1).
   Now assume that n >= 2.

1. Look at the exponent field of each mpfr_t input (this is fast); also
   track the signs of Inf's and zeros. This step allows one to:
   * detect the singular cases, i.e. when there is either at least a NaN
     or an Inf in the inputs (in case of a NaN or two Inf's of different
     signs, we can immediately return NaN), or all inputs are zeros;
   * determine the maximum exponent maxexp of the regular inputs and
     the number n' of regular inputs.

2. Compute the truncated sum in two's complement by taking into account
   the part of the inputs in a window [minexp,maxexp[ (we choose not to
   include maxexp in order to match the floating-point representation
   chosen in MPFR, where the significand is in [1/2,1[; this also means
   that minexp will be the next maxexp, unless there is a hole in the
   inputs, as explained below).
   Due to the accumulation of values and the choice of two's complement
   (a way to represent the sign), we will need some extra bits to avoid
   overflows. The absolute value of the sum is less than n' * 2^maxexp,
   taken up to ceil(log2(n')) extra bits, and one needs one more bit to
   be able to determine the sign, so that cq = ceil(log2(n')) + 1 extra
   bits will be considered.
   Moreover we will choose minexp = maxexp - output_prec - dq, where dq
   will be chosen around log2(n') to take into account the accumulation
   of errors, i.e. from everything less significant than minexp. The
   final choice for dq should be done after testing: in practice, one
   may need a little more precision to handle partial cancellation at
   this iteration, but important cancellation will always lead to other
   iterations. For instance, we will choose the smallest dq > log2(n')
   such that the window bitsize, which is maxexp + cq - minexp, i.e.
   cq + output_prec + dq, is a multiple of GMP_NUMB_BITS.
   In the same loop over the inputs, determine the maximum exponent
   maxexp2 of the numbers formed starting with the most significant
   represented bit that has been ignored: one will get either minexp
   (if an input has been truncated at this iteration) or the maximum
   exponent of the numbers that have been completely ignored.
   Note: for simplicity of the maxexp2 computation, any value > minexp
   will be regarded as equivalent to minexp during the iteration.

3. If applicable (see below), add both windows.

4. Determine the number of cancelled bits.

5. If the truncated sum is 0, reiterate at (2) with maxexp = maxexp2,
   rounded above to a multiple of GMP_NUMB_BITS (see below).

6. If the error is too large, shift the truncated sum to the left of
   the window (most significant part), and reiterate at (2) with a
   second window (least significant part). Using a second window is
   useful to avoid carry propagation (potentially for each add/sub)
   to the full window: it is expected that in general, the second
   window will be small (small cancellation). The cumulated size of
   both windows should be no more than: output_prec + k * log2(n'),
   where k is a small constant.

7. If only the sign of the error term is unknown, reiterate at (2)
   to compute it using a second window where output_prec = 0, i.e.
   around maxexp + cq to maxexp - dq.
   Note: the sign of the error term is needed to round the result in
   the right direction and/or to determine the ternary value.

8. Copy the rounded result to the destination, and exit with the
   correct ternary value.

An example:

              [1]    [2]   A  [3]
  u0 = *****   |      |    .   |
  u1 =   ******|**    |    .   |
  u2 = ********|***** |    .   |
  u3 =         | *****|    .   |
  u4 =         |      |    ****|**
  u5 =         |      |     ***|*

At iteration 1, minexp is determined to be [1]; thus u0, a part of u1,
and a part of u2 are taken into account for the truncated sum. Then it
appears that an important cancellation occurred, and another step (2)
is needed. Since u1 was truncated, the new maxexp will be minexp, i.e.
[1]. At iteration 2, minexp is determined to be [2]; thus a part of u1,
a part of u2, and u3 are taken into account for the truncated sum. Now
assume that on this example, the error is small enough, but its sign is
unknown. Thus another step (2) is needed, with the conditions of (7).
Since no numbers were truncated at the previous iteration, maxexp is
the maximum exponent of the remaining numbers, here the one of u4, and
minexp is determined to be [3]. Assume that the sign of the error can
be determined now, so that we can return the rounded result with the
ternary value.

As a bonus, this will also solve overflow, underflow and normalization
issues, since everything is done in fixed point and the output exponent
will be considered only at the end (early overflow detection could also
be done).

Both windows will occupy the same area, but the limit between them
(thus their sizes) will possibly change at some iterations. At worst,
the size of the window for (2) will be around output_prec + 2*log2(n')
at the first iterations (as long as a full cancellation is possible),
then will become something intermediate at some iteration, then will
be around 2*log2(n') at the next iterations.

A limb L = *zp in memory will generally contain a part of a significand.
One can define its exponent ze, such that the actual value of this limb
is L * 2^(ze-GMP_NUMB_BITS), i.e. ze is its exponent where the limb is
regarded as a number in [1/2,1[. If an array of limbs zp[] is regarded
as a significand in [1/2,1[, then the exponent of its actual value is
also ze.

Variables:
  tp: pointer to a temporary area that will contain a shifted value.
  xp: pointer to the main window.
  yp: pointer to the secondary window (yp = xp at iteration 1).
  ts: the maximum size of the temporary area, in limbs; then xp = tp+ts.
  ws: the maximum size of both main & secondary windows, in limbs.
  xs: number of limbs of the main window already determined;
      the secondary window will point at xp+xs.
  ys: number of limbs of the secondary window; check that xs + ys <= ws.
  rn: number n' of inputs that are regular numbers (regular inputs).
  logn: ceil(log2(rn)).
  cq: logn + 1.
  dq: > logn, such that the window bitsize is a multiple of GMP_NUMB_BITS.
  sq: output precision (precision of the sum).
  maxexp: the maximum exponent of the bit-window of the inputs that is
          taken into account (for the current iteration), excluded.
  minexp: the minimum exponent of the bit-window of the inputs that is
          taken into account (for the current iteration), included.

*** TODO ***
It seems that carry propagation (mpn_add_1 & mpn_sub_1 in the code) is
most often limited. But consider the following cases, where all inputs
have the minimal precision 2, and the output precision is p:
  u0 = 1
  u_i = (-1)^i * 2^(-p) for i > 0
Here long carry propagation will occur for each addition of the initial
iteration, so that the complexity will be O(n*p) instead of O(n+p) if
we choose to delay carry propagation (however such a choice may slower
the average case and take more memory, such as around 3*p instead of
2*p). Using a second window when a new iteration is needed can help in
some cases, such as:
  u0 = 2^q with some q > 0
  u1 = 1
  u2 = -2^q
  u_i = (-1)^i * 2^(-p) for i > 2
but such examples are very specific cases, and as seen with the first
example, a better way (e.g. delaying carry propagation) must be found
if one wants to avoid long carry propagation in all cases.

--------

Note: see the following paper and its references:
http://www.eecs.berkeley.edu/~hdnguyen/public/papers/ARITH21_Fast_Sum.pdf
VL: This is very different:
          In MPFR             In the paper & references
    arbitrary precision            fixed precision
     correct rounding        just reproducible rounding
    integer operations        floating-point operations
        sequencial             parallel (& sequential)
*/

/* I would really like to use "mpfr_srcptr const []" but the norm is buggy:
   it doesn't automatically cast a "mpfr_ptr []" to "mpfr_srcptr const []"
   if necessary. So the choice are:
     mpfr_s **                : OK
     mpfr_s *const*           : OK
     mpfr_s **const           : OK
     mpfr_s *const*const      : OK
     const mpfr_s *const*     : no
     const mpfr_s **const     : no
     const mpfr_s *const*const: no
   VL: this is not a bug, but a feature. See the reason here:
     http://c-faq.com/ansi/constmismatch.html
*/
static void heap_sort (mpfr_srcptr *const, unsigned long, mpfr_srcptr *);
static void count_sort (mpfr_srcptr *const, unsigned long, mpfr_srcptr *,
                        mpfr_exp_t, mpfr_uexp_t);

/* Either sort the tab in perm and returns 0
   Or returns 1 for +INF, -1 for -INF and 2 for NAN.
   Also set *maxprec to the maximal precision of tab[0..n-1] and of the
   initial value of *maxprec.
*/
int
mpfr_sum_sort (mpfr_srcptr *const tab, unsigned long n, mpfr_srcptr *perm,
               mpfr_prec_t *maxprec)
{
  mpfr_exp_t min, max;
  mpfr_uexp_t exp_num;
  unsigned long i;
  int sign_inf;

  sign_inf = 0;
  min = MPFR_EMIN_MAX;
  max = MPFR_EMAX_MIN;
  for (i = 0; i < n; i++)
    {
      if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (tab[i])))
        {
          if (MPFR_IS_NAN (tab[i]))
            return 2; /* Return NAN code */
          else if (MPFR_IS_INF (tab[i]))
            {
              if (sign_inf == 0) /* No previous INF */
                sign_inf = MPFR_SIGN (tab[i]);
              else if (sign_inf != MPFR_SIGN (tab[i]))
                return 2; /* Return NAN */
            }
        }
      else
        {
          MPFR_ASSERTD (MPFR_IS_PURE_FP (tab[i]));
          if (MPFR_GET_EXP (tab[i]) < min)
            min = MPFR_GET_EXP(tab[i]);
          if (MPFR_GET_EXP (tab[i]) > max)
            max = MPFR_GET_EXP(tab[i]);
        }
      if (MPFR_PREC (tab[i]) > *maxprec)
        *maxprec = MPFR_PREC (tab[i]);
    }
  if (MPFR_UNLIKELY (sign_inf != 0))
    return sign_inf;

  exp_num = max - min + 1;
  /* FIXME : better test */
  if (exp_num > n * MPFR_INT_CEIL_LOG2 (n))
    heap_sort (tab, n, perm);
  else
    count_sort (tab, n, perm, min, exp_num);
  return 0;
}

#define GET_EXP1(x) (MPFR_IS_ZERO (x) ? min : MPFR_GET_EXP (x))
/* Performs a count sort of the entries */
static void
count_sort (mpfr_srcptr *const tab, unsigned long n,
            mpfr_srcptr *perm, mpfr_exp_t min, mpfr_uexp_t exp_num)
{
  unsigned long *account;
  unsigned long target_rank, i;
  MPFR_TMP_DECL(marker);

  /* Reserve a place for potential 0 (with EXP min-1)
     If there is no zero, we only lose one unused entry */
  min--;
  exp_num++;

  /* Performs a counting sort of the entries */
  MPFR_TMP_MARK (marker);
  account = (unsigned long *) MPFR_TMP_ALLOC (exp_num * sizeof *account);
  for (i = 0; i < exp_num; i++)
    account[i] = 0;
  for (i = 0; i < n; i++)
    account[GET_EXP1 (tab[i]) - min]++;
  for (i = exp_num - 1; i >= 1; i--)
    account[i - 1] += account[i];
  for (i = 0; i < n; i++)
    {
      target_rank = --account[GET_EXP1 (tab[i]) - min];
      perm[target_rank] = tab[i];
    }
  MPFR_TMP_FREE (marker);
}


#define GET_EXP2(x) (MPFR_IS_ZERO (x) ? MPFR_EMIN_MIN : MPFR_GET_EXP (x))

/* Performs a heap sort of the entries */
static void
heap_sort (mpfr_srcptr *const tab, unsigned long n, mpfr_srcptr *perm)
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
          if (GET_EXP2 (perm[pere]) > GET_EXP2 (perm[i]))
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
                  if (GET_EXP2(perm[fils_droit]) < GET_EXP2(perm[fils_gauche]))
                    fils_indigne = fils_droit;
                  else
                    fils_indigne = fils_gauche;

                  if (GET_EXP2 (perm[i]) > GET_EXP2 (perm[fils_indigne]))
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
                  if (GET_EXP2 (perm[i]) > GET_EXP2 (perm[fils_gauche]))
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
 * intermediate size set to F. Return non-zero if at least one of
 * the operations is inexact (thus 0 implies that the sum is exact).
 * Internal use function.
 */
static int
sum_once (mpfr_ptr ret, mpfr_srcptr *const tab, unsigned long n, mpfr_prec_t F)
{
  mpfr_t sum;
  unsigned long i;
  int error_trap;

  MPFR_ASSERTD (n >= 2);

  mpfr_init2 (sum, F);
  error_trap = mpfr_set (sum, tab[0], MPFR_RNDN);
  for (i = 1; i < n - 1; i++)
    {
      MPFR_ASSERTD (!MPFR_IS_NAN (sum) && !MPFR_IS_INF (sum));
      if (mpfr_add (sum, sum, tab[i], MPFR_RNDN))
        error_trap = 1;
    }
  if (mpfr_add (ret, sum, tab[n - 1], MPFR_RNDN))
    error_trap = 1;
  mpfr_clear (sum);
  return error_trap;
}

/* The following function will disappear in the final code. */
static int
mpfr_sum_old (mpfr_ptr ret, mpfr_ptr *const tab_p, unsigned long n, mpfr_rnd_t rnd)
{
  mpfr_t cur_sum;
  mpfr_prec_t prec;
  mpfr_srcptr *perm, *const tab = (mpfr_srcptr *) tab_p;
  int k, error_trap;
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_TMP_DECL (marker);

  if (MPFR_UNLIKELY (n <= 1))
    {
      if (n < 1)
        {
          MPFR_SET_ZERO (ret);
          MPFR_SET_POS (ret);
          return 0;
        }
      else
        return mpfr_set (ret, tab[0], rnd);
    }

  /* Sort and treat special cases */
  MPFR_TMP_MARK (marker);
  perm = (mpfr_srcptr *) MPFR_TMP_ALLOC (n * sizeof *perm);
  prec = MPFR_PREC (ret);
  error_trap = mpfr_sum_sort (tab, n, perm, &prec);
  /* Check if there was a NAN or a INF */
  if (MPFR_UNLIKELY (error_trap != 0))
    {
      MPFR_TMP_FREE (marker);
      if (error_trap == 2)
        {
          MPFR_SET_NAN (ret);
          MPFR_RET_NAN;
        }
      MPFR_SET_INF (ret);
      MPFR_SET_SIGN (ret, error_trap);
      MPFR_RET (0);
    }

  /* Initial precision is max(prec(ret),prec(tab[0]),...,prec(tab[n-1])) */
  k = MPFR_INT_CEIL_LOG2 (n) + 1;
  prec += k + 2;
  mpfr_init2 (cur_sum, prec);

  /* Ziv Loop */
  MPFR_SAVE_EXPO_MARK (expo);
  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      error_trap = sum_once (cur_sum, perm, n, prec + k);
      if (MPFR_LIKELY (error_trap == 0 ||
                       (!MPFR_IS_ZERO (cur_sum) &&
                        mpfr_can_round (cur_sum, prec - 2,
                                        MPFR_RNDN, rnd, MPFR_PREC (ret)))))
        break;
      MPFR_ZIV_NEXT (loop, prec);
      mpfr_set_prec (cur_sum, prec);
    }
  MPFR_ZIV_FREE (loop);
  MPFR_TMP_FREE (marker);

  if (mpfr_set (ret, cur_sum, rnd))
    error_trap = 1;
  mpfr_clear (cur_sum);

  MPFR_SAVE_EXPO_FREE (expo);
  if (mpfr_check_range (ret, 0, rnd))
    error_trap = 1;
  return error_trap; /* It doesn't return the ternary value */
}

int
mpfr_sum (mpfr_ptr sum, mpfr_ptr *const p, unsigned long n, mpfr_rnd_t rnd)
{
  MPFR_LOG_FUNC
    (("n=%lu rnd=%d", n, rnd),
     ("sum[%Pu]=%.*Rg", mpfr_get_prec (sum), mpfr_log_prec, sum));

  if (MPFR_UNLIKELY (n <= 1))
    {
      if (n < 1)
        {
          MPFR_SET_ZERO (sum);
          MPFR_SET_POS (sum);
          MPFR_RET (0);
        }
      else
        return mpfr_set (sum, p[0], rnd);
    }
  else
    {
      mp_limb_t *tp;  /* pointer to a temporary area (fixed) */
      mp_limb_t *xp;  /* pointer to the main window (fixed) */
      mp_size_t ts;   /* maximum size of the temporary area, in limbs */
      mp_size_t ws;   /* maximum size of both main & secondary windows */
      mp_size_t xs;   /* #limbs of the main window already determined */
      mp_size_t ys;   /* number of limbs of the secondary window */
      mpfr_exp_t maxexp;
      unsigned long i, rn;
      int logn;       /* ceil(log2(rn)) */
      int cq, dq;
      mpfr_prec_t sq;
      MPFR_TMP_DECL (marker);

      /* Pre-iteration (Step 1) */
      {
        /* sign of infinities and zeros (0: currently unknown) */
        int sign_inf = 0, sign_zero = 0;

        rn = 0;  /* will be the number of regular inputs */
        maxexp = MPFR_EXP_MIN;  /* max(Empty), <= any valid exponent */
        for (i = 0; i < n; i++)
          {
            if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (p[i])))
              {
                if (MPFR_IS_NAN (p[i]))
                  {
                    /* The current value p[i] is NaN. Then the sum is NaN. */
                  nan:
                    MPFR_SET_NAN (sum);
                    MPFR_RET_NAN;
                  }
                else if (MPFR_IS_INF (p[i]))
                  {
                    /* The current value p[i] is an infinity.
                       There are two cases:
                       1. This is the first infinity value (sign_inf == 0).
                          Then set sign_inf to its sign, and go on.
                       2. All the infinities found until now have the same
                          sign sign_inf. If this new infinity has a different
                          sign, then return NaN immediately, else go on. */
                    if (sign_inf == 0)
                      sign_inf = MPFR_SIGN (p[i]);
                    else if (MPFR_SIGN (p[i]) != sign_inf)
                      goto nan;
                  }
                else if (MPFR_UNLIKELY (rn == 0))
                  {
                    /* The current value p[i] is a zero. The code below
                       matters only when all values found until now are
                       zeros, otherwise it is harmless (the test rn == 0
                       above is just a minor optimization).
                       Here we track the sign of the zero result when all
                       inputs are zeros: if all zeros have the same sign,
                       the result will have this sign, otherwise (i.e. if
                       there is at least a zero of each sign), the sign of
                       the zero result depends only on the rounding mode
                       (note that this choice is sticky when new zeros are
                       considered). */
                    MPFR_ASSERTD (MPFR_IS_ZERO (p[i]));
                    if (sign_zero == 0)
                      sign_zero = MPFR_SIGN (p[i]);
                    else if (MPFR_SIGN (p[i]) != sign_zero)
                      sign_zero = rnd == MPFR_RNDD ? -1 : 1;
                  }
              }
            else
              {
                /* The current value p[i] is a regular number. */
                mpfr_exp_t e = MPFR_GET_EXP (p[i]);
                if (e > maxexp)
                  maxexp = e;  /* maximum exponent found until now */
                rn++;  /* current number of regular inputs */
              }
          }

        MPFR_LOG_MSG (("rn=%lu sign_inf=%d sign_zero=%d\n",
                       rn, sign_inf, sign_zero));

        /* At this point the result cannot be NaN (this case has already
           been filtered out). */

        if (MPFR_UNLIKELY (sign_inf != 0))
          {
            /* At least one infinity, and all of them have the same sign
               sign_inf. The sum is the infinity of this sign. */
            MPFR_SET_INF (sum);
            MPFR_SET_SIGN (sum, sign_inf);
            MPFR_RET (0);
          }

        /* At this point, all the inputs are finite numbers. */

        if (MPFR_UNLIKELY (rn == 0))
          {
            /* All the numbers were zeros (and there is at least one).
               The sum is zero with sign sign_zero. */
            MPFR_ASSERTD (sign_zero != 0);
            MPFR_SET_ZERO (sum);
            MPFR_SET_SIGN (sum, sign_zero);
            MPFR_RET (0);
          }

      } /* End of the pre-iteration (Step 1) */

      /* Generic case: all the inputs are finite numbers, with at least
         a regular number. */

      /* rn is the number of regular inputs (the singular ones will be
         ignored). Compute logn = ceil(log2(rn)). */
      logn = MPFR_INT_CEIL_LOG2 (rn);

      MPFR_LOG_MSG (("logn=%d maxexp=%" MPFR_EXP_FSPEC "d\n",
                     logn, (mpfr_eexp_t) maxexp));

      sq = MPFR_GET_PREC (sum);
      cq = logn + 1;

      /* First determine the maximum size of the main window; it will also
         be the size of the secondary window at the first iteration. */
      ys = MPFR_PREC2LIMBS (cq + sq + logn + 1);

      /* We now choose dq so that cq + sq + dq = maximum bitsize of the
         main window. */
      dq = (mpfr_prec_t) ys * GMP_NUMB_BITS - (cq + sq);

      /* Take the secondary window into account. Mainly for Step 7. */
      /* TODO: check this in the final code... */
      ws = ys + MPFR_PREC2LIMBS (cq + dq);

      /* An input block will have up to sq + dq bits, and its shifted
         value (to be correctly aligned) may take GMP_NUMB_BITS - 1
         additional bits. */
      ts = MPFR_PREC2LIMBS (sq + dq + GMP_NUMB_BITS - 1);

      MPFR_TMP_MARK (marker);

      tp = MPFR_TMP_LIMBS_ALLOC (ts + ws);
      xp = tp + ts;
      xs = 0;   /* initial size of the main window: 0 */

      do
        {
          mp_limb_t *yp = xp + xs;  /* pointer to the secondary window */
          mpfr_exp_t minexp = maxexp - (sq + dq);
          mpfr_exp_t maxexp2 = MPFR_EXP_MIN;  /* < any valid exponent */

          MPFR_ASSERTD (xs + ys <= ws);

          /* TODO: this may not be necessary if the condition of Step 5
             is satisfied. */
          MPN_ZERO (yp, ys);  /* zero the secondary window */

          for (i = 0; i < n; i++)
            if (! MPFR_IS_SINGULAR (p[i]))
              {
                mp_limb_t *vp;
                mp_size_t vs;
                mpfr_exp_t pexp, vd;

                pexp = MPFR_GET_EXP (p[i]);

                vp = MPFR_MANT (p[i]);
                vs = MPFR_PREC2LIMBS (MPFR_GET_PREC (p[i]));
                vd = pexp - vs * GMP_NUMB_BITS - minexp;
                /* vd is the exponent of the least significant represented
                   bit of p[i] (including the trailing bits at 0) minus the
                   exponent of the least significant bit of the window. */

                if (vd < 0)
                  {
                    mp_size_t vds;

                    /* At least one bit of p[i] will be ignored. Let's
                       update maxexp2. Note: any value > minexp will be
                       changed to minexp after the iteration if maxexp2
                       needs to be taken into account. */
                    if (pexp > maxexp2)
                      maxexp2 = pexp;

                    /* We need to ignore the least |vd| significant bits
                       of p[i]. First, let's ignore the least
                       vds = |vd| / GMP_NUMB_BITS limbs. */
                    vd = - vd;
                    vds = vd / GMP_NUMB_BITS;
                    vs -= vds;
                    if (vs <= 0)  /* No overlapping between p[i] */
                      continue;   /* and the window.             */
                    vp += vds;
                    vd -= vds * GMP_NUMB_BITS;
                    MPFR_ASSERTD (vd >= 0 && vd < GMP_NUMB_BITS);

                    if (vd != 0)
                      {
                        mpn_rshift (tp, vp, vs > ys ? (vs = ys, ys + 1) : vs,
                                    vd);
                        vp = tp;
                      }
                    else if (vs > ys)
                      vs = ys;
                    MPFR_ASSERTD (vs <= ys);

                    if (MPFR_IS_POS (p[i]))
                      {
                        mp_limb_t carry;

                        carry = mpn_add_n (yp, yp, vp, vs);
                        if (carry != 0 && vs < ys)
                          mpn_add_1 (yp + vs, yp + vs, ys - vs, carry);
                      }
                    else
                      {
                        mp_limb_t borrow;

                        borrow = mpn_sub_n (yp, yp, vp, vs);
                        if (borrow != 0 && vs < ys)
                          mpn_sub_1 (yp + vs, yp + vs, ys - vs, borrow);
                      }
                  }
                else
                  {
                    mp_limb_t *dp;
                    mp_size_t ds, vds;

                    /* We need to ignore the least vd significant bits
                       of the window. First, let's ignore the least
                       vds = vd / GMP_NUMB_BITS limbs. -> (dp,ds) */
                    vds = vd / GMP_NUMB_BITS;
                    ds = ys - vds;
                    if (ds <= 0)
                      continue;
                    dp = yp + vds;
                    vd -= vds * GMP_NUMB_BITS;
                    MPFR_ASSERTD (vd >= 0 && vd < GMP_NUMB_BITS);

                    if (vs > ds)
                      vs = ds;
                    if (vd != 0)
                      {
                        mpn_lshift (tp, vp, vs, vd);
                        vp = tp;
                      }

                    if (MPFR_IS_POS (p[i]))
                      {
                        mp_limb_t carry;

                        carry = mpn_add_n (dp, dp, vp, vs);
                        if (carry != 0 && vs < ds)
                          mpn_add_1 (dp + vs, dp + vs, ds - vs, carry);
                      }
                    else
                      {
                        mp_limb_t borrow;

                        borrow = mpn_sub_n (dp, dp, vp, vs);
                        if (borrow != 0 && vs < ds)
                          mpn_sub_1 (dp + vs, dp + vs, ds - vs, borrow);
                      }
                  }
              }

          /* TODO: implement steps 3 to 8. */





          /* Obsolete code, to be updated... */

          /* The truncated sum has been computed. Let's determine its
             sign, absolute value, and the number of significant bits
             (starting from the most significant 1).
             Note: the mpn_neg may be useless, but it doesn't waste
             much time, it's simpler and it makes the code shorter
             than analyzing different cases. */
          wi = ys - 1;
          neg = MPFR_LIMB_MSB (yp[wi]) != 0;
          /* FIXME: Do not do the mpn_neg inside the loop since in
             case of another iteration, we want to keep the number
             in the same representation, but shifted to the left of
             the window. */
          if (neg)
            mpn_neg (yp, yp, ys);
          signif = wq;
          for (; wi >= 0; wi--)
            {
              if (yp[wi] != 0)
                {
                  int cnt;
                  count_leading_zeros (cnt, yp[wi]);
                  signif -= cnt;
                  break;
                }
              signif -= GMP_NUMB_BITS;
            }
          MPFR_LOG_MSG (("neg=%d signif=%Pu\n", neg, signif));

          if (MPFR_UNLIKELY (signif == 0))
            {
              /* There may be a big hole between numbers! We need to
                 determine a new maximum exponent among the numbers
                 with at least a represented bit of exponent <= minexp.
                 Note that there may be no such numbers, in which case
                 the exact sum is 0. */

            }

        }
      while (0); /* the condition will be determined later */

      MPFR_TMP_FREE (marker);

      /* ... */

      return mpfr_sum_old (sum, p, n, rnd);
    }
}
