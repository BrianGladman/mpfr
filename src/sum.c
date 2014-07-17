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

/*
Preliminary note: The previous mpfr_sum algorithm was a high-level one
and could take too much memory and too much time in cases with numbers
of different magnitudes and cancellation. Here we will use functions
from the mpn layer of GMP only (no mpfr_add), as the algorithm is based
on additions by blocks, involving virtual splitting of MPFR numbers.

Let sq be the output precision.

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

2. A truncated sum will be computed in two's complement, stored in a
   fixed memory area, called "accumulator", whose characteristics are
   determined here.

   Only a fixed part of the inputs will be taken into account for the
   truncated sum: the bits whose exponent is in some window (interval)
   denoted [minexp,maxexp[.

   We choose not to include maxexp in the interval in order to match the
   floating-point representation chosen in MPFR, where the significand
   is in [1/2,1[; this also means that the current minexp will be maxexp
   at the next iteration, unless there is a "hole" between the inputs,
   as explained below.

   Due to the accumulation of values and the choice of two's complement
   (a way to represent the sign), we will need some extra bits to avoid
   overflows. The absolute value of the sum is less than n' * 2^maxexp,
   taken up to ceil(log2(n')) extra bits, and one needs one more bit to
   be able to determine the sign, so that cq = ceil(log2(n')) + 1 extra
   bits will be considered.

   For the other side, we define minexp = maxexp - sq - dq, where dq
   will be chosen around log2(n') to take into account the accumulation
   of errors, i.e. from everything less significant than minexp. The
   final choice for dq should be done after testing: in practice, one
   may need a little more precision to handle partial cancellation at
   this iteration, but important cancellation will always lead to other
   iterations. For instance, we will choose the smallest dq > log2(n')
   such that the window bitsize, which is maxexp + cq - minexp, i.e.
   cq + sq + dq, is a multiple of GMP_NUMB_BITS.

3. Compute the truncated sum in two's complement by taking into account
   the part of the inputs in the window [minexp,maxexp[.

   In the same loop over the inputs, determine the maximum exponent
   maxexp2 of the numbers formed starting with the most significant
   represented bit that has been ignored: one will get either minexp
   (if an input has been truncated at this iteration) or the maximum
   exponent of the numbers that have been completely ignored.

4. Determine the number of cancelled bits.

5. If the truncated sum (i.e. the value of the accumulator) is 0, then
   reiterate at (3) with maxexp = maxexp2.

6. If, because of cancellation, the correctly rounded sum cannot be
   determined possibly by using the sign of the error term, then shift
   the truncated sum to the left boundary (most significant part) of
   the accumulator, and reiterate at (3), where:
   * the shift count is determined in such a way to avoid overflows
     at the next iteration, i.e. to be able to retrieve the sum with
     its sign;
   * the new value of cq is implied by this shift count and maxexp2
     (the value of maxexp at the next iteration).
   Note: It is expected that in general, the cancellation is small,
   so that the new additions in (3) will occur only in a small part
   of the accumulator. But see below about long carry propagation.

7. Copy the rounded significand to the destination; one goal is to free
   the accumulator in case there is work to do in (8).

8. If the sign of the error term is needed to determine the returned
   result (correctly rounded sum and ternary value), then reiterate
   at (3) to compute it using a specific window as if sq were 0, i.e.
   around maxexp + cq to maxexp - dq.

9. Correct the significand if need be (+ or - 1 ulp), determine the
   exponent, and exit with the correct ternary value.

Accumulator at the first iteration:

                 [------]---------------------------------------]
                    cq  ^- maxexp        sq + dq        minexp -^

If there is a second iteration with maxexp2 = minexp - 3 and a shift
of 20 bits:
                                             <---- 20 zeros ---->
                 [---------------------------0000000000000000000]
                                       maxexp -^        minexp -^

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
appears that an important cancellation occurred, and another step (3)
is needed. Since u1 was truncated, the new maxexp will be minexp, i.e.
[1]. At iteration 2, minexp is determined to be [2]; thus a part of u1,
a part of u2, and u3 are taken into account for the truncated sum. Now
assume that on this example, the error is small enough, but its sign is
unknown. Thus another step (3) is needed, with the conditions of (7).
Since no numbers were truncated at the previous iteration, maxexp is
the maximum exponent of the remaining numbers, here the one of u4, and
minexp is determined to be [3]. Assume that the sign of the error can
be determined now, so that we can return the rounded result with the
ternary value.

As a bonus, this will also solve overflow, underflow and normalization
issues, since everything is done in fixed point and the output exponent
will be considered only at the end (early overflow detection could also
be done).

A limb L = *zp in memory will generally contain a part of a significand.
One can define its exponent ze, such that the actual value of this limb
is L * 2^(ze-GMP_NUMB_BITS), i.e. ze is its exponent where the limb is
regarded as a number in [1/2,1[. If an array of limbs zp[] is regarded
as a significand in [1/2,1[, then the exponent of its actual value is
also ze.

Variables:
  tp: pointer to a temporary area that will contain a shifted value.
  wp: pointer to the accumulator.
  ts: the size of the temporary area, in limbs; then wp = tp + ts.
  ws: the size of the accumulator, in limbs.
  rn: number n' of inputs that are regular numbers (regular inputs).
  logn: ceil(log2(rn)).
  cq: logn + 1.
  sq: output precision (precision of the sum).
  dq: > logn, such that cq + sq + dq is a multiple of GMP_NUMB_BITS.
  maxexp: the maximum exponent of the bit-window of the inputs that is
          taken into account (for the current iteration), excluded.
  minexp: the minimum exponent of the bit-window of the inputs that is
          taken into account (for the current iteration), included.

Note 1: Data locality can be improved after the first iteration if the
shifted values are stored at the end of the temporary area instead of
the beginning. The reason is that only the least significant part of
the accumulator will be used once a good approximate value of the sum
is known, and the accumulator lies just after the temporary area. But
the gain would probably not be noticeable in practice.

Note 2: At step (4), it was considered to determine the number of
cancelled limbs instead of the number of cancelled bits in order to
avoid a non-trivial shift at step (6), making this step a bit faster.
However this choice would have required a larger value of dq (with
an increment of up to something like GMP_NUMB_BITS - 1) at the first
iteration to be able to handle a number of cancelled bits just below
GMP_NUMB_BITS (similar problem to the one of small leading digit in
high-radix floating-point representations).

Note 3: Compared with Ziv's strategy, the issue here is not really
exact values that are very close to a breakpoint, but cancellation.
Moreover, we do not need to recompute everything at each iteration.
The issue with the value very close to a breakpoint actually occurs
at step (8); however, contrary to the usual case, we do not want to
reiterate with more precision as this could take too much time and
memory. Indeed, due to a possible "hole" between the inputs, the
distance between the exact value and a breakpoint can be extremely
small.

*** To be considered in future versions ***
It seems that carry propagation (mpn_add_1 & mpn_sub_1 in the code) is
most often limited. But consider the following cases, where all inputs
have the minimal precision 2, and the output precision is p:
  u0 = 1
  u_i = (-1)^i * 2^(-p) for i > 0
Here long carry propagation will occur for each addition of the initial
iteration, so that the complexity will be O(n*p) instead of O(n+p) if
we choose to delay carry propagation (however such a choice may slower
the average case and take more memory, such as around 3*p instead of
2*p).
When a new iteration is needed due to cancellation, a second accumulator
was considered in some early version of the algorithm: the temporary
results of the computations during the new iteration would be stored in
this second accumulator, which would generally be small, thus limiting
carry propagation; this method is actually equivalent to delaying carry
propagation. It could help in some cases, such as:
  u0 = 2^q with some q > 0
  u1 = 1
  u2 = -2^q
  u_i = (-1)^i * 2^(-p) for i > 2
but such examples are very specific cases, and as seen with the first
example, a better way must be chosen if avoiding long carry propagation
is regarded as important (in all cases). Moreover, while the use of two
accumulators does not take much more memory (since both accumulators can
occupy the same area, with a flexible limit between them), it probably
makes the code a bit more complex, and noticeably slower if n is small.

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
      mp_limb_t *tp;  /* pointer to a temporary area */
      mp_limb_t *wp;  /* pointer to the accumulator */
      mp_size_t ts;   /* size of the temporary area, in limbs */
      mp_size_t ws;   /* size of the accumulator, in limbs */
      mpfr_exp_t maxexp;
      unsigned long i, rn;
      int logn;       /* ceil(log2(rn)) */
      int cq, dq;
      mpfr_prec_t sq;
      int inex;
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

      /* Step 2: set up some variables and the accumulator. */

      /* rn is the number of regular inputs (the singular ones will be
         ignored). Compute logn = ceil(log2(rn)). */
      logn = MPFR_INT_CEIL_LOG2 (rn);

      MPFR_LOG_MSG (("logn=%d maxexp=%" MPFR_EXP_FSPEC "d\n",
                     logn, (mpfr_eexp_t) maxexp));

      sq = MPFR_GET_PREC (sum);
      cq = logn + 1;

      /* First determine the size of the accumulator. */
      ws = MPFR_PREC2LIMBS (cq + sq + logn + 1);

      /* We now choose dq so that cq + sq + dq = accumulator bitsize. */
      dq = (mpfr_prec_t) ws * GMP_NUMB_BITS - (cq + sq);

      /* An input block will have up to sq + dq bits, and its shifted
         value (to be correctly aligned) may take GMP_NUMB_BITS - 1
         additional bits. */
      ts = MPFR_PREC2LIMBS (sq + dq + GMP_NUMB_BITS - 1);

      MPFR_TMP_MARK (marker);

      /* TODO: one may need a bit more memory later for Step 6.
         Should it be allocated here? */
      tp = MPFR_TMP_LIMBS_ALLOC (ts + ws);
      wp = tp + ts;

      MPN_ZERO (wp, ws);  /* zero the accumulator */

      while (1)
        {
          mpfr_exp_t minexp = maxexp - (sq + dq);
          mpfr_exp_t maxexp2 = MPFR_EXP_MIN;  /* < any valid exponent */

          /* Step 3: compute the truncated sum. */

          for (i = 0; i < n; i++)
            if (! MPFR_IS_SINGULAR (p[i]))
              {
                mp_limb_t *vp;
                mp_size_t vs;
                mpfr_exp_t pe, vd;
                mpfr_prec_t pq;

                pe = MPFR_GET_EXP (p[i]);
                pq = MPFR_GET_PREC (p[i]);

                vp = MPFR_MANT (p[i]);
                vs = MPFR_PREC2LIMBS (pq);
                vd = pe - vs * GMP_NUMB_BITS - minexp;
                /* vd is the exponent of the least significant represented
                   bit of p[i] (including the trailing bits, whose value
                   is 0) minus the exponent of the least significant bit
                   of the accumulator. */

                if (vd < 0)
                  {
                    mp_size_t vds;
                    int tr;

                    /* This covers the following cases:
                     *     [-+- accumulator ---]
                     *   [---|----- p[i] ------|--]
                     *       |   [----- p[i] --|--]
                     *       |                 |[----- p[i] -----]
                     *       |                 |    [----- p[i] -----]
                     *     maxexp           minexp
                     */

                    if (pe <= minexp)
                      {
                        /* p[i] is entirely after the LSB of the accumulator,
                           so that it will be ignored at this iteration. */
                        if (pe > maxexp2)
                          maxexp2 = pe;
                        continue;
                      }

                    /* If some significant bits of p[i] are after the LSB
                       of the accumulator, then maxexp2 will necessarily
                       be minexp. */
                    if (MPFR_LIKELY (pe - pq < minexp))
                      maxexp2 = minexp;

                    /* We need to ignore the least |vd| significant bits
                       of p[i]. First, let's ignore the least
                       vds = |vd| / GMP_NUMB_BITS limbs. */
                    vd = - vd;
                    vds = vd / GMP_NUMB_BITS;
                    vs -= vds;
                    MPFR_ASSERTD (vs > 0);  /* see pe <= minexp test above */
                    vp += vds;
                    vd -= vds * GMP_NUMB_BITS;
                    MPFR_ASSERTD (vd >= 0 && vd < GMP_NUMB_BITS);

                    if (pe > maxexp)
                      {
                        vs -= (pe - maxexp) / GMP_NUMB_BITS;
                        tr = (pe - maxexp) % GMP_NUMB_BITS;
                      }
                    else
                      tr = 0;

                    if (vd != 0)
                      {
                        MPFR_ASSERTD (vs <= ts);
                        mpn_rshift (tp, vp, vs, vd);
                        vp = tp;
                        tr += vd;
                        if (tr >= GMP_NUMB_BITS)
                          {
                            vs--;
                            tr -= GMP_NUMB_BITS;
                          }
                        MPFR_ASSERTD (tr >= 0 && tr < GMP_NUMB_BITS);
                        if (tr != 0)
                          {
                            tp[vs-1] &=
                              MPFR_LIMB_MASK (GMP_NUMB_BITS - tr);
                            tr = 0;
                          }
                        /* Truncation has now been taken into account. */
                        MPFR_ASSERTD (tr == 0);
                      }

                    MPFR_ASSERTD (vs <= ws);

                    if (tr != 0)
                      {
                        /* We can't truncate the most significant limb of
                           the input. So, let's ignore it now. It will be
                           taken into account after the addition. */
                        vs--;
                      }

                    if (MPFR_IS_POS (p[i]))
                      {
                        mp_limb_t carry;

                        carry = mpn_add_n (wp, wp, vp, vs);
                        if (tr != 0)
                          carry += vp[vs] &
                            MPFR_LIMB_MASK (GMP_NUMB_BITS - tr);
                        if (carry != 0 && vs < ws)
                          mpn_add_1 (wp + vs, wp + vs, ws - vs, carry);
                      }
                    else
                      {
                        mp_limb_t borrow;

                        borrow = mpn_sub_n (wp, wp, vp, vs);
                        if (tr != 0)
                          borrow += vp[vs] &
                            MPFR_LIMB_MASK (GMP_NUMB_BITS - tr);
                        if (borrow != 0 && vs < ws)
                          mpn_sub_1 (wp + vs, wp + vs, ws - vs, borrow);
                      }
                  }
                else  /* vd >= 0 */
                  {
                    mp_limb_t *dp;
                    mp_size_t ds, vds;
                    int tr;

                    /* This covers the following cases:
                     *               [-+- accumulator ---]
                     *   [- p[i] -]    |
                     *             [---|-- p[i] ------]  |
                     *          [------|-- p[i] ---------]
                     *                 |   [- p[i] -]    |
                     *               maxexp           minexp
                     */

                    /* We need to ignore the least vd significant bits
                       of the accumulator. First, let's ignore the least
                       vds = vd / GMP_NUMB_BITS limbs. -> (dp,ds) */
                    vds = vd / GMP_NUMB_BITS;
                    ds = ws - vds;
                    if (ds <= 0)
                      continue;
                    dp = wp + vds;
                    vd -= vds * GMP_NUMB_BITS;
                    MPFR_ASSERTD (vd >= 0 && vd < GMP_NUMB_BITS);

                    if (pe > maxexp)
                      {
                        vs -= (pe - maxexp) / GMP_NUMB_BITS;
                        tr = (pe - maxexp) % GMP_NUMB_BITS;
                        if (tr > vd || (vd != 0 && tr == vd))
                          {
                            vs--;
                            tr -= GMP_NUMB_BITS;
                          }
                      }
                    else
                      tr = 0;
                    MPFR_ASSERTD (tr >= 1 - GMP_NUMB_BITS && tr <= vd);

                    if (vd != 0)
                      {
                        MPFR_ASSERTD (vs + 1 <= ts);
                        tp[vs] = mpn_lshift (tp, vp, vs, vd);
                        MPFR_ASSERTD (vd - tr > 0);
                        if (vd - tr < GMP_NUMB_BITS)
                          tp[vs] &= MPFR_LIMB_MASK (vd - tr);
                        vp = tp;
                        tr = 0;
                      }

                    if (MPFR_IS_POS (p[i]))
                      {
                        mp_limb_t carry;

                        carry = mpn_add_n (dp, dp, vp, vs);
                        if (tr != 0)
                          carry += vp[vs] & MPFR_LIMB_MASK (- tr);
                        if (carry != 0 && vs < ds)
                          mpn_add_1 (dp + vs, dp + vs, ds - vs, carry);
                      }
                    else
                      {
                        mp_limb_t borrow;

                        borrow = mpn_sub_n (dp, dp, vp, vs);
                        if (tr != 0)
                          borrow += vp[vs] & MPFR_LIMB_MASK (- tr);
                        if (borrow != 0 && vs < ds)
                          mpn_sub_1 (dp + vs, dp + vs, ds - vs, borrow);
                      }
                  }
              }  /* end of the iteration (Step 3) */

          {
            mpfr_prec_t cancel;  /* number of cancelled bits */
            mp_size_t wi;        /* index in the accumulator */
            mp_limb_t msl;       /* most significant limb */
            mpfr_exp_t e;        /* temporary exponent of the result */

            /* Step 4: determine the number of cancelled bits. */

            cancel = 0;
            wi = ws - 1;
            MPFR_ASSERTD (wi >= 0);
            msl = wp[wi];

            /* Limbs whose bits are identical (000...00 or 111...11). */
            if (MPFR_UNLIKELY (msl == MPFR_LIMB_ZERO || msl == MPFR_LIMB_MAX))
              {
                while (wi >= 0 && wp[wi] == msl)
                  {
                    cancel += GMP_NUMB_BITS;
                    wi--;
                  }

                if (wi < 0 && msl == MPFR_LIMB_ZERO)
                  {
                    /* Step 5: the truncated sum is zero. Reiterate with
                     * maxexp = maxexp2. Note: we do not need to zero the
                     * accumulator since it is already 0 in this case.
                     */
                    maxexp = maxexp2;
                    continue;
                  }
              }

            /* Let's count the number of identical leading bits of
               the next limb, if there is one. */
            if (MPFR_LIKELY (wi >= 0))
              {
                int cnt;

                msl = wp[wi];
                if (msl & MPFR_LIMB_HIGHBIT)
                  msl ^= MPFR_LIMB_MAX;
                count_leading_zeros (cnt, msl);
                cancel += cnt;
              }

            /* Step 6 */

            e = maxexp + cq - cancel;

            /* The truncated sum is in the binade [2^(e-1),2^e]
               (closed on both ends due to two's complement).
               The error is less than 2^(minexp+logn). */


          }
        }

      MPFR_TMP_FREE (marker);
      return mpfr_check_range (sum, inex, rnd);
    }
}
