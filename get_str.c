/* mpfr_get_str -- output a floating-point number to a string

Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.
Contributed by Alain Delplanque and Paul Zimmermann.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <string.h> /* For strlen */
#include <limits.h> /* For CHAR_BIT, LONG_MAX and LONG_MIN */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

static int mpfr_get_str_aux (char *const, mp_exp_t *const, mp_limb_t *const,
                       mp_size_t, mp_exp_t, long, int, size_t, mp_rnd_t);

static const char num_to_text[] = "0123456789abcdefghijklmnopqrstuvwxyz";

/* copy most important limbs of {op, n2} in {rp, n1} */
/* if n1 > n2 put 0 in low limbs of {rp, n1} */
#define MPN_COPY2(rp, n1, op, n2) \
  if ((n1) <= (n2)) \
    { \
      MPN_COPY ((rp), (op) + (n2) - (n1), (n1)); \
    } \
  else \
    { \
      MPN_COPY ((rp) + (n1) - (n2), (op), (n2)); \
      MPN_ZERO ((rp), (n1) - (n2)); \
    }

#define MPFR_ROUND_FAILED 3

/* Input: an approximation r*2^f of an real Y, with |r*2^f-Y| <= 2^(e+f).
   Returns if possible in the string s the mantissa corresponding to
   the integer nearest to Y, within the direction rnd, and returns the
   exponent in exp.
   n is the number of limbs of r.
   e represents the maximal error in the approximation of Y
      (e < 0 iff the approximation is exact, i.e. r*2^f = Y).
   b is the wanted base (2 <= b <= 36).
   m is the number of wanted digits in the mantissa.
   rnd is the rounding mode.
   It is assumed that b^(m-1) <= Y < b^(m+1), thus the returned value
   satisfies b^(m-1) <= rnd(Y) < b^(m+1).

   Rounding may fail for two reasons:
   - the error is too large to determine the integer N nearest to Y
   - either the number of digits of N in base b is too large (m+1),
     N=2*N1+(b/2) and the rounding mode is to nearest. This can
     only happen when b is even.

   Return value:
   - the direction of rounding (-1, 0, 1) if rounding is possible
   - -MPFR_ROUND_FAILED if rounding not possible because m+1 digits
   - MPFR_ROUND_FAILED otherwise (too large error)
*/
static int
mpfr_get_str_aux (char *const str, mp_exp_t *const exp, mp_limb_t *const r,
                  mp_size_t n, mp_exp_t f, long e, int b, size_t m,
                  mp_rnd_t rnd)
{
  int dir;                  /* direction of the rounded result */
  mp_limb_t ret = 0;        /* possible carry in addition */
  mp_size_t i0, j0;         /* number of limbs and bits of Y */
  unsigned char *str1;      /* string of m+2 characters */
  size_t size_s1;           /* length of str1 */
  mp_rnd_t rnd1;
  size_t i;
  int exact = (e < 0);
  MPFR_TMP_DECL(marker);

  /* if f > 0, then the maximal error 2^(e+f) is larger than 2 so we can't
     determine the integer Y */
  MPFR_ASSERTN(f <= 0);
  /* if f is too small, then r*2^f is smaller than 1 */
  MPFR_ASSERTN(f > (-n * BITS_PER_MP_LIMB));

  MPFR_TMP_MARK(marker);

  /* R = 2^f sum r[i]K^(i)
     r[i] = (r_(i,k-1)...r_(i,0))_2
     R = sum r(i,j)2^(j+ki+f)
     the bits from R are referenced by pairs (i,j) */

  /* check if is possible to round r with rnd mode
     where |r*2^f-Y| <= 2^(e+f)
     the exponent of R is: f + n*BITS_PER_MP_LIMB
     we must have e + f == f + n*BITS_PER_MP_LIMB - err
     err = n*BITS_PER_MP_LIMB - e
     R contains exactly -f bits after the integer point:
     to determine the nearest integer, we thus need a precision of
     n * BITS_PER_MP_LIMB + f */

  if (exact || mpfr_can_round_raw (r, n, (mp_size_t) 1,
            n * BITS_PER_MP_LIMB - e, GMP_RNDN, rnd, n * BITS_PER_MP_LIMB + f))
    {
      /* compute the nearest integer to R */

      /* bit of weight 0 in R has position j0 in limb r[i0] */
      i0 = (-f) / BITS_PER_MP_LIMB;
      j0 = (-f) % BITS_PER_MP_LIMB;

      ret = mpfr_round_raw (r + i0, r, n * BITS_PER_MP_LIMB, 0,
                            n * BITS_PER_MP_LIMB + f, rnd, &dir);
      MPFR_ASSERTD(dir != MPFR_ROUND_FAILED);

      /* warning: mpfr_round_raw_generic returns MPFR_EVEN_INEX (2) or
         -MPFR_EVEN_INEX (-2) in case of even rounding */

      if (ret) /* Y is a power of 2 */
        {
          if (j0)
            r[n - 1] = MPFR_LIMB_HIGHBIT >> (j0 - 1);
          else /* j0=0, necessarily i0 >= 1 otherwise f=0 and r is exact */
            {
              r[n - 1] = ret;
              r[--i0] = 0; /* set to zero the new low limb */
            }
        }
      else /* shift r to the right by (-f) bits (i0 already done) */
        {
          if (j0)
            mpn_rshift (r + i0, r + i0, n - i0, j0);
        }

      /* now the rounded value Y is in {r+i0, n-i0} */

      /* convert r+i0 into base b */
      str1 = (unsigned char*) MPFR_TMP_ALLOC (m + 3); /* need one extra character for mpn_get_str */
      size_s1 = mpn_get_str (str1, b, r + i0, n - i0);

      /* round str1 */
      MPFR_ASSERTN(size_s1 >= m);
      *exp = size_s1 - m; /* number of superfluous characters */

      /* if size_s1 = m + 2, necessarily we have b^(m+1) as result,
         and the result will not change */

      /* so we have to double-round only when size_s1 = m + 1 and
         (i) the result is inexact
         (ii) or the last digit is non-zero */
      if ((size_s1 == m + 1) && ((dir != 0) || (str1[size_s1 - 1] != 0)))
        {
          /* rounding mode */
          rnd1 = rnd;

          /* round to nearest case */
          if (rnd == GMP_RNDN)
            {
              if (2 * str1[size_s1 - 1] == b)
                {
                  if (dir == 0 && exact) /* exact: even rounding */
                    {
                      rnd1 = ((str1[size_s1-2] & 1) == 0)
                        ? GMP_RNDD : GMP_RNDU;
                    }
                  else
                    {
                      /* otherwise we cannot round correctly: for example
                         if b=10, we might have a mantissa of
                         xxxxxxx5.00000000 which can be rounded to nearest
                         to 8 digits but not to 7 */
                      dir = -MPFR_ROUND_FAILED;
                      MPFR_ASSERTD(dir != MPFR_EVEN_INEX);
                      goto free_and_return;
                    }
                }
              else if (2 * str1[size_s1 - 1] < b)
                rnd1 = GMP_RNDD;
              else
                rnd1 = GMP_RNDU;
            }

          /* now rnd1 is either GMP_RNDD or GMP_RNDZ -> truncate
                             or GMP_RDNU -> round towards infinity */

          /* round away from zero */
          if (rnd1 == GMP_RNDU)
            {
              if (str1[size_s1 - 1] != 0)
                {
                  /* the carry cannot propagate to the whole string, since
                     Y = x*b^(m-g) < 2*b^m <= b^(m+1)-b
                     where x is the input float */
                  MPFR_ASSERTN(size_s1 >= 2);
                  i = size_s1 - 2;
                  while (str1[i] == b - 1)
                    {
                      MPFR_ASSERTD(i > 0);
                      str1[i--] = 0;
                    }
                  str1[i]++;
                }
              dir = 1;
            }
          /* round toward zero (truncate) */
          else
            dir = -1;
        }

      /* copy str1 into str and convert to ASCII */
      for (i = 0; i < m; i++)
        str[i] = num_to_text[(int) str1[i]];
        str[m] = 0;
    }
  /* mpfr_can_round_raw failed: rounding is not possible */
  else
    {
      dir = MPFR_ROUND_FAILED; /* should be different from MPFR_EVEN_INEX */
      MPFR_ASSERTD(dir != MPFR_EVEN_INEX);
    }

 free_and_return:
  MPFR_TMP_FREE(marker);

  return dir;
}

/* mpfr_l2b[b-2][0] is a 23-bit upper approximation to log(b)/log(2),
   mpfr_l2b[b-2][1] is a 76-bit upper approximation to log(2)/log(b),
   if not null. The array is initialized with null pointers. */
mpfr_ptr MPFR_THREAD_ATTR mpfr_l2b[BASE_MAX-1][2] = { { NULL } };

/* returns ceil(e * log2(b)^((-1)^i)), or ... + 1 */
static mp_exp_t
ceil_mul (mp_exp_t e, int beta, int i)
{
  mpfr_ptr p;
  mpfr_t t;
  mp_exp_t r;

  p = mpfr_l2b[beta-2][i];
  if (p == NULL)
    {
      mpfr_l2b[beta-2][i] = p =
        (mpfr_ptr) (*__gmp_allocate_func) (sizeof (mpfr_t));
      if (i == 0)
        {
          mpfr_init2 (p, 23);
          mpfr_set_ui (p, beta, GMP_RNDU);
          mpfr_log2 (p, p, GMP_RNDU);
        }
      else
        {
          mpfr_init2 (p, 77);
          mpfr_set_ui (p, beta, GMP_RNDD);
          mpfr_log2 (p, p, GMP_RNDD);
          mpfr_ui_div (p, 1, p, GMP_RNDU);
        }
    }
  mpfr_init2 (t, sizeof (mp_exp_t) * CHAR_BIT);
  mpfr_set_exp_t (t, e, GMP_RNDU);
  mpfr_mul (t, t, p, GMP_RNDU);
  r = mpfr_get_exp_t (t, GMP_RNDU);
  mpfr_clear (t);
  return r;
}

/* prints the mantissa of x in the string s, and writes the corresponding
   exponent in e.
   x is rounded with direction rnd, m is the number of digits of the mantissa,
   b is the given base (2 <= b <= 36).

   Return value:
   if s=NULL, allocates a string to store the mantissa, with
   m characters, plus a final '\0', plus a possible minus sign
   (thus m+1 or m+2 characters).

   Important: when you call this function with s=NULL, don't forget to free
   the memory space allocated, with free(s, strlen(s)).
*/
char*
mpfr_get_str (char *s, mp_exp_t *e, int b, size_t m, mpfr_srcptr x, mp_rnd_t rnd)
{
  int exact;                      /* exact result */
  mp_exp_t exp, g;
  mp_exp_t prec; /* precision of the computation */
  long err;
  mp_limb_t *a;
  mp_exp_t exp_a;
  mp_limb_t *result;
  mp_limb_t *xp, *x1;
  mp_limb_t *reste;
  size_t nx, nx1;
  size_t n, i;
  char *s0;
  int neg;
  int ret; /* return value of mpfr_get_str_aux */
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_TMP_DECL(marker);

  /* if exact = 1 then err is undefined */
  /* otherwise err is such that |x*b^(m-g)-a*2^exp_a| < 2^(err+exp_a) */

  /* is the base valid? */
  if (b < 2 || b > 36)
    return NULL;

  if (MPFR_UNLIKELY (MPFR_IS_NAN (x)))
    {
      if (s == NULL)
        s = (char *) (*__gmp_allocate_func) (6);
      strcpy (s, "@NaN@");
      return s;
    }

  neg = MPFR_SIGN(x) < 0; /* 0 if positive, 1 if negative */

  if (MPFR_UNLIKELY (MPFR_IS_INF (x)))
    {
      if (s == NULL)
        s = (char *) (*__gmp_allocate_func) (neg + 6);
      strcpy (s, (neg) ? "-@Inf@" : "@Inf@");
      return s;
    }

  MPFR_SAVE_EXPO_MARK (expo);  /* needed for ceil_mul (at least) */

  if (m == 0)
    {

      /* take at least 1 + ceil(n*log(2)/log(b)) digits, where n is the
         number of bits of the mantissa, to ensure back conversion from
         the output gives the same floating-point.

         Warning: if b = 2^k, this may be too large. The worst case is when
         the first base-b digit contains only one bit, so we get
         1 + ceil((n-1)/k) = 2 + floor((n-2)/k) instead.
      */
      m = 1 + ceil_mul (IS_POW2(b) ? MPFR_PREC(x) - 1: MPFR_PREC(x), b, 1);
      if (m < 2)
        m = 2;
    }

  /* the code below for non-power-of-two bases works for m=1 */
  MPFR_ASSERTN (m >= 2 || (IS_POW2(b) == 0 && m >= 1));

  /* x is a floating-point number */

  if (MPFR_IS_ZERO(x))
    {
      if (s == NULL)
        s = (char*) (*__gmp_allocate_func) (neg + m + 1);
      s0 = s;
      if (neg)
        *s++ = '-';
      memset (s, '0', m);
      s[m] = '\0';
      *e = 0; /* a bit like frexp() in ISO C99 */
      MPFR_SAVE_EXPO_FREE (expo);
      return s0; /* strlen(s0) = neg + m */
    }

  if (s == NULL)
    s = (char*) (*__gmp_allocate_func) (neg + m + 1);
  s0 = s;
  if (neg)
    *s++ = '-';

  xp = MPFR_MANT(x);

  if (IS_POW2(b))
    {
      int pow2;
      mp_exp_t f, r;
      mp_limb_t *x1;
      mp_size_t nb;
      int inexp;

      count_leading_zeros (pow2, (mp_limb_t) b);
      pow2 = BITS_PER_MP_LIMB - pow2 - 1; /* base = 2^pow2 */

      /* set MPFR_EXP(x) = f*pow2 + r, 1 <= r <= pow2 */
      f = (MPFR_GET_EXP (x) - 1) / pow2;
      r = MPFR_GET_EXP (x) - f * pow2;
      if (r <= 0)
        {
          f --;
          r += pow2;
        }

      /* the first digit will contain only r bits */
      prec = (m - 1) * pow2 + r; /* total number of bits */
      n = (prec - 1) / BITS_PER_MP_LIMB + 1;

      MPFR_TMP_MARK (marker);
      x1 = (mp_limb_t*) MPFR_TMP_ALLOC((n + 1) * sizeof (mp_limb_t));
      nb = n * BITS_PER_MP_LIMB - prec;
      /* round xp to the precision prec, and put it into x1
         put the carry into x1[n] */
      if ((x1[n] = mpfr_round_raw (x1, xp, MPFR_PREC(x),
                                  MPFR_IS_STRICTNEG(x),
                                   prec, rnd, &inexp)))
        {
          /* overflow when rounding x: x1 = 2^prec */
          if (r == pow2)    /* prec = m * pow2,
                               2^prec will need (m+1) digits in base 2^pow2 */
            {
              /* divide x1 by 2^pow2, and increase the exponent */
              mpn_rshift (x1, x1, n + 1, pow2);
              f ++;
            }
          else /* 2^prec needs still m digits, but x1 may need n+1 limbs */
            n ++;
        }

      /* it remains to shift x1 by nb limbs to the right, since mpn_get_str
         expects a right-normalized number */
      if (nb != 0)
        {
          mpn_rshift (x1, x1, n, nb);
          /* the most significant word may be zero */
          if (x1[n - 1] == 0)
            n --;
        }

      mpn_get_str ((unsigned char*) s, b, x1, n);
      for (i=0; i<m; i++)
        s[i] = num_to_text[(int) s[i]];
      s[m] = 0;

      /* the exponent of s is f + 1 */
      *e = f + 1;

      MPFR_TMP_FREE(marker);
      MPFR_SAVE_EXPO_FREE (expo);
      return (s0);
    }

  /* if x < 0, reduce to x > 0 */
  if (neg)
    rnd = MPFR_INVERT_RND(rnd);

  g = ceil_mul (MPFR_GET_EXP (x) - 1, b, 1);
  exact = 1;
  prec = ceil_mul (m, b, 0) + 1;
  exp = ((mp_exp_t) m < g) ? g - (mp_exp_t) m : (mp_exp_t) m - g;
  prec += MPFR_INT_CEIL_LOG2 (prec); /* number of guard bits */
  if (exp != 0) /* add maximal exponentiation error */
    prec += 3 * (mp_exp_t) MPFR_INT_CEIL_LOG2 (exp);

  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      MPFR_TMP_MARK(marker);

      exact = 1;

      /* number of limbs */
      n = 1 + (prec - 1) / BITS_PER_MP_LIMB;

      /* a will contain the approximation of the mantissa */
      a = (mp_limb_t*) MPFR_TMP_ALLOC (n * sizeof (mp_limb_t));

      nx = 1 + (MPFR_PREC(x) - 1) / BITS_PER_MP_LIMB;

      if ((mp_exp_t) m == g) /* final exponent is 0, no multiplication or
                                division to perform */
        {
          if (nx > n)
            exact = mpn_scan1 (xp, 0) >= (nx - n) * BITS_PER_MP_LIMB;
          err = !exact;
          MPN_COPY2 (a, n, xp, nx);
          exp_a = MPFR_GET_EXP (x) - n * BITS_PER_MP_LIMB;
        }
      else if ((mp_exp_t) m > g) /* we have to multiply x by b^exp */
        {
          /* a2*2^exp_a =  b^e */
          err = mpfr_mpn_exp (a, &exp_a, b, exp, n);
          /* here, the error on a is at most 2^err ulps */
          exact = (err == -1);

          /* x = x1*2^(n*BITS_PER_MP_LIMB) */
          x1 = (nx >= n) ? xp + nx - n : xp;
          nx1 = (nx >= n) ? n : nx; /* nx1 = min(n, nx) */

          /* test si exact */
          if (nx > n)
            exact = (exact &&
                     ((mpn_scan1 (xp, 0) >= (nx - n) * BITS_PER_MP_LIMB)));

          /* we loose one more bit in the multiplication,
             except when err=0 where we loose two bits */
          err = (err <= 0) ? 2 : err + 1;

          /* result = a * x */
          result = (mp_limb_t*) MPFR_TMP_ALLOC ((n + nx1) * sizeof (mp_limb_t));
          mpn_mul (result, a, n, x1, nx1);
          exp_a += MPFR_GET_EXP (x);
          if (mpn_scan1 (result, 0) < (nx1 * BITS_PER_MP_LIMB))
            exact = 0;

          /* normalize a and truncate */
          if ((result[n + nx1 - 1] & MPFR_LIMB_HIGHBIT) == 0)
            {
              mpn_lshift (a, result + nx1, n , 1);
              a[0] |= result[nx1 - 1] >> (BITS_PER_MP_LIMB - 1);
              exp_a --;
            }
          else
            MPN_COPY (a, result + nx1, n);
        }
      else
        {
          /* a2*2^exp_a =  b^e */
          err = mpfr_mpn_exp (a, &exp_a, b, exp, n);
          exact = (err == -1);

          /* allocate memory for x1, result and reste */
          x1 = (mp_limb_t*) MPFR_TMP_ALLOC (2 * n * sizeof (mp_limb_t));
          result = (mp_limb_t*) MPFR_TMP_ALLOC ((n + 1) * sizeof (mp_limb_t));
          reste = (mp_limb_t*) MPFR_TMP_ALLOC (n * sizeof (mp_limb_t));

          /* initialize x1 = x */
          MPN_COPY2 (x1, 2 * n, xp, nx);
          if ((exact) && (nx > 2 * n) &&
              (mpn_scan1 (xp, 0) < (nx - 2 * n) * BITS_PER_MP_LIMB))
            exact = 0;

          /* result = x / a */
          mpn_tdiv_qr (result, reste, 0, x1, 2 * n, a, n);
          exp_a = MPFR_GET_EXP (x) - exp_a - 2 * n * BITS_PER_MP_LIMB;

          /* test if division was exact */
          if (exact)
            exact = mpn_popcount (reste, n) == 0;

          /* normalize the result and copy into a */
          if (result[n] == 1)
            {
              mpn_rshift (a, result, n, 1);
              a[n - 1] |= MPFR_LIMB_HIGHBIT;;
              exp_a ++;
            }
          else
            MPN_COPY (a, result, n);

          err = (err == -1) ? 2 : err + 2;
        }

      /* check if rounding is possible */
      if (exact)
        err = -1;
      ret = mpfr_get_str_aux (s, e, a, n, exp_a, err, b, m, rnd);
      if (ret == MPFR_ROUND_FAILED)
        {
          /* too large error: increment the working precision */
          MPFR_ZIV_NEXT (loop, prec);
        }
      else if (ret == -MPFR_ROUND_FAILED)
        {
          /* too many digits in mantissa: exp = |m-g| */
          if ((mp_exp_t) m > g) /* exp = m - g, multiply by b^exp */
            {
              g++;
              exp --;
            }
          else /* exp = g - m, divide by b^exp */
            {
              g++;
              exp ++;
            }
        }
      else
        break;

      MPFR_TMP_FREE(marker);
    }
  MPFR_ZIV_FREE (loop);

  *e += g;

  MPFR_TMP_FREE(marker);
  MPFR_SAVE_EXPO_FREE (expo);
  return s0;
}

void mpfr_free_str (char *str)
{
   (*__gmp_free_func) (str, strlen (str) + 1);
}
