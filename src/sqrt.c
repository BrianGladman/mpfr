/* mpfr_sqrt -- square root of a floating-point number

Copyright 1999-2016 Free Software Foundation, Inc.
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

/*
*** Note ***

  The code of mpfr_sqrt1 and/or functions it calls depends on
  implementation-defined features of the C standard. See lines with a
  cast to mp_limb_signed_t, associated with a shift to the right '>>'.
  Such features are known to behave as the code below expects with GCC,
  according to the GCC manual (excerpt from the GCC 4.9.3 manual):

    4 C Implementation-defined behavior
    ***********************************

    4.5 Integers
    ============

    * 'The result of, or the signal raised by, converting an integer to
      a signed integer type when the value cannot be represented in an
      object of that type (C90 6.2.1.2, C99 and C11 6.3.1.3).'

      For conversion to a type of width N, the value is reduced modulo
      2^N to be within range of the type; no signal is raised.

    * 'The results of some bitwise operations on signed integers (C90
      6.3, C99 and C11 6.5).'

      Bitwise operators act on the representation of the value including
      both the sign and value bits, where the sign bit is considered
      immediately above the highest-value value bit.  Signed '>>' acts
      on negative numbers by sign extension.

  It is not known whether it works with other compilers. Thus this code
  is currently enabled only when __GNUC__ is defined (which includes
  compilers that declare a compatibility with GCC). A configure test
  might be an alternative solution (but without any guarantee, in case
  the result may also depend on the context).

  Warning! The right shift of a negative value corresponds to an integer
  division by a power of two, with rounding toward negative.

  TODO: Complete the comments when a right shift of a negative value
  may be involved, so that the rounding toward negative appears in the
  proof. There has been at least an error with a proof of a bound!
*/

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#if !defined(MPFR_GENERIC_ABI) && (GMP_NUMB_BITS == 32 || GMP_NUMB_BITS == 64)

/* The tables T1[] and T2[] below were generated using the Sage code below,
   with T1,T2 = bipartite(4,4,4,12,16). Note: we would get a slightly smaller
   error using an approximation of the form T1[a,b] * (1 + T2[a,b]), but this
   would make the code more complex, and the approximation would not always fit
   on p2 bits (assuming p2 >= p1).
# bi-partite method, using a table T1 of p1 bits, and T2 of p2 bits
# approximate sqrt(a:b:c) with T1[a,b] + T2[a,c]
def bipartite(pa,pb,pc,p1,p2):
   T1 = dict()
   T2 = dict()
   maxerr = maxnum = 0
   for a in range(2^(pa-2),2^pa):
      for b in range(2^pb):
         A1 = (a*2^pb+b)/2^(pa+pb)
         A2 = (a*2^pb+b+1)/2^(pa+pb)
         X = (1/sqrt(A1) + 1/sqrt(A2))/2
         X = round(X*2^p1)/2^p1
         T1[a,b] = X*2^p1
         maxnum = max(maxnum, abs(T1[a,b]))
         # print "a=", a, "b=", b, "T1=", x
         maxerr = max(maxerr, n(abs(1/sqrt(A1) - X)))
         maxerr = max(maxerr, n(abs(1/sqrt(A2) - X)))
   print "maxerr1 =", maxerr, log(maxerr)/log(2.0), "maxnum=", maxnum
   maxerr = maxnum = 0
   for a in range(2^(pa-2),2^pa):
      for c in range(2^pc):
         Xmin = infinity
         Xmax = -infinity
         for b in range(2^pb):
            A = (a*2^(pb+pc)+b*2^pc+c)/2^(pa+pb+pc)
            X = 1/sqrt(A) - T1[a,b]/2^p1
            X = round(X*2^p2)/2^p2
            Xmin = min (Xmin, X)
            Xmax = max (Xmax, X)
            A = (a*2^(pb+pc)+b*2^pc+c+1)/2^(pa+pb+pc)
            X = 1/sqrt(A) - T1[a,b]/2^p1
            X = round(X*2^p2)/2^p2
            Xmin = min (Xmin, X)
            Xmax = max (Xmax, X)
         T2[a,c] = round((Xmin + Xmax)/2*2^p2)
         maxnum = max(maxnum, abs(T2[a,c]))
         # print "a=", a, "c=", c, "T2=", T2[a,c]
         for b in range(2^pb):
            A = (a*2^(pb+pc)+b*2^pc+c)/2^(pa+pb+pc)
            X = 1/sqrt(A)
            maxerr = max(maxerr, n(abs(X - (T1[a,b]/2^p1 + T2[a,c]/2^p2))))
            A = (a*2^(pb+pc)+b*2^pc+c+1)/2^(pa+pb+pc)
            X = 1/sqrt(A)
            maxerr = max(maxerr, n(abs(X - (T1[a,b]/2^p1 + T2[a,c]/2^p2))))
   print "maxerr2 =", maxerr, log(maxerr)/log(2.0), "maxnum=", maxnum
   return [T1[a,b] for a in range(2^(pa-2),2^pa) for b in range(2^pb)], \
          [T2[a,c] for a in range(2^(pa-2),2^pa) for c in range(2^pc)]
*/

static const unsigned short T1[] = {8160, 8098, 8037, 7977, 7919, 7861, 7805, 7751, 7697, 7644, 7593, 7542, 7493, 7445, 7397, 7350, 7304, 7260, 7215, 7172, 7129, 7088, 7047, 7006, 6966, 6927, 6889, 6851, 6814, 6778, 6742, 6706, 6671, 6637, 6603, 6570, 6537, 6505, 6473, 6442, 6411, 6381, 6351, 6321, 6292, 6263, 6235, 6206, 6179, 6152, 6125, 6098, 6072, 6046, 6020, 5995, 5970, 5946, 5921, 5897, 5874, 5850, 5827, 5804, 5781, 5759, 5737, 5715, 5693, 5672, 5651, 5630, 5609, 5589, 5569, 5549, 5529, 5509, 5490, 5471, 5452, 5433, 5415, 5396, 5378, 5360, 5342, 5324, 5307, 5290, 5273, 5256, 5239, 5222, 5206, 5189, 5173, 5157, 5141, 5125, 5110, 5094, 5079, 5064, 5049, 5034, 5019, 5004, 4990, 4975, 4961, 4947, 4933, 4919, 4905, 4892, 4878, 4865, 4851, 4838, 4825, 4812, 4799, 4786, 4773, 4761, 4748, 4736, 4724, 4711, 4699, 4687, 4675, 4663, 4652, 4640, 4628, 4617, 4605, 4594, 4583, 4572, 4561, 4550, 4539, 4528, 4517, 4506, 4496, 4485, 4475, 4464, 4454, 4444, 4434, 4423, 4413, 4403, 4394, 4384, 4374, 4364, 4355, 4345, 4335, 4326, 4317, 4307, 4298, 4289, 4280, 4271, 4262, 4253, 4244, 4235, 4226, 4217, 4208, 4200, 4191, 4183, 4174, 4166, 4157, 4149, 4141, 4132, 4124, 4116, 4108, 4100};

static const short T2[] = {420, 364, 308, 252, 196, 140, 84, 28, -31, -87, -141, -194, -248, -302, -356, -410, 307, 267, 226, 185, 145, 104, 63, 21, -24, -65, -105, -145, -185, -224, -263, -303, 237, 205, 174, 142, 111, 79, 48, 16, -16, -45, -75, -105, -136, -167, -198, -229, 187, 161, 136, 110, 85, 61, 36, 12, -15, -41, -67, -92, -117, -142, -167, -192, 159, 138, 117, 96, 76, 54, 33, 12, -10, -30, -51, -71, -92, -113, -134, -154, 130, 112, 95, 77, 60, 43, 26, 9, -10, -28, -46, -64, -81, -98, -116, -133, 115, 100, 85, 70, 54, 39, 23, 8, -7, -22, -37, -52, -66, -82, -96, -111, 99, 86, 73, 60, 46, 32, 19, 6, -8, -21, -34, -47, -60, -73, -86, -100, 86, 75, 63, 52, 40, 28, 17, 5, -7, -19, -31, -43, -54, -66, -78, -90, 77, 66, 56, 46, 35, 25, 16, 6, -5, -15, -26, -36, -47, -57, -67, -77, 70, 60, 51, 42, 33, 23, 14, 5, -5, -14, -24, -33, -43, -52, -62, -72, 65, 57, 49, 40, 31, 23, 14, 6, -3, -12, -20, -28, -37, -45, -54, -62};

/* For GMP_NUMB_BITS=32: return a (1+20)-bit approximation x0 of 2^36/sqrt(a0).
   For GMP_NUMB_BITS=64: return a (1+40)-bit approximation x0 of 2^72/sqrt(a0).
   Assume a0 >= 2^(GMP_NUMB_BITS-2), and GMP_NUMB_BITS = 32 or 64.
   Must ensure: x0 <= 2^36/sqrt(a0) for GMP_NUMB_BITS=32,
                x0 <= 2^72/sqrt(a0) for GMP_NUMB_BITS=64. */
static mp_limb_t
mpn_rsqrtrem1 (mp_limb_t a0)
{
  mp_limb_t a = a0 >> (GMP_NUMB_BITS - 4);
  mp_limb_t b = (a0 >> (GMP_NUMB_BITS - 8)) & 0xf;
  mp_limb_t c = (a0 >> (GMP_NUMB_BITS - 12)) & 0xf;
  mp_limb_t x0, t;


  MPFR_STAT_STATIC_ASSERT (GMP_NUMB_BITS == 32 || GMP_NUMB_BITS == 64);

  x0 = ((mp_limb_t) T1[(a-4)*16+b] << 4) + T2[(a-4)*16+c];
  /* now x0 is a 16-bit approximation, with maximal error 2^(-9.46):
     -2^(-9.46) <= x0/2^16 - 1/sqrt(a/2^4) <= 2^(-9.46) */

#if GMP_NUMB_BITS == 32
  x0 >>= 1; /* reduce approximation to 1+15 bits */
  /* In principle, a0>>10 can have up to 22 bits, and (x0^2)>>12 can have up to
     20 bits, thus the product can have up to 42 bits, but since a0*x0^2 is very
     near 2^62, thus (a0>>10)*(x0^2)>>12 is very near 2^40, and when we reduce
     it mod 2^32 and interpret as a signed number in [-2^31, 2^31-1], we get
     the correct remainder (a0>>10)*(x0^2)>>12 - 2^40 */
  t = (mp_limb_signed_t) ((a0 >> 10) * ((x0 * x0) >> 12)) >> 8;
  /* |t| < 6742843 <= 2^23 - 256 (by exhaustive search) */
  t = (mp_limb_signed_t) t >> 8; /* now |t| < 2^15, thus |x0*t| < 2^31 */
  x0 = (x0 << 5) - ((mp_limb_signed_t) (x0 * t) >> 20);

  /* by exhaustive search on all possible values of a0, we get:
     -1.61 <= x0 - 2^36/sqrt(a0) <= 3.11 thus
     -2^(-19.3) < -1.54e-6 <= x0/2^20 - 2^16/sqrt(a0) <= 2.97e-6 < 2^(-18.3)
  */

  return x0 - 4; /* ensures x0 <= 2^36/sqrt(a0) */
#else /* GMP_NUMB_BITS = 64 */
  {
    mp_limb_t a1 = a0 >> (GMP_NUMB_BITS - 32);
    /* a1 has 32 bits, thus a1*x0^2 has 64 bits */
    t = (mp_limb_signed_t) (-a1 * (x0 * x0)) >> 32;
    /* t has 32 bits now */
  }
  x0 = (x0 << 15) + ((mp_limb_signed_t) (x0 * t) >> (17+1));
  /* now x0 is a 31-bit approximation (32 bits, 1 <= x0/2^31 <= 2),
     with maximal error 2^(-19.19) */

  {
    mp_limb_t a1 = a0 >> (GMP_NUMB_BITS - 40); /* a1 has 40 bits */
    t = (x0 * x0) >> 22; /* t has 40 bits */
    /* a1 * t has 80 bits, but we know the upper 19 bits cancel with 1 */
    t = (mp_limb_signed_t) (-a1 * t) >> 31;
    /* it remains 49 bits in theory, but t has only 31 bits at most */
  }
  x0 = (x0 << 9) + ((mp_limb_signed_t) (x0 * t) >> (31+1+49-40));

  /* now x0 is a 1+40-bit approximation,
     more precisely we have (experimentally):
     -2^(-38.2) < -3.16e-12 <= x0/2^40 - 2^32/sqrt(a0) <= 3.84e-12 < 2^(-37.9)
  */
  return x0 - 5; /* ensures x0 <= 2^72/sqrt(a0) */
#endif
}

/* Given as input np[0] and np[1], with B/4 <= np[1] (where B = 2^GMP_NUMB_BITS),
   mpn_sqrtrem2 returns a value x, 0 <= x <= 1, and stores values s in sp[0] and
   r in rp[0] such that:

   n := np[1]*B + np[0] = s^2 + x*B + r, with n < (s+1)^2

   or equivalently x*B + r <= 2*s.

   This comment is for GMP_NUMB_BITS=64 for simplicity, but the code is valid
   for any even value of GMP_NUMB_BITS.
   The algorithm used is the following, and uses Karp-Markstein's trick:
   - start from x, a 41-bit approximation of 2^72/sqrt(n1), with x <= 2^72/sqrt(n1)
   - y = floor(n1*x/2^72), which is a 32-bit approximation of sqrt(n1)
   - t = n1 - y^2
   - u = (x * t) >> 33
   - y = (y << 32) + u

   Proof:
   * we know that Newton's iteration for the reciprocal square root,

                x' = x + (x/2) (1 - a*x^2),                         (1)

     if evaluated with infinite precision, always produces x' <= 1/sqrt(a).
     See for example Lemma 3.14 in "Modern Computer Arithmetic" by Brent and
     Zimmermann.

   * if we multiply both sides of equation (1) by a, and write y0 = a*x
     and y0' = a*x', we get:

                y0' = y0 + (x/2)*(a - y0^2)                         (2)

     and since x' <= 1/sqrt(a), y0' <= sqrt(a).

   * now assume we replace y0 in (2) by y = y0 - e for some e >= 0.
     Equation (2) becomes:

                y' = y0 + (x/2)*(a - y0^2) - e*(1-x*y0+x*e/2)       (3)

     Since y0 = a*x, the term 1-x*y0+x*e/2 equals 1-a*x^2+x*e/2,
     which is non-negative since x <= 1/sqrt(a). Thus y' <= y0' <= sqrt(a).

     In practice, we should ensure y <= n1*x/2^64, so that t and u are
     non-negative, so that all right shifts round towards -infinity.
*/
static mp_limb_t
mpn_sqrtrem2 (mpfr_limb_ptr sp, mpfr_limb_ptr rp, mpfr_limb_srcptr np)
{
  mp_limb_t x, y, t, high, low;

  x = mpn_rsqrtrem1 (np[1]);
  /* we must have x^2*n1 <= 2^72 for GMP_NUMB_BITS=32
                         <= 2^144 for GMP_NUMB_BITS=64 */

#if GMP_NUMB_BITS == 32
  MPFR_ASSERTD ((double) x * (double) x * (double) np[1]
                < 4722366482869645213696.0);
  /* x is an approximation of 2^36/sqrt(n1), x has 1+20 bits */

  /* We know x/2^20 <= 2^16/sqrt(n1) + 2^(-18)
     thus n1*x/2^36 <= sqrt(n1) + 2^(-18)*n1/2^16 <= sqrt(n1) + 2^(-2). */

  /* compute y = floor(np[1]*x/2^36), cutting the upper 24 bits of n1 in two
     parts of 12 and 11 bits, which can be multiplied by x without overflow
     (warning: if we take 12 bits from low, it might overflow with x */
  high = np[1] >> 20; /* upper 12 bits from n1 */
  MPFR_ASSERTD((double) high * (double) x < 4294967296.0);
  low = (np[1] >> 9) & 0x7ff; /* next 11 bits */
  MPFR_ASSERTD((double) low * (double) x < 4294967296.0);
  y = high * x + ((low * x) >> 11); /* y approximates n1*x/2^20 */
  y = (y - 0x4000) >> 16; /* the constant 0x4000 = 2^(36-2-20) takes
                             into account the 2^(-2) error above, to ensure
                             y <= sqrt(n1) */
#else /* GMP_NUMB_BITS = 64 */
  MPFR_ASSERTD ((double) x * (double) x * (double) np[1]
                < 2.2300745198530623e43);
  /* x is an approximation of 2^72/sqrt(n1), x has 1+40 bits */

  /* We know x/2^40 <= 2^32/sqrt(n1) + 3.9e-12 <= 2^32/sqrt(n1) + 2^(-37)
     thus n1*x/2^72 <= sqrt(n1) + 2^(-37)*n1/2^32 <= sqrt(n1) + 2^(-5). */

  /* compute y = floor(np[1]*x/2^72), cutting the upper 48 bits of n1 in two
     parts of 24 and 23 bits, which can be multiplied by x without overflow
     (warning: if we take 24 bits from low, it might overflow with x) */
  high = np[1] >> 40; /* upper 24 bits from n1 */
  MPFR_ASSERTD((double) high * (double) x < 18446744073709551616.0);
  low = (np[1] >> 17) & 0x7fffff; /* next 23 bits */
  MPFR_ASSERTD((double) low * (double) x < 18446744073709551616.0);
  y = high * x + ((low * x) >> 23); /* y approximates n1*x/2^40 */
  y = (y - 0x8000000) >> 32; /* the constant 0x8000000 = 2^(72-5-40) takes
                                into account the 2^(-5) error above, to ensure
                                y <= sqrt(n1) */
#endif
  /* y is an approximation of sqrt(n1), with y <= sqrt(n1) */

  t = np[1] - y * y;
  MPFR_ASSERTD((mp_limb_signed_t) t >= 0);

#if GMP_NUMB_BITS == 32
  /* we now compute t = (x * t) >> (20 + 1), but x * t might have
     more than 32 bits thus we cut t in two parts of 12 and 11 bits */
  high = t >> 11;
  low = t & 0x7ff;
  MPFR_ASSERTD((double) high * (double) x < 4294967296.0);
  MPFR_ASSERTD((double) low * (double) x < 4294967296.0);
  t = high * x + ((low * x) >> 11); /* approximates t*x/2^11 */

  y = (y << (GMP_NUMB_BITS / 2)) + (t >> 10);
#else
  /* we now compute t = (x * t) >> (40 + 1), but x * t might have
     more than 64 bits thus we cut t in two parts of 24 and 23 bits */
  high = t >> 23;
  low = t & 0x7fffff;
  MPFR_ASSERTD((double) high * (double) x < 18446744073709551616.0);
  MPFR_ASSERTD((double) low * (double) x < 18446744073709551616.0);
  t = high * x + ((low * x) >> 23); /* approximates t*x/2^23 */

  y = (y << (GMP_NUMB_BITS / 2)) + (t >> 18);
#endif

  /* the correction code below assumes y >= 2^(GMP_NUMB_BITS - 1) */
  if (y < (MPFR_LIMB_ONE << (GMP_NUMB_BITS - 1)))
    y = MPFR_LIMB_ONE << (GMP_NUMB_BITS - 1);

  umul_ppmm (x, t, y, y);
  MPFR_ASSERTD(x < np[1] || (x == np[1] && t <= np[0])); /* y should not be too large */
  sub_ddmmss (x, t, np[1], np[0], x, t);

  /* Remainder x*2^GMP_NUMB_BITS+t should be <= 2*y (which implies x <= 1).
     If x = 1, it suffices to check t > 2*y mod 2^GMP_NUMB_BITS. */
  while (x > 1 || (x == 1 && t > 2 * y))
    {
      /* (y+1)^2 = y^2 + 2*y + 1 */
      x -= 1 + (t < (2 * y + 1));
      t -= 2 * y + 1;
      y ++;
      /* GMP_NUMB_BITS=32: average number of loops observed is 0.982,
         max is 3 (for n1=1277869843, n0=3530774227);
         GMP_NUMB_BITS=64: average number of loops observed is 0.593,
         max is 3 (for n1=4651405965214438987, n0=18443926066120424952).
      */
    }

  sp[0] = y;
  rp[0] = t;
  return x;
}

/* Special code for prec(r), prec(u) < GMP_NUMB_BITS. We cannot have
   prec(u) = GMP_NUMB_BITS here, since when the exponent of u is odd,
   we need to shift u by one bit to the right without losing any bit. */
static int
mpfr_sqrt1 (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t p = MPFR_GET_PREC(r);
  mpfr_prec_t exp_u = MPFR_EXP(u), exp_r, sh = GMP_NUMB_BITS - p;
  mp_limb_t u0, r0, rb, sb, mask, sp[2];
  mpfr_limb_ptr rp = MPFR_MANT(r);

  /* first make the exponent even */
  u0 = MPFR_MANT(u)[0];
  if (((unsigned int) exp_u & 1) != 0)
    {
      u0 >>= 1;
      exp_u ++;
    }
  MPFR_ASSERTD (((unsigned int) exp_u & 1) == 0);
  exp_r = exp_u / 2;

  /* then compute the integer square root of u0*2^GMP_NUMB_BITS */
  sp[0] = 0;
  sp[1] = u0;
  sb |= mpn_sqrtrem2 (&r0, &sb, sp);

  rb = r0 & (MPFR_LIMB_ONE << (sh - 1));
  mask = MPFR_LIMB_MASK(sh);
  sb |= (r0 & mask) ^ rb;
  rp[0] = r0 & ~mask;

  /* rounding */

  /* Note: if 1 and 2 are in [emin,emax], no overflow nor underflow
     is possible */
  if (MPFR_UNLIKELY (exp_r > __gmpfr_emax))
    return mpfr_overflow (r, rnd_mode, 1);

  /* See comments in mpfr_div_1 */
  if (MPFR_UNLIKELY (exp_r < __gmpfr_emin))
    {
      if (rnd_mode == MPFR_RNDN)
        {
          if ((exp_r == __gmpfr_emin - 1) && (rp[0] == ~mask) && rb)
            goto rounding; /* no underflow */
          if (exp_r < __gmpfr_emin - 1 || (rp[0] == MPFR_LIMB_HIGHBIT && sb == 0))
            rnd_mode = MPFR_RNDZ;
        }
      else if (!MPFR_IS_LIKE_RNDZ(rnd_mode, 0))
        {
          if ((exp_r == __gmpfr_emin - 1) && (rp[0] == ~mask) && (rb | sb))
            goto rounding; /* no underflow */
        }
      return mpfr_underflow (r, rnd_mode, 1);
    }

 rounding:
  MPFR_EXP (r) = exp_r;
  if (rb == 0 && sb == 0)
    {
      MPFR_ASSERTD(exp_r >= __gmpfr_emin);
      MPFR_ASSERTD(exp_r <= __gmpfr_emax);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (rb && sb == 0 &&
                      (rp[0] & (MPFR_LIMB_ONE << sh)) == 0))
        goto truncate;
      else
        goto add_one_ulp;
    }
  else if (MPFR_IS_LIKE_RNDZ(rnd_mode, 0))
    {
    truncate:
      MPFR_ASSERTD(exp_r >= __gmpfr_emin);
      MPFR_ASSERTD(exp_r <= __gmpfr_emax);
      MPFR_RET(-1);
    }
  else /* round away from zero */
    {
    add_one_ulp:
      rp[0] += MPFR_LIMB_ONE << sh;
      if (rp[0] == 0)
        {
          rp[0] = MPFR_LIMB_HIGHBIT;
          if (MPFR_UNLIKELY(exp_r + 1 > __gmpfr_emax))
            return mpfr_overflow (r, rnd_mode, 1);
          MPFR_ASSERTD(exp_r + 1 <= __gmpfr_emax);
          MPFR_ASSERTD(exp_r + 1 >= __gmpfr_emin);
          MPFR_SET_EXP (r, exp_r + 1);
        }
      MPFR_RET(1);
    }
}

#if GMP_NUMB_BITS == 64
/* For GMP_NUMB_BITS=64: return a (1+80)-bit approximation x = xp[1]*B+xp[0]
   of 2^144/sqrt(ap[1]*B+ap[0]).
   Assume ap[1] >= B/4, thus sqrt(ap[1]*B+ap[0]) >= B/2, thus x <= 2^81. */
static void
mpn_rsqrtrem2 (mpfr_limb_ptr xp, mpfr_limb_srcptr ap)
{
  mp_limb_t t1, t0, u1, u0, r0, r1, r2;

  MPFR_STAT_STATIC_ASSERT (GMP_NUMB_BITS == 64);

  xp[1] = mpn_rsqrtrem1 (ap[1]);

  /* now we should compute x + (x/2) * (1 - a*x^2), where the upper ~40 bits of
     a*x^2 should cancel with 1, thus we only need the following 40 bits */

  /* xp[1] has 1+40 bits, with xp[1] <= 2^72/sqrt(ap[1]) */
  umul_ppmm (t1, t0, xp[1], xp[1]);

  /* now t1 has at most 18 bits, with least 16 bits being its fractional value */

  u1 = ap[1] >> 48;
  u0 = (ap[1] << 16) | (ap[0] >> 48);

  /* u1 has the 16 most significant bits of a, and u0 the next 64 bits */

  /* we want t1*u1 << 48 + (t1*u0+t0*u1) >> 16 + t0*u0 >> 80
              [32 bits]       [64 bits]            [48 bits]
     but since we know the upper ~40 bits cancel with 1, we can ignore t1*u1. */
  umul_ppmm (r2, r1, t0, u0);
  r0 = t1 * u0; /* we ignore the upper 16 bits of the product */
  r1 = t0 * u1; /* we ignore the upper 16 bits of the product */
  r0 = r0 + r1 + r2;

  /* the upper ~8 bits of r0 should cancel with 1, we are interested in the next
     40 bits */

  umul_ppmm (t1, t0, xp[1], -r0);

  /* we should now add t1 >> 33 at xp[1] */
  xp[1] += t1 >> 33;
  xp[0] = t1 << 31;
  /* shift by 24 bits to the right, since xp[1] has 24 leading zeros,
     and we now expect 48 */
  xp[0] = (xp[1] << 40) | (xp[0] >> 24);
  xp[1] = xp[1] >> 24;
}

#if 0
/* this code seems to be slightly slower than the one below, for precision
   113 bits on a 64-bit machine, and in addition it requires GMP internals
   (for __gmpn_invert_limb) */
static mp_limb_t
mpn_sqrtrem4 (mpfr_limb_ptr sp, mpfr_limb_ptr rp, mpfr_limb_srcptr ap)
{
  mp_limb_t x, r, t;

  x = mpn_sqrtrem2 (sp + 1, rp + 1, ap + 2);

  /* now a1 = s1^2 + x*B + r1, with a1 = {ap+2,2}, s1 = sp[1], r1 = rp[1],
     with 0 <= x*B + r1 <= 2*s1 */

  /* use Newton's iteration for the square root, x' = 1/2 * (x + a/x), which
     rewrites as:

     x' = x + (a - x^2)/(2x)

     Since x' - sqrt(a) = (x-sqrt(a))^2/(2x), we have x' >= sqrt(a).
     If we round (a - x^2)/(2x) downwards, then we will still get
     an upper approximation of floor(sqrt(a)). */

  t = rp[1] << (GMP_NUMB_BITS - 1);
  rp[1] = (x << (GMP_NUMB_BITS - 1)) | (rp[1] >> 1); /* rp[1] = floor((a-x^2)/2) */
  /* Since x*B + r1 <= 2*s1, we now have rp[1] <= s1, and since __udiv_qrnnd_preinv
     requires its 3rd argument to be smaller than its 5th argument, we must
     distinguish the case rp[1] == s1. */
  if (MPFR_UNLIKELY(rp[1] == sp[1]))
    {
      /* Necessarily t=0 since x*B + r1 <= 2*s1 initially.
         Taking sp[0] = 2^GMP_NUMB_BITS (if representable) would be too large
         since it would mean sp[1] = sp[1]+1, and the remainder rp[1] would
         become negative. */
      sp[0] = MPFR_LIMB_MAX;
      r = sp[1];
    }
  else
    {
      r = __gmpn_invert_limb (sp[1]);
      __udiv_qrnnd_preinv (sp[0], r, rp[1], t, sp[1], r);
    }

  /* now sp[0] = floor((a-x^2)/(2x)) = floor((x*B+r1)/2/s1) */

  /* r1*B = 2*s1*s0 + 2*r
     a - s^2 = a1*B^2 + a0 - (s1*B+s0)^2
     = (s1^2+r1)*B^2 + a0 - s1^2*B^2 - 2*s1*s0*B - s0^2
     = r1*B^2 + a0 - 2*s1*s0*B - s0^2
     = 2*r*B + a0 - s0^2 */
  umul_ppmm (rp[1], rp[0], sp[0], sp[0]); /* s0^2 */
  t = -mpn_sub_n (rp, ap, rp, 2); /* a0 - s0^2 */
  t += 2 * (r >> (GMP_NUMB_BITS - 1));
  r <<= 1;
  rp[1] += r;
  t += rp[1] < r;
  if ((mp_limb_signed_t) t < 0)
    {
      mpn_sub_1 (sp, sp, 2, 1);
      t += mpn_addmul_1 (rp, sp, 2, 2);
      t += mpn_add_1 (rp, rp, 2, 1);
    }
  return t;
}
#else
/* Given as input ap[0-3], with B/4 <= ap[3] (where B = 2^GMP_NUMB_BITS),
   mpn_sqrtrem4 returns a value x, 0 <= x <= 1, and stores values s in sp[0-1] and
   r in rp[0-1] such that:

   n := ap[3]*B^3 + ap[2]*B^2 + ap[1]*B + ap[0] = s^2 + x*B + r, with n < (s+1)^2

   or equivalently x*B + r <= 2*s.

   This code currently assumes GMP_NUMB_BITS = 64, and takes on average
   135.68 cycles on an Intel i5-6500 for precision 113 bits */
static mp_limb_t
mpn_sqrtrem4 (mpfr_limb_ptr sp, mpfr_limb_ptr rp, mpfr_limb_srcptr ap)
{
  mp_limb_t x[2], t1, t0, r2, r1, h, l, u2, u1, b[4];

  MPFR_STAT_STATIC_ASSERT(GMP_NUMB_BITS == 64);

  mpn_rsqrtrem2 (x, ap + 2);

  /* x[1]*B+x[0] is a 80-bit approximation of 2^144/sqrt(ap[3]*B+ap[2]),
     and should be smaller */

  /* first compute y0 = a*x with at least 80 bits of precision */

  t1 = ap[3] >> 48;
  t0 = (ap[3] << 16) | (ap[2] >> 48);

  /* now t1:t0 is a (16+64)-bit approximation of a,
     (x1*B+x0) * (t1*B+t0) = (x1*t1)*B^2 + (x1*t0+x0*t1)*B + x0*t0 */
  r2 = x[1] * t1; /* r2 has 32 bits */
  umul_ppmm (h, r1, x[1], t0);
  r2 += h;
  umul_ppmm (h, l, x[0], t1);
  r1 += l;
  r2 += h + (r1 < l);
  umul_ppmm (h, l, x[0], t0);
  r1 += h;
  r2 += (r1 < h);

  /* r2 has 32 bits, r1 has 64 bits, thus we have 96 bits in total, we put 64
     bits in r2 and 16 bits in r1 */
  r2 = (r2 << 32) | (r1 >> 32);
  r1 = (r1 << 32) >> 48;

  /* we consider y0 = r2*2^16 + r1, which has 80 bits, and should be smaller than
     2^16*sqrt(ap[3]*B+ap[2]) */

  /* Now compute y0 + (x/2)*(a - y0^2), which should give ~160 correct bits.
     Since a - y0^2 has its ~80 most significant bits that cancel, it remains
     only ~48 bits. */

  /* now r2:r1 represents y0, with r2 of 64 bits and r1 of 16 bits,
     and we compute y0^2, whose upper ~80 bits should cancel with a:
     y0^2 = r2^2*2^32 + 2*r2*r1*2^16 + r1^2. */
  t1 = r2 * r2; /* we can simply ignore the upper 64 bits of r2^2 */
  umul_ppmm (h, l, r2, r1);
  t0 = l << 49; /* takes into account the factor 2 in 2*r2*r1 */
  u1 = (r1 * r1) << 32; /* temporary variable */
  t0 += u1;
  t1 += ((h << 49) | (l >> 15)) + (t0 < u1); /* takes into account the factor 2 */

  /* now t1:t0 >> 32 equals y0^2 mod 2^96, since y0 has 160 bits, we should shift
     t1:t0 by 64 bits to the right */
  t0 = ap[2] - t1 - (t0 != 0); /* we round downwards to get a lower approximation
                                  of sqrt(a) at the end */

  /* now t0 equals ap[3]*B+ap[2] - ceil(y0^2/2^32) */

  umul_ppmm (u2, u1, x[1], t0);
  umul_ppmm (h,  l,  x[0], t0);
  u1 += h;
  u2 += (u1 < h);

  /* divide by 2 to take into account the factor 1/2 in (x/2)*(a - y0^2) */
  u1 = (u2 << 63) | (u1 >> 1);
  u2 = u2 >> 1;

  /* u2:u1 approximates (x/2)*(ap[3]*B+ap[2] - y0^2/2^32) / 2^64,
     and should be smaller */

  r1 <<= 48; /* put back the most significant bits of r1 in place */

  /* add u2:u1 >> 16 to y0 */
  sp[0] = r1 + ((u2 << 48) | (u1 >> 16));
  sp[1] = r2 + (u2 >> 16) + (sp[0] < r1);

  mpn_mul_n (b, sp, sp, 2);
  b[2] = ap[2] - b[2] - mpn_sub_n (rp, ap, b, 2);

  /* invariant: the remainder {ap, 4} - {sp, 2}^2 is b[2]*B^2 + {rp, 2} */

  t0 = mpn_lshift (b, sp, 2, 1);

  /* Warning: the initial {sp, 2} might be < 2^127, thus t0 might be 0. */

  /* invariant: 2*{sp,2} = t0*B + {b, 2} */

  /* While the remainder is greater than 2*s we should subtract 2*s+1 to the
     remainder, and add 1 to the square root. This loop seems to be executed
     at most twice. */
  while (b[2] > t0 || (b[2] == t0 &&
                       (rp[1] > b[1] || (rp[1] == b[1] && rp[0] > b[0]))))
    {
      /* subtract 2*s to b[2]*B^2 + {rp, 2} */
      b[2] -= t0 + mpn_sub_n (rp, rp, b, 2);
      /* subtract 1 to b[2]*B^2 + {rp, 2}: b[2] -= mpn_sub_1 (rp, rp, 2, 1) */
      if (rp[0]-- == 0)
        b[2] -= (rp[1]-- == 0);
      /* add 1 to s */
      mpn_add_1 (sp, sp, 2, 1);
      /* add 2 to t0*B + {b, 2}: t0 += mpn_add_1 (b, b, 2, 2) */
      b[0] += 2;
      if (b[0] == 0)
        t0 += (++b[1] == 0);
    }

  return b[2];
}
#endif

/* Special code for GMP_NUMB_BITS < prec(r) < 2*GMP_NUMB_BITS,
   and GMP_NUMB_BITS < prec(u) <= 2*GMP_NUMB_BITS.
   This code should work for any value of GMP_NUMB_BITS, but since mpn_sqrtrem4
   currently assumes GMP_NUMB_BITS=64, it only works for GMP_NUMB_BITS=64. */
static int
mpfr_sqrt2 (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t p = MPFR_GET_PREC(r);
  mpfr_limb_ptr up = MPFR_MANT(u), rp = MPFR_MANT(r);
  mp_limb_t np[4], tp[2], rb, sb, mask;
  mpfr_prec_t exp_u = MPFR_EXP(u), exp_r, sh = 2 * GMP_NUMB_BITS - p;

  if (((unsigned int) exp_u & 1) != 0)
    {
      np[3] = up[1] >> 1;
      np[2] = (up[1] << (GMP_NUMB_BITS - 1)) | (up[0] >> 1);
      np[1] = up[0] << (GMP_NUMB_BITS - 1);
      exp_u ++;
    }
  else
    {
      np[3] = up[1];
      np[2] = up[0];
      np[1] = 0;
    }
  MPFR_ASSERTD (((unsigned int) exp_u & 1) == 0);
  exp_r = exp_u / 2;

  np[0] = 0;
  sb = mpn_sqrtrem4 (rp, tp, np);
  sb |= tp[0] | tp[1];
  rb = rp[0] & (MPFR_LIMB_ONE << (sh - 1));
  mask = MPFR_LIMB_MASK(sh);
  sb |= (rp[0] & mask) ^ rb;
  rp[0] = rp[0] & ~mask;

  /* rounding */
  if (MPFR_UNLIKELY (exp_r > __gmpfr_emax))
    return mpfr_overflow (r, rnd_mode, 1);

  /* See comments in mpfr_div_1 */
  if (MPFR_UNLIKELY (exp_r < __gmpfr_emin))
    {
      if (rnd_mode == MPFR_RNDN)
        {
          if (exp_r == __gmpfr_emin - 1 && (rp[1] == MPFR_LIMB_MAX &&
                                            rp[0] == ~mask) && rb)
            goto rounding; /* no underflow */
          if (exp_r < __gmpfr_emin - 1 || (rp[1] == MPFR_LIMB_HIGHBIT &&
                                           rp[0] == MPFR_LIMB_ZERO && sb == 0))
            rnd_mode = MPFR_RNDZ;
        }
      else if (!MPFR_IS_LIKE_RNDZ(rnd_mode, 0))
        {
          if (exp_r == __gmpfr_emin - 1 && (rp[1] == MPFR_LIMB_MAX &&
                                            rp[0] == ~mask) && (rb | sb))
            goto rounding; /* no underflow */
        }
      return mpfr_underflow (r, rnd_mode, 1);
    }

 rounding:
  MPFR_EXP (r) = exp_r;
  if (rb == 0 && sb == 0)
    {
      MPFR_ASSERTD(exp_r >= __gmpfr_emin);
      MPFR_ASSERTD(exp_r <= __gmpfr_emax);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (rb && sb == 0 &&
                      (rp[0] & (MPFR_LIMB_ONE << sh)) == 0))
        goto truncate;
      else
        goto add_one_ulp;
    }
  else if (MPFR_IS_LIKE_RNDZ(rnd_mode, 0))
    {
    truncate:
      MPFR_ASSERTD(exp_r >= __gmpfr_emin);
      MPFR_ASSERTD(exp_r <= __gmpfr_emax);
      MPFR_RET(-1);
    }
  else /* round away from zero */
    {
    add_one_ulp:
      rp[0] += MPFR_LIMB_ONE << sh;
      rp[1] += rp[0] == 0;
      if (rp[1] == 0)
        {
          rp[1] = MPFR_LIMB_HIGHBIT;
          if (MPFR_UNLIKELY(exp_r + 1 > __gmpfr_emax))
            return mpfr_overflow (r, rnd_mode, 1);
          MPFR_ASSERTD(exp_r + 1 <= __gmpfr_emax);
          MPFR_ASSERTD(exp_r + 1 >= __gmpfr_emin);
          MPFR_SET_EXP (r, exp_r + 1);
        }
      MPFR_RET(1);
    }
}
#endif /* GMP_NUMB_BITS == 64 */

#endif /* !defined(MPFR_GENERIC_ABI) && (GMP_NUMB_BITS == 32 || GMP_NUMB_BITS == 64) */

int
mpfr_sqrt (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mp_size_t rsize; /* number of limbs of r (plus 1 if exact limb multiple) */
  mp_size_t rrsize;
  mp_size_t usize; /* number of limbs of u */
  mp_size_t tsize; /* number of limbs of the sqrtrem remainder */
  mp_size_t k;
  mp_size_t l;
  mpfr_limb_ptr rp, rp0;
  mpfr_limb_ptr up;
  mpfr_limb_ptr sp;
  mp_limb_t sticky0; /* truncated part of input */
  mp_limb_t sticky1; /* truncated part of rp[0] */
  mp_limb_t sticky;
  int odd_exp;
  int sh; /* number of extra bits in rp[0] */
  int inexact; /* return ternary flag */
  mpfr_exp_t expr;
  MPFR_TMP_DECL(marker);

  MPFR_LOG_FUNC
    (("x[%Pu]=%.*Rg rnd=%d", mpfr_get_prec (u), mpfr_log_prec, u, rnd_mode),
     ("y[%Pu]=%.*Rg inexact=%d",
      mpfr_get_prec (r), mpfr_log_prec, r, inexact));

  if (MPFR_UNLIKELY(MPFR_IS_SINGULAR(u)))
    {
      if (MPFR_IS_NAN(u))
        {
          MPFR_SET_NAN(r);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_ZERO(u))
        {
          /* 0+ or 0- */
          MPFR_SET_SAME_SIGN(r, u);
          MPFR_SET_ZERO(r);
          MPFR_RET(0); /* zero is exact */
        }
      else
        {
          MPFR_ASSERTD(MPFR_IS_INF(u));
          /* sqrt(-Inf) = NAN */
          if (MPFR_IS_NEG(u))
            {
              MPFR_SET_NAN(r);
              MPFR_RET_NAN;
            }
          MPFR_SET_POS(r);
          MPFR_SET_INF(r);
          MPFR_RET(0);
        }
    }
  if (MPFR_UNLIKELY(MPFR_IS_NEG(u)))
    {
      MPFR_SET_NAN(r);
      MPFR_RET_NAN;
    }
  MPFR_SET_POS(r);

  /* See the note at the beginning of this file about __GNUC__. */
#if !defined(MPFR_GENERIC_ABI) && defined(__GNUC__) && \
    (GMP_NUMB_BITS == 32 || GMP_NUMB_BITS == 64)
  if (MPFR_GET_PREC (r) < GMP_NUMB_BITS && MPFR_GET_PREC (u) < GMP_NUMB_BITS)
    return mpfr_sqrt1 (r, u, rnd_mode);
#endif

#if !defined(MPFR_GENERIC_ABI) && GMP_NUMB_BITS == 64
  if (GMP_NUMB_BITS < MPFR_GET_PREC (r) && MPFR_GET_PREC (r) < 2*GMP_NUMB_BITS
      && MPFR_LIMB_SIZE(u) == 2)
    return mpfr_sqrt2 (r, u, rnd_mode);
#endif

  MPFR_TMP_MARK (marker);
  MPFR_UNSIGNED_MINUS_MODULO (sh, MPFR_GET_PREC (r));
  if (sh == 0 && rnd_mode == MPFR_RNDN)
    sh = GMP_NUMB_BITS; /* ugly case */
  rsize = MPFR_LIMB_SIZE(r) + (sh == GMP_NUMB_BITS);
  /* rsize is the number of limbs of r + 1 if exact limb multiple and rounding
     to nearest, this is the number of wanted limbs for the square root */
  rrsize = rsize + rsize;
  usize = MPFR_LIMB_SIZE(u); /* number of limbs of u */
  rp0 = MPFR_MANT(r);
  rp = (sh < GMP_NUMB_BITS) ? rp0 : MPFR_TMP_LIMBS_ALLOC (rsize);
  up = MPFR_MANT(u);
  sticky0 = MPFR_LIMB_ZERO; /* truncated part of input */
  sticky1 = MPFR_LIMB_ZERO; /* truncated part of rp[0] */
  odd_exp = (unsigned int) MPFR_GET_EXP (u) & 1;
  inexact = -1; /* return ternary flag */

  sp = MPFR_TMP_LIMBS_ALLOC (rrsize);

  /* copy the most significant limbs of u to {sp, rrsize} */
  if (MPFR_LIKELY(usize <= rrsize)) /* in case r and u have the same precision,
                                       we have indeed rrsize = 2 * usize */
    {
      k = rrsize - usize;
      if (MPFR_LIKELY(k))
        MPN_ZERO (sp, k);
      if (odd_exp)
        {
          if (MPFR_LIKELY(k))
            sp[k - 1] = mpn_rshift (sp + k, up, usize, 1);
          else
            sticky0 = mpn_rshift (sp, up, usize, 1);
        }
      else
        MPN_COPY (sp + rrsize - usize, up, usize);
    }
  else /* usize > rrsize: truncate the input */
    {
      k = usize - rrsize;
      if (odd_exp)
        sticky0 = mpn_rshift (sp, up + k, rrsize, 1);
      else
        MPN_COPY (sp, up + k, rrsize);
      l = k;
      while (sticky0 == MPFR_LIMB_ZERO && l != 0)
        sticky0 = up[--l];
    }

  /* sticky0 is non-zero iff the truncated part of the input is non-zero */

  tsize = mpn_sqrtrem (rp, NULL, sp, rrsize);

  /* a return value of zero in mpn_sqrtrem indicates a perfect square */
  sticky = sticky0 || tsize != 0;

  /* truncate low bits of rp[0] */
  sticky1 = rp[0] & ((sh < GMP_NUMB_BITS) ? MPFR_LIMB_MASK(sh)
                     : MPFR_LIMB_MAX);
  rp[0] -= sticky1;

  sticky = sticky || sticky1;

  expr = (MPFR_GET_EXP(u) + odd_exp) / 2;  /* exact */

  if (rnd_mode == MPFR_RNDZ || rnd_mode == MPFR_RNDD || sticky == MPFR_LIMB_ZERO)
    {
      inexact = (sticky == MPFR_LIMB_ZERO) ? 0 : -1;
      goto truncate;
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      /* if sh < GMP_NUMB_BITS, the round bit is bit (sh-1) of sticky1
                  and the sticky bit is formed by the low sh-1 bits from
                  sticky1, together with the sqrtrem remainder and sticky0. */
      if (sh < GMP_NUMB_BITS)
        {
          if (sticky1 & (MPFR_LIMB_ONE << (sh - 1)))
            { /* round bit is set */
              if (sticky1 == (MPFR_LIMB_ONE << (sh - 1)) && tsize == 0
                  && sticky0 == 0)
                goto even_rule;
              else
                goto add_one_ulp;
            }
          else /* round bit is zero */
            goto truncate; /* with the default inexact=-1 */
        }
      else /* sh = GMP_NUMB_BITS: the round bit is the most significant bit
              of rp[0], and the remaining GMP_NUMB_BITS-1 bits contribute to
              the sticky bit */
        {
          if (sticky1 & MPFR_LIMB_HIGHBIT)
            { /* round bit is set */
              if (sticky1 == MPFR_LIMB_HIGHBIT && tsize == 0 && sticky0 == 0)
                goto even_rule;
              else
                goto add_one_ulp;
            }
          else /* round bit is zero */
            goto truncate; /* with the default inexact=-1 */
        }
    }
  else /* rnd_mode=GMP_RDNU, necessarily sticky <> 0, thus add 1 ulp */
    goto add_one_ulp;

 even_rule: /* has to set inexact */
  if (sh < GMP_NUMB_BITS)
    inexact = (rp[0] & (MPFR_LIMB_ONE << sh)) ? 1 : -1;
  else
    inexact = (rp[1] & MPFR_LIMB_ONE) ? 1 : -1;
  if (inexact == -1)
    goto truncate;
  /* else go through add_one_ulp */

 add_one_ulp:
  inexact = 1; /* always here */
  if (sh == GMP_NUMB_BITS)
    {
      rp ++;
      rsize --;
      sh = 0;
    }
  /* now rsize = MPFR_LIMB_SIZE(r) */
  if (mpn_add_1 (rp0, rp, rsize, MPFR_LIMB_ONE << sh))
    {
      expr ++;
      rp0[rsize - 1] = MPFR_LIMB_HIGHBIT;
    }
  goto end;

 truncate: /* inexact = 0 or -1 */
  if (sh == GMP_NUMB_BITS)
    MPN_COPY (rp0, rp + 1, rsize - 1);

 end:
  /* Do not use MPFR_SET_EXP because the range has not been checked yet. */
  MPFR_ASSERTN (expr >= MPFR_EMIN_MIN && expr <= MPFR_EMAX_MAX);
  MPFR_EXP (r) = expr;
  MPFR_TMP_FREE(marker);

  return mpfr_check_range (r, inexact, rnd_mode);
}
