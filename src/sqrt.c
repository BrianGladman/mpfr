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
   with T1,T2 = bipartite(4,4,4,16,16). Note: we would get a slightly smaller
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

static const mp_limb_t T1[] = {130566, 129565, 128587, 127631, 126696, 125781, 124886, 124009, 123151, 122311, 121487, 120680, 119888, 119112, 118351, 117604, 116871, 116152, 115446, 114753, 114072, 113402, 112745, 112099, 111464, 110839, 110225, 109621, 109027, 108442, 107867, 107301, 106743, 106194, 105654, 105122, 104597, 104081, 103572, 103070, 102576, 102089, 101608, 101134, 100667, 100207, 99752, 99304, 98861, 98425, 97994, 97569, 97149, 96735, 96326, 95922, 95523, 95129, 94740, 94356, 93976, 93601, 93230, 92864, 92502, 92144, 91790, 91441, 91095, 90753, 90415, 90081, 89750, 89423, 89100, 88780, 88463, 88150, 87840, 87534, 87230, 86930, 86633, 86339, 86048, 85759, 85474, 85191, 84912, 84635, 84360, 84088, 83819, 83553, 83289, 83027, 82768, 82512, 82257, 82005, 81756, 81508, 81263, 81020, 80780, 80541, 80304, 80070, 79837, 79607, 79379, 79152, 78928, 78705, 78484, 78265, 78048, 77833, 77619, 77408, 77197, 76989, 76782, 76577, 76374, 76172, 75972, 75773, 75576, 75381, 75187, 74994, 74803, 74613, 74425, 74239, 74053, 73869, 73687, 73505, 73325, 73147, 72969, 72793, 72619, 72445, 72273, 72102, 71932, 71763, 71596, 71429, 71264, 71100, 70937, 70776, 70615, 70455, 70297, 70139, 69983, 69828, 69673, 69520, 69368, 69216, 69066, 68917, 68768, 68621, 68475, 68329, 68184, 68041, 67898, 67756, 67615, 67475, 67336, 67197, 67060, 66923, 66787, 66652, 66518, 66384, 66252, 66120, 65989, 65858, 65729, 65600};

static const mp_limb_signed_t T2[] = {415, 360, 304, 249, 194, 139, 85, 29, -34, -89, -143, -198, -252, -307, -361, -415, 304, 264, 224, 183, 143, 102, 62, 21, -24, -64, -105, -145, -185, -225, -265, -305, 236, 204, 173, 142, 111, 79, 48, 17, -18, -50, -81, -112, -143, -174, -205, -236, 190, 164, 139, 114, 89, 64, 39, 13, -14, -39, -65, -90, -114, -140, -165, -189, 157, 136, 115, 94, 74, 53, 32, 11, -12, -33, -53, -74, -95, -116, -137, -157, 133, 115, 97, 80, 63, 45, 27, 9, -10, -27, -45, -62, -80, -97, -115, -133, 114, 99, 84, 69, 54, 38, 23, 8, -9, -24, -39, -53, -68, -84, -99, -114, 99, 86, 73, 60, 46, 33, 20, 7, -8, -21, -34, -47, -60, -73, -86, -99, 88, 76, 65, 53, 41, 30, 18, 6, -7, -18, -30, -41, -53, -64, -76, -87, 78, 68, 57, 47, 36, 26, 16, 5, -6, -16, -26, -37, -47, -58, -68, -78, 70, 61, 51, 42, 33, 24, 14, 5, -5, -15, -24, -33, -42, -52, -61, -70, 63, 55, 47, 38, 30, 21, 13, 4, -5, -13, -22, -30, -38, -47, -55, -63};

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

  x0 = T1[((a - 4) << 4) + b] + T2[((a - 4) << 4) + c];
  /* now x0 is a 16-bit approximation, with maximal error 2^(-9.48):
     -2^(-9.48) <= x0/2^16 - 1/sqrt(a0/2^64) <= 2^(-9.48)
     The worst case is obtained for floor(a0/2^52) = 1265. */

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
     with maximal error 2^(-19.19):
     -2^(-19.22) <= x0/2^31 - 1/sqrt(a0/2^64) <= 2^(-19.22).
     The worst case is attained for a1 = 1326448625, with x0 = 3864240837. */

  {
    mp_limb_t a1 = a0 >> (GMP_NUMB_BITS - 41); /* 2^39 <= a1 < 2^41 */
    t = (x0 * x0) >> 23; /* t < 2^41 */
    /* a1 * t has at most 82 bits, but we know the upper 19 bits cancel with 1 */
    t = (mp_limb_signed_t) (-a1 * t) >> 31;
    /* it remains 49 bits in theory, but t has only 31 bits at most */
    /* the error on t is at most 1 (from the truncation >> 31,
       plus 2^(41-31) from the truncation of a0 into a1 multiplied by t,
       and another 2^(41-31) from the truncation of x0^2 multiplied by a1,
       thus at most 2^11+1 (in fact since t/2^41 is about 2^41/a1, the maximal
       error is when a0 is near 2^64 or 2^62, and is 1.25*2^10 + 1) */
  }
  x0 = (x0 << 9) + ((mp_limb_signed_t) (x0 * t) >> (31+1+49-40));
  /* The error on x0 is at most (1.25*2^10+1)*2^(32-41) from the previous error on t,
     plus 1 from the truncation by >> (31+1+49-40), thus at most 2.5+2^(-9)+1.
     Since the mathematical error is bounded by 2^(-37.85) <= 4.44 ulps [but this
     error can only make x0 smaller], we have at the end:
     -(3.5+2^(-9))/2^40 - 2^(-37.85) <= x0/2^40 - 2^32/sqrt(a0) <= (3.5+2^(-9))/2^40:

     -2^(-37.01) < x0/2^40 - 2^32/sqrt(a0) < 2^(-38.19)
  */
  x0 = x0 - 4; /* ensures x0 <= 2^72/sqrt(a0) */
  /* x0-2 fails for example for a0=4967732205162787840 (gives x0=2118754313105) */
  /* now -2^(-36.30) < x0/2^40 - 2^32/sqrt(a0) < 0 */
  return x0;
#endif
}

#if GMP_NUMB_BITS == 64
/* The Taylor coefficient of order 0 of sqrt(i/2^8+x) is
   U[i-64][0]/2^64 + U[i-64][1]/2^128, then the Taylor coefficient of order j is
   (up to sign) U[i-64][j+1]/2^(64-8*j).
   The maximal number of bits is:
   j=1:64 j=2:56 j=3:49 j=4:43 j=5:36 j=6:30 j=7:23
   The sign is implicit: u[j] < 0 for j even except j=0.
   The maximal error is < .927e-21 (attained for i=64). */

#include "sqrt_tab.h"

/* Return an approximation of sqrt(2^64*n), with 2^62 <= n < 2^64,
   and error < 1 ulp (in unknown direction).
   We use a Taylor polynomial of degree 7. */
static mp_limb_t
mpn_sqrtrem2_approx (mp_limb_t n)
{
  int i = n >> 56;
  mp_limb_t x, h, l;
  const mp_limb_t *u;

  x = n << 8;
  u = U[i - 64];
  umul_ppmm (h, l, u[8], x);
  /* the truncation error on h is at most 1 here */
  umul_ppmm (h, l, u[7] - h, x);
  /* the truncation error on h is at most 2 */
  umul_ppmm (h, l, u[6] - h, x);
  /* the truncation error on h is at most 3 */
  umul_ppmm (h, l, u[5] - h, x);
  /* the truncation error on h is at most 4 */
  umul_ppmm (h, l, u[4] - h, x);
  /* the truncation error on h is at most 5 */
  umul_ppmm (h, l, u[3] - h, x);
  /* the truncation error on h is at most 6 */
  umul_ppmm (h, l, u[2] - h, x >> 8); /* here we shift by 8 since u[0] has the
                                         same weight 1/2^64 as u[2], the truncation
                                         error on h + l/2^64 is at most 6/2^8 */
  add_ssaaaa (h, l, h, l, u[0], u[1]);
  /* Since the above addition is exact, the truncation error on h + l/2^64
     is still at most 6/2^8. Together with the mathematical error < .927e-21*2^64,
     the total error on h + l/2^64 is < 0.0406 */
  return h + (l >> 63); /* round to nearest */
}

/* put in rh,rl the upper 2 limbs of the product xh,xl * yh,yl,
   with error less than 3 ulps */
#define umul_ppmm2(rh,rl,xh,xl,yh,yl)    \
  {                                      \
    mp_limb_t _h, _l;                    \
    umul_ppmm (rh, rl, xh, yh);          \
    umul_ppmm (_h, _l, xh, yl);          \
    rl += _h;                            \
    rh += (rl < _h);                     \
    umul_ppmm (_h, _l, xl, yh);          \
    rl += _h;                            \
    rh += (rl < _h);                     \
   }

/* Put in rp[1]*2^64+rp[0] an approximation of sqrt(2^128*n),
   with 2^126 <= n := np[1]*2^64 + np[0] < 2^128.
   We use a Taylor polynomial of degree 14.
   The coefficients of degree 0 to 9 are represented by two 64-bit limbs
   (most significant first), the remaining coefficients by one 64-bit limb,
   thus the degree-0 coefficient is u[0]/2^64 + u[1]/2^128,
   the degree-1 coefficient is u[2]/2^64 + u[3]/2^128,
   ...,
   the degree-9 coefficient is u[18]/2^64 + u[19]/2^128,
   the degree-10 coefficient is u[20]/2^64,
   ...,
   the degree-14 coefficient is u[24]/2^64. */
static void
mpn_sqrtrem4_approx (mpfr_limb_ptr rp, mpfr_limb_srcptr np)
{
  int i = np[1] >> 56;
  mp_limb_t xh, xl, h, l, yh, yl;
  const mp_limb_t *u;

  xh = (np[1] << 8) | (np[0] >> 56);
  xl = np[0] << 8;
  u = V[i - 64];
  umul_ppmm (h, l, u[24],     xh);
  umul_ppmm (h, l, u[23] - h, xh);
  umul_ppmm (h, l, u[22] - h, xh);
  umul_ppmm (h, l, u[21] - h, xh);
  umul_ppmm (h, l, u[20] - h, xh);
  /* now we have to deal with two limbs */
  sub_ddmmss (yh, yl, u[18], u[19], 0, h);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  sub_ddmmss (yh, yl, u[16], u[17], h, l);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  sub_ddmmss (yh, yl, u[14], u[15], h, l);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  sub_ddmmss (yh, yl, u[12], u[13], h, l);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  sub_ddmmss (yh, yl, u[10], u[11], h, l);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  sub_ddmmss (yh, yl, u[8], u[9], h, l);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  sub_ddmmss (yh, yl, u[6], u[7], h, l);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  sub_ddmmss (yh, yl, u[4], u[5], h, l);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  sub_ddmmss (yh, yl, u[2], u[3], h, l);
  umul_ppmm2 (h, l, yh, yl, xh, xl);
  add_ssaaaa (rp[1], rp[0], u[0], u[1], h, l);
}
#endif /* GMP_NUMB_BITS == 64 */

#if GMP_NUMB_BITS == 32
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
  /* x is an approximation of 2^72/sqrt(n1), x < 2^41 */

  /* We know 2^32/sqrt(n1) - 2^(-36.30) <= x/2^40 <= 2^32/sqrt(n1) */

  /* compute y = floor(np[1]*x/2^72): we cut the upper 48 bits of n1 in two parts of
     24 and 23 bits, which can be multiplied by x without overflow (warning: if we
     take 24 bits from low, it might overflow with x). We could simply write
     umul_ppmm (high, low, np[1], x) followed by y = high >> 8 but the following
     code is faster. */
  high = np[1] >> 40; /* upper 24 bits from n1 */
  MPFR_ASSERTD((double) high * (double) x < 18446744073709551616.0);
  low = (np[1] >> 17) & 0x7fffff; /* next 23 bits */
  MPFR_ASSERTD((double) low * (double) x < 18446744073709551616.0);
  y = high * x + ((low * x) >> 23); /* y approximates n1*x/2^40 */
  y = y >> 32;
#endif

  /* Now y is an approximation of sqrt(n1), with y <= sqrt(n1). The errors are:
     (a) mathematical error: according to Lemma 3.14 from "Modern Computer
         Arithmetic", it is bounded by 3/2*x^3/theta^4*2^(-2*36.30), where
         x <= theta + 2^(-36.30), thus x/theta <= 1 + 2^(-36.30), thus the
         mathematical error is bounded by 2^(-72.01), thus 2^(-40.01) for a
         value normalized by 2^32;
     (b) truncation error in Newton's floor(np[1]*x/2^72): it is bounded by
         2^17*x/2^72 for the neglected part of np[1], and 1 in (low * x) >> 23.
     Thus we have:
     sqrt(np[1]) - 1 - 2^(-13.99) < y <= sqrt(np[1]) */

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
  if (MPFR_UNLIKELY(y < MPFR_LIMB_HIGHBIT))
    y = MPFR_LIMB_HIGHBIT;

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
#endif

/* Special code for prec(r), prec(u) < GMP_NUMB_BITS. We cannot have
   prec(u) = GMP_NUMB_BITS here, since when the exponent of u is odd,
   we need to shift u by one bit to the right without losing any bit. */
static int
mpfr_sqrt1 (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t p = MPFR_GET_PREC(r);
  mpfr_prec_t exp_u = MPFR_EXP(u), exp_r, sh = GMP_NUMB_BITS - p;
  mp_limb_t u0, r0, rb, sb, mask = MPFR_LIMB_MASK(sh);
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
#if GMP_NUMB_BITS == 64
  r0 = mpn_sqrtrem2_approx (u0);
  sb = 1; /* when we can round correctly with the approximation, the sticky bit
             is non-zero */

  /* Since the exact square root is in [r0 - 0.5406, r0 + 0.5406], we can round
     correctly except when the last sh-1 bits of r0 are 000...000. */
  if (MPFR_UNLIKELY((r0 & (mask >> 1)) == 0))
    {
      umul_ppmm (rb, sb, r0, r0);
      /* for the exact square root, we should have 0 <= (u0-rb)*2^64 - sb <= 2*r0 */
      if (rb > u0 || (rb == u0 && sb > 0)) /* r0 is too large */
        {
          r0 --;
          umul_ppmm (rb, sb, r0, r0);
        }
      /* if u0 <= rb + 1, then (u0-rb)*2^64 - sb <= 2^64 <= 2*r0
         if u0 >= rb + 3, then (u0-rb)*2^64 - sb > 2*2*64 > 2*r0 */
      else if (u0 > rb + 2 || (u0 == rb + 2 && -sb > 2 * r0))
        {
          r0 ++;
          umul_ppmm (rb, sb, r0, r0);
        }
      sub_ddmmss (rb, sb, u0, 0, rb, sb);
      /* now we should have rb*2^64 + sb <= 2*r0 */
      MPFR_ASSERTN(rb == 0 || (rb == 1 && sb <= 2 * r0));
      sb = rb | sb;
    }
#else /* 32-bit variant. FIXME: write mpn_sqrtrem2_approx for GMP_NUMB_BITS=32,
         then use the above code and get rid of mpn_sqrtrem2 */
  {
    mp_limb_t sp[2];
    sp[1] = u0;
    sp[0] = 0;
    sb |= mpn_sqrtrem2 (&r0, &sb, sp);
  }
#endif

  rb = r0 & (MPFR_LIMB_ONE << (sh - 1));
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

#if 1
  np[0] = 0;
  sb = mpn_sqrtrem4 (rp, tp, np);
  sb |= tp[0] | tp[1];
  rb = rp[0] & (MPFR_LIMB_ONE << (sh - 1));
  mask = MPFR_LIMB_MASK(sh);
  sb |= (rp[0] & mask) ^ rb;
  rp[0] = rp[0] & ~mask;
#else
  mpn_sqrtrem4_approx (rp, np);
#endif

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
