/* random_deviate routines for mpfr_erandom and mpfr_nrandom.
   If MPFR_ALT_RANDOM is defined, include also
   mpfr_urandom_alt (rop, state, rnd_mode) -- mimic the functionality of
   mpfr_urandom.

Copyright 2013 Free Software Foundation, Inc.
Contributed by Charles Karney <charles@karney.com>, SRI International.

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
 * A mpfr_random_deviate represents the initial portion e bits of a random
 * deviate uniformly distributed in (0,1) as
 *
 *  typedef struct {
 *    unsigned long e;            // bits in the fraction 
 *    unsigned long h;            // the high W bits of the fraction
 *    mpz_t f;                    // the rest of the fraction
 *  } mpfr_random_deviate_t[1];
 *
 * e is always a multiple of RANDOM_CHUNK.  The first RANDOM_CHUNK bits, the
 * high fraction, are held in an unsigned long, h, and the rest are held in an
 * mpz_t, f.  The data in h is undefined if e == 0 and, similarly the data in f
 * is undefined if e <= RANDOM_CHUNK.
 */

#include <random_deviate.h>

/*
 * RANDOM_CHUNK can be picked in the range 1 <= RANDOM_CHUNK <= 64.  Low values
 * of RANDOM_CHUNK are good for testing, since they are more likely to make
 * bugs obvious.  For portability, pick RANDOM_CHUNK <= 32 (since an unsigned
 * long may only hold 32 bits).  For reproducibility across platforms,
 * standardize on RANDOM_CHUNK = 32.
 *
 * When RANDOM_CHUNK = 32, this representation largely avoids manipulating
 * mpz's (until the final cast to an mpfr is done).  In addition
 * mpfr_random_deviate_less usually entails just a single comparison of
 * unsigned longs.  In this way, we can stick with the published interface for
 * extracting portions of an mpz (namely through mpz_tstbit) without hurting
 * efficiency.
 */
#if !defined(RANDOM_CHUNK)
/* note: for MPFR, we could use RANDOM_CHUNK = 32 or 64 according to the
   number of bits per limb, but we use 32 everywhere to get reproducible
   results on 32-bit and 64-bit computers */
#define RANDOM_CHUNK 32     /* Require 1 <= RANDOM_CHUNK <= 32; recommend 32 */
#endif

/*
 * If RANDOM_SIMPLE is 1, then the high fraction is not used allowing for
 * considerably simpler implementations.  These should give precisely the same
 * results as for non-RANDOM_SIMPLE and act as a cross check.  When requesting
 * results with a low precision the non-RANDOM_SIMPLE implementation is about
 * 1.8 times faster.
 */
#if !defined(RANDOM_SIMPLE)
#define RANDOM_SIMPLE 0         /* Use 0 except if testing */
#endif

static const unsigned W = RANDOM_CHUNK;

/* allocate and set to (0,1) */
void mpfr_random_deviate_init(mpfr_random_deviate_t x) {
  mpz_init(x->f); x->e = 0;
}

/* reset to (0,1) */
void mpfr_random_deviate_reset(mpfr_random_deviate_t x) {
#if RANDOM_SIMPLE
  mpz_set_ui(x->f, 0);
#endif
  x->e = 0;
}

/* deallocate */
void mpfr_random_deviate_clear(mpfr_random_deviate_t x) { mpz_clear(x->f); }

/* swap two random deviates */
void mpfr_random_deviate_swap(mpfr_random_deviate_t x,
                              mpfr_random_deviate_t y) {
  unsigned long t;
  t = x->e; x->e = y->e; y->e = t;
#if !RANDOM_SIMPLE
  t = x->h; x->h = y->h; y->h = t;
#endif
  mpz_swap(x->f, y->f);
}

/* ensure x has at least k bits */
static void random_deviate_generate(mpfr_random_deviate_t x, unsigned long k,
                                    gmp_randstate_t r, mpz_t t) {
  /* if t is non-null, it is used as a temporary */
#if RANDOM_SIMPLE
  if (t) {
    if (x->e >= k) return;
    /* force computation of first W bits separately so the results match the
     * non-simple case */
    random_deviate_generate(x, W, r, 0);
    /* passed a mpz_t so compute needed bits in one call to mpz_urandomb */
    k = ((k + (W-1)) / W) * W;  /* Round up to multiple of W */
    k -= x->e;                  /* The number of new bits */
    if (k == 0) return; /* Generation of first W bits might mean we can exit */
    mpz_urandomb(t, r, k);
    mpz_mul_2exp(x->f, x->f, k);
    mpz_add(x->f, x->f, t);
    x->e += k;
  } else {
    /* no mpz_t so compute the bits W at a time via gmp_urandomb_ui */
    while (x->e < k) {
      unsigned long w = gmp_urandomb_ui(r, W);
      mpz_mul_2exp(x->f, x->f, W);
      mpz_add_ui(x->f, x->f, w);
      x->e += W;
    }
  }
#else
  if (x->e >= k) return;
  if (x->e == 0) {
      x->h = gmp_urandomb_ui(r, W); /* Generate the high fraction */
      x->e = W;
      if (x->e >= k) return;    /* Maybe that's it? */
  }
  if (t) {
    /* passed a mpz_t so compute needed bits in one call to mpz_urandomb */
    k = ((k + (W-1)) / W) * W;  /* Round up to multiple of W */
    k -= x->e;                  /* The number of new bits */
    mpz_urandomb(x->e == W ? x->f : t, r, k); /* Copy directly to x->f? */
    if (x->e > W) {
      mpz_mul_2exp(x->f, x->f, k);
      mpz_add(x->f, x->f, t);
    }
    x->e += k;
  } else {
    /* no mpz_t so compute the bits W at a time via gmp_urandomb_ui */
    while (x->e < k) {
      unsigned long w = gmp_urandomb_ui(r, W);
      if (x->e == W)
        mpz_set_ui(x->f, w);
      else {
        mpz_mul_2exp(x->f, x->f, W);
        mpz_add_ui(x->f, x->f, w);
      }
      x->e += W;
    }
  }
#endif
}

/*
 * return index [-1..63] of highest bit set.  Return -1 if x = 0, 63 is if x =
 * ~0.  (From Algorithms for programmers by Joerg Arndt.)
 */
static int highest_bit_idx(unsigned long x) {
  int r = 0;
  if (x == 0) return -1;
  /* handle 64-bit unsigned longs in a way that doesn't trigger warnings when
   * they are only 32-bits */
  if (x & ~0xffffffffUL) { x >>= 16; x >>= 16; r +=32; }
  if (x &  0xffff0000UL) { x >>= 16; r += 16; }
  if (x &  0x0000ff00UL) { x >>=  8; r +=  8; }
  if (x &  0x000000f0UL) { x >>=  4; r +=  4; }
  if (x &  0x0000000cUL) { x >>=  2; r +=  2; }
  if (x &  0x00000002UL) {           r +=  1; }
  return r;
}

/* return position of leading bit, counting from 1 */
static long random_deviate_leading_bit(mpfr_random_deviate_t x,
                                       gmp_randstate_t r) {
#if !RANDOM_SIMPLE
  random_deviate_generate(x, W, r, 0);
  if (x->h) return W - highest_bit_idx(x->h);
  random_deviate_generate(x, 2 * W, r, 0);
#endif
  while (mpz_sgn(x->f) == 0)
    random_deviate_generate(x, x->e + W, r, 0);
  return x->e + 1 - mpz_sizeinbase(x->f, 2);
}

/* return kth bit of fraction, representing 2^-k */
int mpfr_random_deviate_tstbit(mpfr_random_deviate_t x, unsigned long k,
                               gmp_randstate_t r) {
  if (k == 0) return 0;
  random_deviate_generate(x, k, r, 0);
#if !RANDOM_SIMPLE
  if (k <= W) return (x->h >> (W - k)) & 1UL;
#endif
  return mpz_tstbit(x->f, x->e - k);
}

/* compare two random deviates, x < y */
int mpfr_random_deviate_less(mpfr_random_deviate_t x, mpfr_random_deviate_t y,
                             gmp_randstate_t r) {
  unsigned long k = 1;
  if (x == y) return 0;
#if !RANDOM_SIMPLE
  random_deviate_generate(x, W, r, 0);
  random_deviate_generate(y, W, r, 0);
  if (x->h != y->h) return x->h < y->h; /* Compare the high fractions */
  k += W;
#endif
  for (; ; ++k) {             /* Compare the rest of the fraction bit by bit */
    int a = mpfr_random_deviate_tstbit(x, k, r);
    int b = mpfr_random_deviate_tstbit(y, k, r);
    if (a != b) return a < b;
  }
}

/* set mpfr_t z = (neg ? -1 : 1) * (n + x) */
int mpfr_random_deviate_value(int neg, unsigned long n,
                              mpfr_random_deviate_t x, mpfr_t z,
                              gmp_randstate_t r, mpfr_rnd_t rnd) {
  /* r is used to add as many bits as necessary to match the precision of z */
  long l;                     /* The leading bit is 2^l */
  long p = mpfr_get_prec(z);  /* Number of bits in result */
  mpz_t t;
  int inex;

  if (n == 0)
    l = -random_deviate_leading_bit(x, r); /* l < 0 */
  else
    l = highest_bit_idx(n); /* l >= 0 */

  /*
   * Leading bit is 2^l; thus the trailing bit in result is 2^(l-p+1) =
   * 2^-(p-l-1).  For the sake of illustration, take l = 0 and p = 4, thus bits
   * through the 1/8 position need to be generated; assume that these bits are
   * 1.010 = 10/8 which represents a deviate in the range (10,11)/8.
   *
   * If the rounding mode is one of RNDZ, RNDU, RNDD, RNDA, we add a 1 bit to
   * the result to give 1.0101 = (10+1/2)/8.  When this is converted to a MPFR
   * the result is rounded to 10/8, 11/8, 10/8, 11/8, respectively, and the
   * inexact flag is set to -1, 1, -1, 1.
   *
   * If the rounding mode is RNDN, an additional random bit must be generated
   * to determine if the result is in (10,10+1/2)/8 or (10+1/2,11)/8.  Assume
   * that this random bit is 0, so the result is 1.0100 = (10+0/2)/8.  Then an
   * additional 1 bit is added to give 1.010101 = (10+1/4)/8.  This last bit
   * avoids the "round ties to even rule" (because there are no ties) and sets
   * the inexact flag so that the result is 10/8 with the inexact flag = 1.
   *
   * Here we always generate at least 2 additional random bits, so that bit
   * position 2^-(p-l+1) is generated.  (The result often contains more random
   * bits than this because random bits are added in batches of W and because
   * additional bits may have been required in the process of generating the
   * random deviate.)  The integer and all the bits in the fraction are then
   * copied into an mpz, the least significant bit is unconditionally set to 1,
   * the sign is set, and the result together with the exponent -x->e is used
   * to generate an mpfr using mpfr_set_z_2exp.
   *
   * If random bits were very expensive, we would only need to generate to the
   * 2^-(p-l-1) bit (no extra bits) for the RNDZ, RNDU, RNDD, RNDA modes and to
   * the 2^-(p-l) bit (1 extra bit) for RNDN.  By always generating 2 bits we
   * save on some bit shuffling when formed the mpz to be converted to an mpfr.
   * The implementation of the RandomNumber class in RandomLib illustrates the
   * more parsimonious approach (which was taken to allow accurate counts of
   * the number of random digits to be made).
   */
  mpz_init(t);
  /*
   * This is the only call to random_deviate_generate where a mpz_t is passed
   * (because an arbitrarily large number of bits may need to be generated).
   */
  if (p - l + 1 > 0)
    random_deviate_generate(x, p + 1 - l, r, t);
#if RANDOM_SIMPLE
  mpz_set_ui(t, n);             /* The integer part */
  mpz_mul_2exp(t, t, x->e);     /* Shift to allow for fraction */
  mpz_add(t, t, x->f);          /* Add high fraction */
#else
  if (n == 0) {
    /* Since the minimum prec is 2 we know that x->h has been generated. */
    mpz_set_ui(t, x->h);        /* Set high fraction */
  } else {
    mpz_set_ui(t, n);           /* The integer part */
    if (x->e > 0) {
      mpz_mul_2exp(t, t, W);    /* Shift to allow for high fraction */
      mpz_add_ui(t, t, x->h);   /* Add high fraction */
    }
  }
  if (x->e > W) {
    mpz_mul_2exp(t, t, x->e - W); /* Shift to allow for low fraction */
    mpz_add(t, t, x->f);          /* Add low fraction */
  }
#endif
  /*
   * We could trim off any excess bits here by shifting rightward.  This is an
   * unnecessary complication.
   */
  mpz_setbit(t, 0);      /* Set the trailing bit so result is always inexact */
  if (neg) mpz_neg(t, t);
  /*
   * Let mpfr_set_z_2exp do all the work of rounding to the requested
   * precision, setting overflow/underflow flags, and returning the right
   * inexact value.
   */
  inex = mpfr_set_z_2exp(z, t, -x->e, rnd);
  mpz_clear(t);
  return inex;
}

#if MPFR_ALT_RANDOM
/*
 * Sampling from the uniform distribution in [0,1].  This duplicates the
 * functionality of mpfr_urandom.  It is included just to exercise
 * mpfr_random_deviate_value.  The relative timings of mpfr_urandom_alt vs
 * mpfr_urandom shows that the overhead of generating the random bits in a
 * mpz_t (via mpfr_random_deviate) and then converting the result to an mpfr_t
 * is about 3/2.
 *
 * It would be possible to generalize this to sample in [m,n], by sampling an
 * integer uniformly in [m,n) and setting the sign and integer in the call to
 * mpfr_random_deviate_value.
 */

/* return a uniform random deviate in [0,1] as a MPFR */
int mpfr_urandom_alt(mpfr_t z, gmp_randstate_t r, mpfr_rnd_t rnd) {
  mpfr_random_deviate_t x;
  int inex;
  mpfr_random_deviate_init(x);
  /*
   * sampling in the integer range [m,n] is illustrated here assuming n > m and
   * that n - m doesn't overflow.
   *
   *   long m = -2, n = 10;
   *   long k = m + gmp_urandomm_ui(r, n - m);
   *   inex = mpfr_random_deviate_value(k < 0, k < 0 ? -(k + 1) : k,
   *                                    x, z, r, rnd);
   */
  inex = mpfr_random_deviate_value(0, 0, x, z, r, rnd);
  mpfr_random_deviate_clear(x);
  return inex;
}
#endif
