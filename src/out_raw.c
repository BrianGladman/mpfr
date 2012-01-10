This is work in progress. Before continuing, please read
http://websympa.loria.fr/wwsympa/arc/mpfr/2011-01/msg00039.html.

/* mpfr_out_raw -- output a floating-point number to binary portable format

Copyright 2011, 2012 Free Software Foundation, Inc.
Contributed by the Arenaire and Caramel projects, INRIA.

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

/* format for out_raw:
   A mpfr_t is represented by up to 3 fields, each one is represented by a
   sequence of 32-bit words. 32-bit words are stored as 4 bytes in little
   endian format:
   (a) a field for sign, precision and other bit-fields:
       - the sign (1 bit)
       - 2 bits for NaN, Inf, 0, normal numbers:
         00: 0
         01: normal
         10: Inf
         11: NaN
       - 1 bit for the exponent encoding (exp_enc) for normal numbers,
         otherwise 0 is stored here (reserved for future extensions)
       - 1 bit for the precision encoding (prec_enc)
       If prec_enc=0, the remaining 27 bits encode the precision (< 2^27)
       If prec_enc=1, the precision is stored in the following 27 bits
       (high part) and then 32 bits (low part). Thus the maximal precision
       is 59 bits.
   (b) (optional) a field for the exponent:
       - if the number is NaN, Inf, 0, this field is empty
       - if exp_enc=0, this field contains one 32-bit (signed) word encoding
         the exponent
       - if exp_enc=1, a first 32-bit word encodes a positive integer m,
         and the following m 32-bit words encode the exponent (in 2-complement
         representation, with least significant words first)
   (c) (optional) a field for the significand:
       - if the number is NaN, Inf, 0, this field is empty
       - otherwise, let p = ceil(prec/32), the significand is represented
         by p consecutive 32-bit words (least significant words first).
         Thus on a little-endian machine the significand can be directly
         copied using memcopy.
   Examples:
   - a normal binary32 IEEE-754 number uses 96 bits: 32 for (a), 32 for (b),
     and 32 for (c);
   - a normal binary64 IEEE-754 number uses 128 bits: 32 for (a), 32 for (b),
     and 64 for (c) (idem for a significand of 64 bits, as in Intel x86
     double-extended format);
   - a normal binary128 IEEE-754 number uses 192 bits: 32 for (a), 32 for (b),
     and 128 for (c).
 */

size_t
mpfr_out_raw (FILE *stream, mpfr_srcptr x)
{
  size_t n; /* number of bytes of the output */
  mpfr_prec_t prec = MPFR_PREC(x);
  int prec_enc, exp_enc = 0;
  unsigned char *s, *t;

  MPFR_ASSERTN (CHAR_BIT == 8);

  prec_enc = (prec >= 134217728);
  n = 4 + 4 * prec_enc;
  if (MPFR_LIKELY (! MPFR_IS_SINGULAR(x)))
    {
      mpfr_exp_t e = MPFR_EXP(x);
      mpfr_prec_t p = (prec - 1) / 32 + 1; /* ceil(prec/32) */

      exp_enc = (e < -2147483648 || 2147483647 < e);
      n += 4 + 4 * exp_enc + 4 * p;
    }
  t = s = (unsigned char *) malloc (n * sizeof (char));
  t[0] = (mpfr_signbit (x) != 0) << 7;
  if (MPFR_IS_NAN(x))
    t[0] += 3 << 5;
  else if (MPFR_IS_INF(x))
    t[0] += 2 << 5;
  else if (MPFR_IS_ZERO(x))
    t[0] += 1 << 5;
  t[0] += exp_enc << 4;
  t[0] += prec_enc << 3;
  if (prec_enc)
    {
      /* FIXME: the shift count may be too large for the type size. */
      /* FIXME: check that the precision is not too large. */
      t[0] += prec >> 56; /* reduction mod 8 is implicit */
      t[1] = prec >> 48;
      t[2] = prec >> 40;
      t[3] = prec >> 32;
      t += 4;
      t[0] = 0;
    }
  t[0] += prec >> 24;  /* reduction mod 256 is implicit */
  t[1] = prec >> 16;
  t[2] = prec >> 8;
  t[3] = prec;
  free (s);
  return n;
}
