/* written by Paul Zimmermann, November 1998-January 1999 */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* Remains to do:
- do not use all bits of b and c when PREC(b)>PREC(a) or PREC(c)>PREC(a)
  [current complexity is O(PREC(b)*PREC(c))]
*/

void 
#if __STDC__
mpfr_mul(mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, unsigned char rnd_mode) 
#else
mpfr_mul(a, b, c, rnd_mode) 
     mpfr_ptr a; 
     mpfr_srcptr b; 
     mpfr_srcptr c; 
     unsigned char rnd_mode;
#endif
{
  unsigned int bn, cn, an, k; int cc;
  mp_limb_t *ap=MANT(a), *bp=MANT(b), *cp=MANT(c), *tmp, b1;
  long int sign_product;
  TMP_DECL(marker); 

  /* deal with NaN and zero */
  if (FLAG_NAN(b) || FLAG_NAN(c)) { SET_NAN(a); return; }
  if (!NOTZERO(b) || !NOTZERO(c)) { SET_ZERO(a); return; }

  sign_product = SIGN(b) * SIGN(c);
  bn = (PREC(b)-1)/mp_bits_per_limb+1; /* number of significant limbs of b */
  cn = (PREC(c)-1)/mp_bits_per_limb+1; /* number of significant limbs of c */
  k = bn+cn; /* effective nb of limbs used by b*c */
  TMP_MARK(marker); 
  tmp = (mp_limb_t*) TMP_ALLOC(k*BYTES_PER_MP_LIMB);

  /* multiplies two mantissa in temporary allocated space */
  b1 = (bn>=cn) ? mpn_mul(tmp, bp, bn, cp, cn) : mpn_mul(tmp, cp, cn, bp, bn);

  /* now tmp[0]..tmp[k-1] contains the product of both mantissa,
     with tmp[k-1]>=2^(mp_bits_per_limb-2) */
  an = (PREC(a)-1)/mp_bits_per_limb+1; /* number of significant limbs of a */
  b1 >>= mp_bits_per_limb-1; /* msb from the product */

  if (b1==0) mpn_lshift(tmp, tmp, k, 1);
  cc = mpfr_round_raw(ap, tmp, PREC(b)+PREC(c), (sign_product<0), PREC(a), rnd_mode);
  if (cc) { /* cc = 1 ==> result is a power of two */
    ap[an-1] = (mp_limb_t) 1 << (BITS_PER_MP_LIMB-1);
  }
  EXP(a) = EXP(b) + EXP(c) + b1 - 1 + cc;
  if (sign_product * SIGN(a)<0) CHANGE_SIGN(a);
  TMP_FREE(marker); 
  return;
}

     
