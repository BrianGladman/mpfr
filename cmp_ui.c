/* (c) PolKA project at Inria Lorraine */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

/* returns a positive value if b>i,
           a negative value if b<i,
           zero if b=i 
*/
int mpfr_cmp_ui ( mpfr_srcptr b, unsigned long int i )
{
  int e, k, bn; mp_limb_t c, *bp;

  if (SIGN(b)<0) return -1;
  else if (!NOTZERO(b)) return((i) ? -1 : 0);
  else { /* b>0 */
    e = EXP(b); /* 2^(e-1) <= b < 2^e */
    if (e>mp_bits_per_limb) return 1;

    c = (mp_limb_t) i;
    count_leading_zeros(k, c);
    k = mp_bits_per_limb - k; /* 2^(k-1) <= i < 2^k */
    if (k!=e) return (e-k);

    /* now k=e */
    c <<= (mp_bits_per_limb-k);
    bn = (PREC(b)-1)/mp_bits_per_limb;
    bp = MANT(b) + bn;
    if (*bp>c) return 1;
    else if (*bp<c) return -1;

    /* most significant limbs agree, check remaining limbs from b */
    while (--bn>=0)
      if (*--bp) return 1;
    return 0;
  }
}

/* returns a positive value if b>i,
           a negative value if b<i,
           zero if b=i 
*/
int mpfr_cmp_si ( mpfr_srcptr b, long int i )
{
  int e, k, bn, si; mp_limb_t c, *bp;

  si = (i<0) ? -1 : 1; /* sign of i */
  if (SIGN(b)*i<0) return SIGN(b); /* both signs differ */
  else if (NOTZERO(b)*i==0) { /* one is zero */
    if (i==0) return ((NOTZERO(b)) ? SIGN(b) : 0);
    else return si; /* b is zero */
      
  }
  else { /* b and i are of same sign */
    e = EXP(b); /* 2^(e-1) <= b < 2^e */
    if (e>mp_bits_per_limb) return si;

    c = (mp_limb_t) ((i<0) ? -i : i);
    count_leading_zeros(k, c);
    k = mp_bits_per_limb - k; /* 2^(k-1) <= i < 2^k */
    if (k!=e) return (si*(e-k));

    /* now k=e */
    c <<= (mp_bits_per_limb-k);
    bn = (PREC(b)-1)/mp_bits_per_limb;
    bp = MANT(b) + bn;
    if (*bp>c) return si;
    else if (*bp<c) return -si;

    /* most significant limbs agree, check remaining limbs from b */
    while (--bn>=0)
      if (*--bp) return si;
    return 0;
  }
}

