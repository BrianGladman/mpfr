
#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void 
#if __STDC__
mpfr_set(mpfr_ptr a, mpfr_srcptr b, unsigned char rnd_mode)
#else
mpfr_set(a, b, rnd_mode) 
     mpfr_ptr a; 
     mpfr_srcptr b; 
     unsigned char rnd_mode;
#endif
{
  int carry, an, preca = PREC(a), sh; mp_limb_t *ap = MANT(a);

  carry = mpfr_round_raw(ap, MANT(b), PREC(b), (SIGN(b)<0), preca, rnd_mode);
  EXP(a) = EXP(b);
  if (carry) {
    an = (preca-1)/BITS_PER_MP_LIMB + 1;
    sh = an * BITS_PER_MP_LIMB - preca;
    if ((*ap >> sh) & 1) {
      fprintf(stderr, "unable to round in mpfr_set\n"); exit(1);
    }
    mpn_rshift(ap, ap, an, 1);
    ap[an-1] |= (mp_limb_t) 1 << (BITS_PER_MP_LIMB-1);
    EXP(a)++;
  }
  if (SIGN(a) != SIGN(b)) CHANGE_SIGN(a);
}
