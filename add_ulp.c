#include "gmp.h"
#include "mpfr.h"

/* sets x to x+sign(x)*ulp(x) */
int mpfr_add_one_ulp(mpfr_ptr x)
{
  int xn, sh; mp_limb_t *xp;

  xn = 1 + (PREC(x)-1)/mp_bits_per_limb;
  sh = xn*mp_bits_per_limb-PREC(x);
  xp = MANT(x);
  if (mpn_add_1(xp, xp, xn, (mp_limb_t)1<<sh)) {
    EXP(x)++;
    mpn_rshift(xp, xp, xn, 1);
    xp[xn-1] += (mp_limb_t)1<<(mp_bits_per_limb-1);
  }
  return 0;
}
