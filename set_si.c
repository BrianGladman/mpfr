#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

void
#if __STDC__
mpfr_set_si(mpfr_ptr x, long i, unsigned char rnd_mode)
#else
mpfr_set_si(x, i, rnd_mode)
     mpfr_ptr x;
     long i;
     unsigned char rnd_mode;
#endif
{
  unsigned long xn, cnt; mp_limb_t ai;

  if (i==0) { SET_ZERO(x); return; }
  xn = (PREC(x)-1)/BITS_PER_MP_LIMB;
  ai = ABS(i); 

  count_leading_zeros(cnt, ai); 

  x -> _mp_d[xn] = ai << cnt;
  /* don't forget to put zero in lower limbs */
  MPN_ZERO(MANT(x), xn);
  x -> _mp_exp = BITS_PER_MP_LIMB - cnt;
  /* warning: don't change the precision of x! */
  if (i*SIGN(x) < 0) CHANGE_SIGN(x);

  return; 
}

void
#if __STDC__
mpfr_set_ui(mpfr_ptr x, unsigned long i, unsigned char rnd_mode)
#else
mpfr_set_ui(x, i, rnd_mode)
     mpfr_ptr x;
     long i;
     unsigned char rnd_mode;
#endif  
{
  unsigned int xn, cnt;

  if (i==0) { SET_ZERO(x); return; }
  xn = (PREC(x)-1)/BITS_PER_MP_LIMB;
  count_leading_zeros(cnt, (mp_limb_t) i); 

  x -> _mp_d[xn] = ((mp_limb_t) i) << cnt; 
  /* don't forget to put zero in lower limbs */
  MPN_ZERO(MANT(x), xn);
  x -> _mp_exp = BITS_PER_MP_LIMB - cnt;
  /* warning: don't change the precision of x! */
  if (i*SIGN(x) < 0) CHANGE_SIGN(x);

  return; 
}

