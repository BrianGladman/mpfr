#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

void
mpfr_set_si(mpfr_ptr x, long i, unsigned char rnd_mode)
{
  unsigned long xn, ai, cnt; 

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
mpfr_set_ui(mpfr_ptr x, unsigned long i, unsigned char rnd_mode)
{
  unsigned int xn, cnt; 

  if (i==0) { SET_ZERO(x); return; }
  xn = (PREC(x)-1)/BITS_PER_MP_LIMB;
  count_leading_zeros(cnt, i); 

  x -> _mp_d[xn] = i << cnt; 
  /* don't forget to put zero in lower limbs */
  MPN_ZERO(MANT(x), xn);
  x -> _mp_exp = BITS_PER_MP_LIMB - cnt;
  /* warning: don't change the precision of x! */
  if (i*SIGN(x) < 0) CHANGE_SIGN(x);

  return; 
}

