#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

void
mpfr_set_si(mpfr_ptr x, long i, unsigned char rnd_mode)
{
  unsigned long xsize, ai, cnt; 

  if (i==0) { SET_ZERO(x); return; }
  xsize = ABSSIZE(x);
  ai = ABS(i); 

  count_leading_zeros(cnt, ai); 

  x -> _mp_d[xsize - 1] = ai << cnt;
  /* don't forget to put zero in lower limbs */
  MPN_ZERO(MANT(x), xsize-1);
  x -> _mp_exp = BITS_PER_MP_LIMB - cnt;
  /* warning: don't change the precision of x! */
  x -> _mp_size = xsize; if (i < 0) CHANGE_SIGN(x);

  return; 
}

void
mpfr_set_ui(mpfr_ptr x, unsigned long i, unsigned char rnd_mode)
{
  unsigned int xsize, cnt; 

  if (i==0) { SET_ZERO(x); return; }
  xsize = ABSSIZE(x);
  count_leading_zeros(cnt, i); 

  x -> _mp_d[xsize - 1] = i << cnt; 
  /* don't forget to put zero in lower limbs */
  MPN_ZERO(MANT(x), xsize-1);
  x -> _mp_exp = BITS_PER_MP_LIMB - cnt;
  /* warning: don't change the precision of x! */
  x -> _mp_size = xsize; 

  return; 
}

