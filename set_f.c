#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

void mpfr_set_f(mpfr_ptr y, mpf_srcptr x, unsigned long prec, char rnd_mode)
{
  mp_limb_t *my, *mx; unsigned long cnt, sx, sy; 

  sx = SIZ (x); sy = ABSSIZE(y);
  my = MANT(y); mx = MANT(x); 

  count_leading_zeros(cnt, mx[sx - 1]);  

  if (sy < sx)
    {
      if (cnt) 
	{ 
	  mpn_lshift(my, mx + 1, sy, cnt); 
	  my [0] |= mx[0] >> (BITS_PER_MP_LIMB - cnt); 
	  EXP(y) = EXP(x)* BITS_PER_MP_LIMB - cnt; 
	}
      else { MPN_COPY(my, mx, sy); EXP(y) = EXP(x) * BITS_PER_MP_LIMB; }
    }
  else
    {
      if (cnt) 
	{ 
	  mpn_lshift(my, mx, sx, cnt); 
	  EXP(y) = EXP(x)* BITS_PER_MP_LIMB - cnt; 
	}
      else { MPN_COPY(my, mx, sy); EXP(y) = EXP(x) * BITS_PER_MP_LIMB; }
    }
  
  mpfr_round(y, rnd_mode, prec); 
}
