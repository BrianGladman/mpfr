#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

void
mpfr_mul_ui(mpfr_ptr y, mpfr_ptr x, unsigned long u, char RND_MODE)
     /* on suppose SIZ(y)=SIZ(x) */
{
  mp_limb_t carry = 0, *my, *old_my; unsigned long c; 
  unsigned long xsize, ysize, cnt; 

  my = MANT(y); 
  ysize = ABSSIZE(y); xsize = ABSSIZE(x); 

  /*  
      if (ysize < xsize) 
      {
      old_my = my; 
      my = _mp_allocate_func(xsize*BYTES_PER_MP_LIMB);
      _mp_free_func(old_my, ysize*BYTES_PER_MP_LIMB); 
      }
  */
  carry = mpn_mul_1(my, MANT(x), xsize, u); 
  count_leading_zeros(cnt, carry); 
      
  c = mpfr_round_raw(my, my, RND_MODE, ysize, PREC(y)-BITS_PER_MP_LIMB+cnt);
  
  /* If cnt = 1111111111111 and c = 1 we shall get depressed */
  if (c && (carry == (1UL << (BITS_PER_MP_LIMB - cnt)) - 1))
    {
      cnt--; 
      mpn_rshift(my, my, ysize, BITS_PER_MP_LIMB - cnt); 
      my[ysize - 1] |= 1 << (BITS_PER_MP_LIMB - 1);
    }
  else
    {
      mpn_rshift(my, my, ysize, BITS_PER_MP_LIMB - cnt); 
      my[ysize - 1] |= (carry << cnt); 
    }
  EXP(y) = EXP(x) + BITS_PER_MP_LIMB - cnt; 
}
