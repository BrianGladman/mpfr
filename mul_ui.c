#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

void
mpfr_mul_ui(mpfr_ptr y, mpfr_srcptr x, unsigned long u, unsigned char RND_MODE)
{
  mp_limb_t carry, *my, *old_my; unsigned long c; 
  unsigned long xsize, ysize, cnt, dif; 
  TMP_DECL(marker);

  TMP_MARK(marker);
  my = MANT(y); 
  ysize = (PREC(y)-1)/BITS_PER_MP_LIMB + 1;
  xsize = (PREC(x)-1)/BITS_PER_MP_LIMB + 1;

  if (ysize < xsize) {
      old_my = my; 
      my = (mp_ptr) TMP_ALLOC (xsize * BYTES_PER_MP_LIMB);
      dif=0;
    }
  else dif=ysize-xsize;

  carry = mpn_mul_1(my+dif, MANT(x), xsize, u);
  MPN_ZERO(my, dif);

  /* WARNING: count_leading_zeros is undefined for carry=0 */
  if (carry) count_leading_zeros(cnt, carry);
  else cnt=BITS_PER_MP_LIMB;
      
  c = mpfr_round_raw(my, my, RND_MODE, xsize, PREC(y)-BITS_PER_MP_LIMB+cnt);
  
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
  if (ysize < xsize) MPN_COPY(old_my, my, ysize);
  TMP_FREE(marker);
}
