#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void
#if __STDC__
mpfr_get_str_raw(char *digit_ptr, mpfr_srcptr x)
#else
mpfr_get_str_raw(digit_ptr, x)
     char *digit_ptr; 
     mpfr_srcptr x; 
#endif
{
  mp_limb_t *mx, wd, t; long ex, sx, k, l, p;

  mx = MANT(x); 
  ex = EXP(x); 
  p = PREC(x); 

  if (SIGN(x) < 0) { *digit_ptr = '-'; digit_ptr++; }
  sprintf(digit_ptr, "0."); digit_ptr += 2; 

  sx = 1+(p-1)/mp_bits_per_limb; /* number of significant limbs */
  for (k = sx - 1; k >= 0 ; k--)
    { 
      wd = mx[k]; 
      t = 1UL << (BITS_PER_MP_LIMB - 1); 
      for (l = BITS_PER_MP_LIMB - 1; l>=0; l--)
	{
	  if (wd & t) 
	    { *digit_ptr = '1'; digit_ptr++; } 
	  else 
	    { *digit_ptr = '0'; digit_ptr++; }
	  t >>= 1; 
	  if (--p==0) { *digit_ptr = '['; digit_ptr++; }
	}
    }
  sprintf(digit_ptr, "]E%ld", ex); 
}
 
void
#if __STDC__
mpfr_print_raw(mpfr_srcptr x)
#else
mpfr_print_raw(x)
     mpfr_srcptr x; 
#endif
{
  char *str; 

  if (FLAG_NAN(x)) printf("NaN");
  else if (!NOTZERO(x)) printf("0");
  else {
     /* 3 char for sign + 0 + binary point
	+ ABSSIZE(x) * BITS_PER_MP_LIMB for mantissa
	+ 2 for brackets in mantissa
	+ 1 for 'E'
	+ 11 for exponent (including sign)
	= 17 + ABSSIZE(x) * BITS_PER_MP_LIMB
      */
     str = (char *) malloc((17 + ABSSIZE(x) * BITS_PER_MP_LIMB)*sizeof(char));
     mpfr_get_str_raw(str, x);

     printf("%s", str); 
     free(str); 
  }
}

