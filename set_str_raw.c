#include <stdio.h>
#include <stdlib.h>
#ifdef HAS_STRING_H
#include <string.h>
#else
#include <strings.h>
#endif
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

/* Currently the number should be of the form +/- xxxx.xxxxxxEyy, with 
   decimal exponent. The mantissa of x is supposed to be large enough 
   to hold all the bits of str. */

void
mpfr_set_str_raw(mpfr_ptr x, char *str)
{
  char *str2, *str0, negative = 0; 
  unsigned long j, l, k = 0, xsize, cnt; mp_limb_t *xp; 
  long expn = 0; char *endstr2;

  xp = x -> _mp_d; 
  xsize = 1 + (PREC(x)-1)/BITS_PER_MP_LIMB;
  str0 = str2 = (char *) malloc(strlen(str)*sizeof(char)); 

  if (*str == '-') { negative = 1; str++; }
  else if (*str == '+') str++;

  while (*str == '0') { str++; } 

  while (*str == '0' || *str == '1')
    { *(str2++) = *(str++); k++; }

  if (*str == '.') 
    {
      str++;
      while (*str == '0' || *str == '1')
	{ *(str2++) = *(str++); }

      if (*str == '[') { while (*str != ']') str++; }
    }
    
  if (*str == '[') { while (*str != ']') str++; }

  if (*str == 'e' || *str == 'E') 
    {
      expn = k + atoi(++str); 
      if (expn < k)
	{
	  fprintf(stderr, "Warning : possible overflow in exponent in Str -> mpfr\n"); 
	}
    }
  else expn=k;

  endstr2 = str2;
  l = (strlen(str0) - 1) / BITS_PER_MP_LIMB + 1; str2 = str0; 

  /* str2[0]..endstr2[-1] contains the mantissa */
  for (k = 1; k <= l; k++)
    {
      j = 0; 
      xp[xsize - k] = 0; 
      while (str2<endstr2 && j < BITS_PER_MP_LIMB)
	{
	  xp[xsize - k] = (xp[xsize - k] << 1) + (*str2 - '0'); 
	  str2++; j++; 
	}
      xp[xsize - k] <<= (BITS_PER_MP_LIMB - j); 
    }

  for (; k <= xsize; k++) { xp[xsize - k] = 0; }

  count_leading_zeros(cnt, xp[xsize - 1]); 
  if (cnt) mpn_lshift(xp, xp, xsize, cnt); 

  x -> _mp_exp = expn - cnt; 
  x -> _mp_size = xsize; if (negative) CHANGE_SIGN(x);

  free(str0); 
  
  /* May change to take into account : syntax errors, overflow in exponent, 
     string truncated because of size of x */
}
