#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void 
#if __STDC__
mpfr_neg(mpfr_ptr a, mpfr_srcptr b, unsigned char rnd_mode)
#else
mpfr_neg(a, b, rnd_mode)
     mpfr_ptr a; 
     mpfr_srcptr b; 
     unsigned char rnd_mode; 
#endif
{
  if (a != b) mpfr_set4(a, b, rnd_mode, -SIGN(b));
  else CHANGE_SIGN(a);
  return; 
}
