#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void 
#if __STDC__
mpfr_neg(mpfr_t a, mpfr_t b, unsigned char rnd_mode)
#else
mpfr_neg(a, b, rnd_mode)
     mpfr_t a; 
     mpfr_t b; 
     unsigned char rnd_mode; 
#endif
{
  CHANGE_SIGN(b);
  if (a != b) {
    mpfr_set(a, b, rnd_mode);
    CHANGE_SIGN(b);
  }
  return; 
}
