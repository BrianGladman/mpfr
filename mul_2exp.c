#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void
#if __STDC__
mpfr_mul_2exp(mpfr_ptr y, mpfr_srcptr x, unsigned long int n, unsigned char rnd_mode)
#else
mpfr_mul_2exp(y, x, n, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     unsigned long int n; 
     unsigned char rnd_mode; 
#endif
{
  /* Important particular case */ 
  if (y != x) mpfr_set(y, x, rnd_mode);
  EXP(y) += n;
  return;
}

