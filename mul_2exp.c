#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void
mpfr_mul_2exp(mpfr_ptr y, mpfr_srcptr x, unsigned long int n, unsigned char rnd_mode)
{
  /* Important particular case */ 
  if (y != x) mpfr_set(y, x, rnd_mode);
  EXP(y) += n;
  return;
}

