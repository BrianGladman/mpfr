#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void mpfr_neg(mpfr_t a, mpfr_t b, unsigned char rnd_mode)
{
  CHANGE_SIGN(b);
  if (a != b) {
    mpfr_set(a, b, rnd_mode);
    CHANGE_SIGN(b);
  }
  return; 
}
