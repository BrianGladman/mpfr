#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void mpfr_set(a, b, rnd_mode) 
mpfr_ptr a; mpfr_srcptr b; unsigned char rnd_mode;
{
  mpfr_round_raw(MANT(a), b->_mp_d, rnd_mode, b->_mp_size, PREC(a));
  EXP(a) = EXP(b);
  if (SIGN(a) != SIGN(b)) CHANGE_SIGN(a);
}
