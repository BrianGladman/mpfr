#include "gmp.h"
#include "gmp-impl.h"

char __gmp_default_rounding_mode = 0;

void
#if __STDC__
mpfr_set_default_rounding_mode (char rnd_mode)
#else
mpfr_set_default_rounding_mode (rnd_mode)
     char rnd_mode;
#endif
{
  __gmp_default_rounding_mode = rnd_mode;
}

