#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

char *
mpfr_print_rnd_mode(unsigned char rnd_mode)
{
  switch (rnd_mode)
    {
    case GMP_RNDD: return("GMP_RNDD"); 
    case GMP_RNDU: return("GMP_RNDU"); 
    case GMP_RNDN: return("GMP_RNDN"); 
    case GMP_RNDZ: return("GMP_RNDZ"); 
    default: return("unknown rounding mode"); 
    }
}
