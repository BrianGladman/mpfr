#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"

main()
{
   mpfr_t x;

   mpfr_init(x);

   mpfr_set_d(x, 1.0, GMP_RNDN);
   mpfr_abs(x, x, GMP_RNDN);
   if (mpfr_get_d(x) != 1.0) {
     fprintf(stderr, "Error in mpfr_abs(1.0)\n"); exit(1);
   }

   mpfr_set_d(x, -1.0, GMP_RNDN);
   mpfr_abs(x, x, GMP_RNDN);
   if (mpfr_get_d(x) != 1.0) {
     fprintf(stderr, "Error in mpfr_abs(-1.0)\n"); exit(1);
   }

   mpfr_clear(x); 
}
