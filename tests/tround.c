#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"

main()
{
   mpfr_t x;

   /* checks that rounds to nearest sets the last
     bit to zero in case of equal distance */
   mpfr_init2(x, 2);
   mpfr_set_d(x, 5.0, 0);
   if (mpfr_get_d(x) != 4.0) { printf("Error in tround: got %1.1f instead of 4.0\n",mpfr_get_d(x)); }

   mpfr_set_d(x, 0.00098539467465030839, 0);

   mpfr_set_d(x, 9.84891017624509146344e-01, GMP_RNDU); 
   if (mpfr_get_d(x) != 1.0) { printf("Error in tround: got %f instead of 1.0\n",mpfr_get_d(x)); exit(1); }

   mpfr_clear(x); 
}
