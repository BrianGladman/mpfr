#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"

main()
{
   mpfr_t x;

   /* checks that rounds to nearest sets the last
     bit to zero in case of equal distance */
   mpfr_init2(x, 59);
   mpfr_set_str_raw(x, "-0.10010001010111000011110010111010111110000000111101100111111E663"); 
   if (mpfr_can_round(x, 54, GMP_RNDZ, GMP_RNDZ, 53) != 0) {
     fprintf(stderr, "Error in mpfr_can_round\n"); exit(1);
   }
}
