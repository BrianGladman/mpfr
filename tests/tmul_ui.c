#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "time.h"

int
main(int argc, char **argv)
{
  mpfr_t x;
  
  mpfr_init2(x, 53);

  /* checks that result is normalized */
  mpfr_set_d(x, 6.93147180559945286227e-01, GMP_RNDZ);
  mpfr_mul_ui(x, x, 1, GMP_RNDZ);
  if (MANT(x)[PREC(x)/BITS_PER_MP_LIMB] >> (BITS_PER_MP_LIMB-1) == 0) {
    fprintf(stderr, "Error in mpfr_mul_ui: result not normalized\n");
    exit(1);
  }

  mpfr_set_d(x, 1.0/3.0, GMP_RNDZ); 
  mpfr_print_raw(x); putchar('\n'); 
  
  mpfr_mul_ui(x, x, 3, GMP_RNDU); 
  mpfr_print_raw(x); putchar('\n'); 

  printf("%f\n", mpfr_get_d(x)); 

  mpfr_clear(x); 
  return(0);
}
