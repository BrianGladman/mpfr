#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include "time.h"

int
main(int argc, char **argv)
{
  mpfr_t x;
  
  mpfr_init2(x, 53);
  mpfr_set_d(x, 1.0/3.0, GMP_RNDZ); 
  mpfr_print_raw(x); putchar('\n'); 
  
  mpfr_mul_ui(x, x, 3, GMP_RNDU); 
  mpfr_print_raw(x); putchar('\n'); 

  printf("%f\n", mpfr_get_d(x)); 

  mpfr_clear(x); 
  return(0);
}
