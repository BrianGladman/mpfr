#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gmp.h"
#include "mpfr.h"

/* checks that x*y gives the same results in double
   and with mpfr with 53 bits of precision */

int
main(argc,argv) int argc; char *argv[];
{
  double x, z; mpfr_t w;

  mpfr_init2(w, 53); 

  srand48(time(NULL)); 
  x = drand48(); 
  mpfr_set_d(w, x, 0);
  mpfr_mul_2exp(w, w, 10, GMP_RNDZ); 
  if (x != (z = mpfr_get_d(w)/1024))
    {
      fprintf(stderr, "%f != %f\n", x, z); 
      return (-1); 
    };
  return (0); 
}

