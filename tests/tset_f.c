#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include "time.h"

int
main()
{
  mpfr_t x; mpf_t y; mpf_t z; unsigned long k, pr; 
  
  mpfr_init2(x, 100);
  mpf_init(y); 

  srandom(time(NULL)); 
  mpf_random2(y, 10, 0); 
  mpfr_set_f(x, y, 53, rand() & 3); 

  mpf_clear(y); mpfr_clear(x); 

  for (k = 1; k <= 100000; k++)
    {
      pr = 1 + (rand()&255); 
      mpf_init2(z, pr);
      mpf_random2(z, z->_mp_prec, 0);
      mpfr_init2(x, pr);
      mpfr_set_f(x, z, pr, 0);
    }
  return(0);
}
