#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void
#if __STDC__
mpfr_init2 (mpfr_t x, unsigned long int p)
#else
mpfr_init2 (x, p)
     mpfr_t x;
     unsigned long int p;
#endif
{
  unsigned long xsize; 

  if (p==0) {
    printf("*** cannot initialize mpfr with precision 0\n"); exit(1);
  }

  xsize = (p - 1)/BITS_PER_MP_LIMB + 1; 

  x -> _mp_prec = p;
  x -> _mp_d = (mp_ptr) (*_mp_allocate_func) 
    (xsize * BYTES_PER_MP_LIMB);
  x -> _mp_size = xsize;
  x -> _mp_exp = 0; /* avoids uninitialized memory reads for zero */
}
