#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void
#if __STDC__
mpfr_set_prec (mpfr_t x, unsigned long int p)
#else
mpfr_set_prec (x, p)
     mpfr_t x;
     unsigned long int p;
#endif
{
  unsigned long xsize;

  if (p==0) {
    printf("*** cannot set precision to 0 bits\n"); exit(1);
  }

  xsize = (p - 1)/BITS_PER_MP_LIMB + 1; /* new limb size */

  if (xsize > ABSSIZE(x)) {
    x -> _mp_d = (mp_ptr) (*_mp_reallocate_func) 
      (x -> _mp_d, ABSSIZE(x)*BYTES_PER_MP_LIMB, xsize * BYTES_PER_MP_LIMB);
    SIZE(x) = xsize; /* new number of allocated limbs */
  }

  x -> _mp_prec = p;
}

unsigned long int
#if __STDC__
mpfr_get_prec (mpfr_t x)
#else
mpfr_get_prec (x)
     mpfr_t x;
#endif
{
  return x -> _mp_prec;
}
