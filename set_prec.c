#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void
#if __STDC__
mpfr_set_prec (mpfr_t x, unsigned long int p, unsigned char rnd_mode)
#else
mpfr_set_prec (x, p, rnd_mode)
     mpfr_t x;
     unsigned long int p;
     unsigned char rnd_mode;
#endif
{
  unsigned long xsize,oldp,oldsize; mp_limb_t *old;

  if (p==0) {
    printf("*** cannot set precision to 0 bits\n"); exit(1);
  }

  xsize = (p - 1)/BITS_PER_MP_LIMB + 1; /* new limb size */

  x -> _mp_d = (mp_ptr) (*_mp_reallocate_func) 
    (x -> _mp_d, ABSSIZE(x)*BYTES_PER_MP_LIMB, xsize * BYTES_PER_MP_LIMB);
  x -> _mp_prec = p;
  SIZE(x) = (SIGN(x)>0) ? xsize : (xsize ^ (1<<31));
}

unsigned long int
#if __STDC__
mpfr_get_prec (mpfr_t x)
#else
mpfr_set_prec (x)
     mpfr_t x;
#endif
{
  return x -> _mp_prec;
}
