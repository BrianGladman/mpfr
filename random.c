#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* Computes a random mpfr in [0, 1[ with precision PREC */

void
mpfr_random(mpfr_ptr x)
{
  mp_limb_t *xp; unsigned long xs, i; 

  EXP(x) = 0; 
  xp = MANT(x); 
  xs = ABSSIZE(x); 

  for (i = 0; i < xs; i++)
    {
      /* random() could be replaced by a homemade random number generator.
	 Indeed, if on Linux random is a good RNG, this is definitely not
	 the case in most Un*xes. */
      xp[i] = random();
    }
  
  mpfr_round_raw(xp, xp, random()&3, SIZE(x), PREC(x)); 
  /* Since the value 1 is forbidden, ignore any possible carry. */
}

void
mpfr_srandom(unsigned long seed)
{
  srandom(seed); 
}
