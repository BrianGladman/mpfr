#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

/* Computes a random mpfr in [0, 1[ with precision PREC */

extern long random _PROTO((void)); 
extern int srandom _PROTO((unsigned int)); 

/* extracted from GNU mpf */
#if defined (__hpux) || defined (__alpha)
/* HPUX lacks random().  DEC OSF/1 1.2 random() returns a double.  */
#define random mrand48
#define srandom srand48
#endif

void
#if __STDC__
mpfr_random(mpfr_ptr x)
#else
mpfr_random(x)
     mpfr_ptr x; 
#endif    
{
  mp_limb_t *xp; unsigned long xs, i, cnt; 

  xp = MANT(x); 
  xs = ABSSIZE(x); 

  for (i = 0; i < xs; i++)
    {
      /* random() c/sh/ould be replaced by a homemade random number generator.
	 Indeed, if on Linux random is a good RNG, this is definitely not
	 the case in most Un*xes. */
      xp[i] = random();
    }
  
  count_leading_zeros(cnt, xp[xs - 1]); 
  if (cnt) mpn_lshift(xp, xp, xs, cnt); 
  EXP(x) = -cnt; 
  mpfr_round_raw(xp, xp, PREC(x), (SIGN(x)<0), PREC(x), random()&3);
  /* ignore any possible carry (this is sheer laziness). */
}

void
mpfr_srandom(unsigned long seed)
{
  srandom(seed); 
}
