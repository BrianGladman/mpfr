#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "longlong.h"

/* #define DEBUG */

void
mpfr_sqrt3 (mpfr_ptr r, mpfr_srcptr u, unsigned char rnd_mode)
{
  mp_ptr up, rp, tmp;
  mp_size_t usize, rrsize;
  mp_size_t rsize;
  mp_size_t prec, err;
  mp_limb_t q_limb;
  long rw, nw; 
  unsigned long cc = 0; 
  char can_round = 0; 
  TMP_DECL (marker); TMP_DECL(marker0); 

  if (FLAG_NAN(u) || SIGN(u) == -1) { SET_NAN(r); return; }
  
  usize = ABSSIZE(u); 
  prec = PREC(r);

  if (!NOTZERO(u))
    {
      EXP(r) = 0; 
      MPN_ZERO(MANT(r), SIZE(r)); 
      return; 
    }

  up = MANT(u);

#ifdef DEBUG
      printf("Entering square root : "); 
      for(k = usize - 1; k >= 0; k--) { printf("%lu ", up[k]); }
      printf(".\n"); 
#endif

  /* Compare the mantissas */
  
  EXP(r) = ((EXP(u) + (EXP(u) & 1)) / 2) ;  
  
  rsize = ((PREC(r) + 2 + (EXP(u) & 1))/BITS_PER_MP_LIMB + 1) << 1; 
  rrsize = (PREC(r) + 2 + (EXP(u) & 1))/BITS_PER_MP_LIMB + 1;
  /* One extra bit is needed in order to get the square root with enough
     precision ; take one extra bit for rrsize in order to solve more 
     easily the problem of rounding to nearest.
     Need to have 2*rrsize = rsize...
     Take one extra bit if the exponent of u is odd since we shall have
     to shift then.
  */

  TMP_MARK(marker0); 
  if (EXP(u) & 1) /* Shift u one bit to the right */
    {
      up = TMP_ALLOC((SIZE(u) + (SIZE(u)*BITS_PER_MP_LIMB == PREC(u)))
		     *BYTES_PER_MP_LIMB);   
  
      /* NE MARCHE QUE SI LA PARTIE NON SIGNIFICATIVE DE u EST A ZERO */
      /* LE CONFIRMER ENCORE UNE FOIS ... FOOLPROOF ? */

      if (mpn_rshift(up, u->_mp_d, SIZE(u), 1))
	up [0] = ((mp_limb_t) 1) << (BITS_PER_MP_LIMB - 1); 
    }

  do
    {
      TMP_MARK (marker);

      err = rsize*BITS_PER_MP_LIMB; 
      if (rsize < usize) { err--; }
      if (err > rrsize * BITS_PER_MP_LIMB) 
	{ err = rrsize * BITS_PER_MP_LIMB; }
      
      tmp = (mp_ptr) TMP_ALLOC (rsize * BYTES_PER_MP_LIMB);  
      rp = (mp_ptr) TMP_ALLOC (rrsize * BYTES_PER_MP_LIMB); 

      if (usize >= rsize) { 
	MPN_COPY (tmp, up + usize - rsize, rsize);
      }
      else { 
	MPN_COPY (tmp + rsize - usize, up, usize);
	MPN_ZERO (tmp, rsize - usize); 
      }

      /* Do the real job */
 
#ifdef DEBUG
      printf("Taking the sqrt of : "); 
      for(k = rsize - 1; k >= 0; k--) { printf("%lu ", tmp[k]); }
      printf(".\n"); 
#endif

      q_limb = mpn_sqrtrem (rp, NULL, tmp, rsize);

#ifdef DEBUG
      printf("The result is : \n"); 
      printf("sqrt : "); 
      for(k = rrsize - 1; k >= 0; k--) { printf("%lu ", rp[k]); }
      printf("(q_limb = %lu)\n", q_limb); 
#endif
      
      can_round = (mpfr_can_round_raw(rp, rrsize, 1, err, 
				      GMP_RNDZ, rnd_mode, PREC(r))); 

      /* If we used all the limbs of both the dividend and the divisor, 
	 then we have the correct RNDZ rounding */

      if (!can_round && (rsize < 2*usize)) 
	{ 
#ifdef DEBUG
	  printf("Increasing the precision.\n"); 
#endif
	  printf("#"); 
	  TMP_FREE(marker); 
	}
    }
  while (!can_round && (rsize < 2*usize) 
	 && (rsize += 2) && (rrsize ++)); 

  if (can_round) 
    {
      cc = mpfr_round_raw(rp, rp, err, 0, PREC(r), rnd_mode);  
      rrsize = (PREC(r) - 1)/BITS_PER_MP_LIMB + 1; 
    }
  else
    /* Use the return value of sqrtrem to decide of the rounding         */
    /* Note that at this point the sqrt has been computed                */
    /* EXACTLY. If rounding = GMP_RNDZ, do nothing [comes from           */
    /* the fact that the exact square root can end with a bunch of ones, */
    /* and in that case we indeed cannot round if we do not know that    */
    /* the computation was exact.                                        */
    switch (rnd_mode)
      {
      case GMP_RNDZ : 
      case GMP_RNDD : break; 

      case GMP_RNDN : 
	/* Not in the situation ...0 111111 */
	rw = (PREC(r) + 1) & (BITS_PER_MP_LIMB - 1);
	if (rw) { rw = BITS_PER_MP_LIMB - rw; nw = 0; } else nw = 1; 
	if ((rp[nw] >> rw) & 1 &&                     /* Not 0111111111 */
	    (q_limb ||                                /* Nonzero remainder */
	    (rw ? (rp[nw] >> (rw - 1)) & 1 : 
	     (rp[nw] >> (BITS_PER_MP_LIMB - 1)) & 1))) /* or even rounding */ 
	  cc = mpn_add_1(rp + nw, rp + nw, rrsize, ((mp_limb_t)1) << rw); 
	break;
 
      case GMP_RNDU : 
	if (q_limb)
	  cc = mpn_add_1(rp, rp, rrsize, 1 << (BITS_PER_MP_LIMB - 
					       (PREC(r) & 
						(BITS_PER_MP_LIMB - 1))));
      }

  if (cc) {
    mpn_rshift(rp, rp, rrsize, 1);
    rp[rrsize-1] |= (mp_limb_t) 1 << (BITS_PER_MP_LIMB-1);
    r->_mp_exp++; 
  }
    
  rp [0] &= ~(((mp_limb_t)1 << (BITS_PER_MP_LIMB - 
		     (PREC(r) & (BITS_PER_MP_LIMB - 1)))) - 1) ; 
  
  rsize = rrsize; 
  rrsize = (PREC(r) - 1)/BITS_PER_MP_LIMB + 1;  
  MPN_COPY(r->_mp_d + SIZE(r) - rrsize, rp + rsize - rrsize, rrsize); 
  TMP_FREE(marker0); TMP_FREE (marker);
}
