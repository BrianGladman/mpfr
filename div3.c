#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "longlong.h"

/* #define DEBUG */

void
mpfr_div3 (mpfr_ptr r, mpfr_srcptr u, mpfr_srcptr v, unsigned char rnd_mode)
{
  mp_srcptr up, vp;
  mp_ptr rp, tp, tmp;
  mp_size_t usize, vsize, rrsize;
  mp_size_t rsize;
  mp_size_t sign_quotient;
  mp_size_t prec, err;
  mp_limb_t q_limb;
  mp_exp_t rexp;
  long k, mult, vn; 
  unsigned long cc = 0; 
  char can_round = 0; 
  TMP_DECL (marker);

  if (FLAG_NAN(u) || FLAG_NAN(v)) { SET_NAN(r); return; }
  
  usize = ABSSIZE(u); 
  vsize = ABSSIZE(v); 
  sign_quotient = (SIGN(u) == SIGN(v) ? 1 : -1); 
  prec = PREC(r);

  if (!NOTZERO(v))
    vsize = 1 / v->_mp_d[vsize - 1];    /* Gestion des infinis ? */
  
  if (MANT(u)[(PREC(u)-1)/BITS_PER_MP_LIMB] == 0)
    {
      r->_mp_exp = 0;
      MPN_ZERO(r->_mp_d, r->_mp_size); 
      return;
    }

  up = u->_mp_d;
  vp = v->_mp_d;

#ifdef DEBUG
      printf("Entering division : "); 
      for(k = usize - 1; k >= 0; k--) { printf("%lu ", up[k]); }
      printf(" by "); 
      for(k = vsize - 1; k >= 0; k--) { printf("%lu ", vp[k]); }
      printf(".\n"); 
#endif

  /* Compare the mantissas */
  mult = mpn_cmp(up, vp, (usize > vsize ? vsize : usize)); 
  if (mult == 0 && vsize > usize)
    {
      vn = vsize - usize; 
      while (vn >= 0) if (vp[vn--]) { mult = 1; break; }
      /* On peut diagnostiquer ici pour pas cher le cas u = v */
    }
  else { mult = (mult < 0 ? 1 : 0); }

  rsize = (prec + 3)/BITS_PER_MP_LIMB + 1; 
  rrsize = prec/BITS_PER_MP_LIMB + 1;
  /* Three extra bits are needed in order to get the quotient with enough
     precision ; take one extra bit for rrsize in order to solve more 
     easily the problem of rounding to nearest. */

  /* ATTENTION, USIZE DOIT RESTER > A VSIZE !!!!!!!! */

  do
    {
      TMP_MARK (marker);

      rexp = u->_mp_exp - v->_mp_exp;
      
      err = rsize*BITS_PER_MP_LIMB; 
      if (rsize < vsize) { err-=2; } 
      if (rsize < usize) { err--; }
      if (err > rrsize * BITS_PER_MP_LIMB) 
	{ err = rrsize * BITS_PER_MP_LIMB; }
      
      tp = (mp_ptr) TMP_ALLOC (rsize * BYTES_PER_MP_LIMB);  
      tmp = (mp_ptr) TMP_ALLOC (rsize * BYTES_PER_MP_LIMB);
      rp = (mp_ptr) TMP_ALLOC (rrsize * BYTES_PER_MP_LIMB); 

      if (vsize >= rsize) { 
	MPN_COPY (tmp, vp + vsize - rsize, rsize);
      }
      else { 
	MPN_COPY (tmp + rsize - vsize, vp, vsize);
	MPN_ZERO (tmp, rsize - vsize); 
      }

      if (usize >= rsize) { 
	MPN_COPY (tp, up + usize - rsize, rsize);
      }
      else {
	MPN_COPY (tp + rsize - usize, up, usize); 
	MPN_ZERO (tp, rsize - usize); 
      }

      /* Do the real job */
 
#ifdef DEBUG
      printf("Dividing : "); 
      for(k = rsize - 1; k >= 0; k--) { printf("%lu ", tp[k]); }
      printf(" by "); 
      for(k = rsize - 1; k >= 0; k--) { printf("%lu ", tmp[k]); }
      printf(".\n"); 
#endif

      q_limb = mpn_divrem (rp, rrsize, tp, rsize, tmp, rsize);

#ifdef DEBUG
      printf("The result is : \n"); 
      printf("Quotient : "); 
      for(k = rrsize - 1; k >= 0; k--) { printf("%lu ", rp[k]); }
      printf("Remainder : "); 
      for(k = rsize - 1; k >= 0; k--) { printf("%lu ", tp[k]); }
      printf("(q_limb = %lu)\n", q_limb); 
#endif
      
      /* msb-normalize the result */
      
      if (q_limb)
	{
	  count_leading_zeros(k, q_limb); 
	  mpn_rshift(rp, rp, rrsize, BITS_PER_MP_LIMB - k); 
	  rp[rrsize - 1] |= (q_limb << k);  
	  rexp += BITS_PER_MP_LIMB - k; 
	}
      else
	{ 
	  count_leading_zeros(k, rp[rrsize - 1]); 
	  if (k) { mpn_lshift(rp, rp, rrsize, k); }
	  rexp -= k; 
	}
      
      can_round = (mpfr_can_round_raw(rp, rrsize, sign_quotient, err, 
				     GMP_RNDN, rnd_mode, prec)
	|| (usize == rsize && vsize == rsize && 
	    mpfr_can_round_raw(rp, rrsize, sign_quotient, err, 
			       GMP_RNDZ, rnd_mode, prec))); 

      /* If we used all the limbs of both the dividend and the divisor, 
	 then we have the correct RNDZ rounding */

      if (!can_round && (rsize < usize || rsize < vsize)) 
	{ 
#ifdef DEBUG
	  printf("Increasing the precision.\n"); 
#endif
	  printf("#"); 
	  TMP_FREE(marker); 
	}
    }
  while (!can_round && (rsize < usize || rsize < vsize) 
	 && (rsize++) && (rrsize++)); 

  /* ON PEUT PROBABLEMENT SE DEBROUILLER DES QUE rsize >= vsize */
  /* MAIS IL FAUT AJOUTER LE BOUT QUI MANQUE DE usize A rsize */
    
  if (can_round) 
    {
      cc = mpfr_round_raw(rp, rp, err, (sign_quotient == -1 ? 1 : 0),
			  prec, rnd_mode);  
      rrsize = (prec - 1)/BITS_PER_MP_LIMB + 1; 
    }
  else
    /* Use the remainder to find out the correct rounding */
    /* Note that at this point the division has been done */
    /* EXACTLY. */
    if ((rnd_mode == GMP_RNDD && sign_quotient == -1) 
	|| (rnd_mode == GMP_RNDU && sign_quotient == 1)
	|| (rnd_mode == GMP_RNDN))
      {	  
	/* We cannot round, so that the last bits of the quotient
	   have to be zero; just look if the remainder is nonzero */
	k = rsize - 1; 
	while (k >= 0) { if (tp[k--]) break; }
	if (k >= 0) /* non-zero remainder */
	  cc = mpn_add_1(rp, rp, rrsize, 1 << (BITS_PER_MP_LIMB - 
					       (prec & 
						(BITS_PER_MP_LIMB - 1))));
      }

  if (sign_quotient == -1) { CHANGE_SIGN(r); } 
  r->_mp_exp = rexp;
  
  if (cc) {
    mpn_rshift(rp, rp, rrsize, 1);
    rp[rrsize-1] |= (mp_limb_t) 1 << (BITS_PER_MP_LIMB-1);
    r->_mp_exp++; 
  }
    
  rp [0] &= ~((1 << (BITS_PER_MP_LIMB - 
		     (prec & (BITS_PER_MP_LIMB - 1)))) - 1) ; 
  
  rsize = rrsize; 
  rrsize = (prec - 1)/BITS_PER_MP_LIMB + 1;  
  MPN_COPY(r->_mp_d, rp + rsize - rrsize, rrsize); 
  TMP_FREE (marker);
}

