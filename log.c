#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"


  /* The computation of log(a) is done using the formula :
     if we want p bits of the result,
                       pi
	  log(a) = ------------  -   m log 2
		    2 AG(1,4/s)


     where s = x 2^m > 2^(p/2)
  */


#define MON_INIT(xp, x, p, s) xp = (mp_ptr) TMP_ALLOC(s*BYTES_PER_MP_LIMB);    x -> _mp_prec = p; x -> _mp_d = xp; x -> _mp_size = s; x -> _mp_exp = 0; 




void 
#if __STDC__
mpfr_log(mpfr_ptr r, mpfr_srcptr a, unsigned char rnd_mode) 
#else
mpfr_log()
     mpfr_ptr r;
     mpfr_srcptr a;
     unsigned char rnd_mode; 
#endif
{
  int p, m, q, bool, err, size;
  mpfr_t cst, rapport, agm, tmp1, tmp2, s, mm;
  mp_limb_t *cstp, *rapportp, *agmp, *tmp1p, *tmp2p, *sp, *mmp;
  double x, ref;
  TMP_DECL(marker);


  /* If a is NaN, the result is NaN */
  if (FLAG_NAN(a)) 
    { SET_NAN(r); return; }

  /* If a is negative or null, the result is NaN */
  if (SIGN(a)<=0)
    { SET_NAN(r); return; }

  /* If a is 1, the result is 0 */
  if (mpfr_cmp_ui_2exp(a,1,0)==0){
    SET_ZERO(r);
    return;
  }


  q=PREC(r);
  
  /* The error due to the lost of bits during the substraction is
     err=lg(p*log(2)/(2*|x-1|)) : greater when x is next to 1 : the
     result is next to 0 but the 2 terms that we substract could be 
     about 10. */

  ref=mpfr_get_d(a)-1.0;
  if (ref<0)
    ref=-ref;
  err=(int) ceil(log(((double) q)*log(2.0)/(2*ref)))+1;
  if (err <0)
    err=1;

  /* The exactness depends on err */
  p=q+11+err;

  bool=1;

  while (bool==1) {
    err=(int) ceil(log(((double) p)*log(2.0)/(2*ref)))+1;
    if (err <0)
      err=1;
    
    /* Calculus of m (depends on p) */
    m=(int) ceil(((double) p)/2.0) -EXP(a)+1;

    /* All the mpfr_t needed have a precision of p */
    TMP_MARK(marker);
    size=(p-1)/BITS_PER_MP_LIMB+1;
    MON_INIT(cstp, cst, p, size);  
    MON_INIT(rapportp, rapport, p, size);
    MON_INIT(agmp, agm, p, size);
    MON_INIT(tmp1p, tmp1, p, size);  
    MON_INIT(tmp2p, tmp2, p, size);  
    MON_INIT(sp, s, p, size);
    MON_INIT(mmp, mm, p, size);

    mpfr_set_si(mm,m,GMP_RNDN);           /* I have m */
    mpfr_set_si(tmp1,1,GMP_RNDN);         /* I have 1 */
    mpfr_set_si(tmp2,4,GMP_RNDN);         /* I have 4 */
    mpfr_mul_2exp(s,a,m,GMP_RNDN);        /* I compute s=a*2^m */ 
    mpfr_div(rapport,tmp2,s,GMP_RNDN);    /* I compute 4/s */
    mpfr_agm(agm,tmp1,rapport,GMP_RNDN);  /* I compute AG(1,4/s) */
    mpfr_mul_2exp(tmp1,agm,1,GMP_RNDN);   /* I compute 2*AG(1,4/s) */
    mpfr_pi(cst, GMP_RNDN);               /* I compute pi */
    mpfr_div(tmp2,cst,tmp1,GMP_RNDN);     /* I compute pi/2*AG(1,4/s) */
    mpfr_log2(cst,GMP_RNDN);              /* I compute log(2) */
    mpfr_mul(tmp1,cst,mm,GMP_RNDN);       /* I compute m*log(2) */
    mpfr_sub(cst,tmp2,tmp1,GMP_RNDN);     /* I compute log(a) */ 
 

    /*printf("avant arrondi : ( %i bits faux)\n",7+err);
      mpfr_print_raw(cst);printf("\n"); */

    
    /* If we can round the result, we set it and go out of the loop */

    if(mpfr_can_round(cst,p-7-err,GMP_RNDN,rnd_mode,q)==1) {
      mpfr_set(r,cst,rnd_mode);
      bool=0;
    }
    /* else we increase the precision */
    else {
      p+=5;
      TMP_FREE(marker);
    }

    /* We clean */
    TMP_FREE(marker);
    
  }
 
}


