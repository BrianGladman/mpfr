#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "mpfr.h"


void mpfr_agm(mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b, unsigned char rnd_mode)
{
  int i, n, p, q;
  mpfr_t u, v, tmpu, tmpv, tmp, zero,prec;


  /* If a or b is NaN, let's return NaN */
  if (FLAG_NAN(a) || FLAG_NAN(b)) 
    { SET_NAN(r); return; }


  /* If a or b is negative, let's return NaN */
  if ((SIGN(a)<0)||(SIGN(b)<0))
    { SET_NAN(r); return; }

  mpfr_init(zero);
  mpfr_set_si(zero,0,GMP_RNDZ);

  
  /* If a or b is 0, let's return 0 set to the precision of the result */
  if ((mpfr_cmp(a,zero)==0)||(mpfr_cmp(b,zero)==0)) 
    {  mpfr_set(r,zero,GMP_RNDZ);
    return;
    }

 /* precision of the following calculus */
  q = PREC(r);
  p = q + 12;
 
  n = floor(log(p)/log(2));


  /* Initialisations */
  mpfr_init2(u,p);
  mpfr_init2(v,p);

  if (mpfr_cmp(a,b) >= 0) { 
    mpfr_set(u,b,GMP_RNDZ);
    mpfr_set(v,a,GMP_RNDU);
  }
  else {
    mpfr_set(u,a,GMP_RNDZ);
    mpfr_set(v,b,GMP_RNDU);
  }

  mpfr_init2(tmp,p);
  mpfr_init2(tmpu,p);
  mpfr_init2(tmpv,p);


  /* Main loop */
  for(i=0;i<n;i++) {
    mpfr_mul(tmp,u,v,GMP_RNDZ);
    mpfr_sqrt(tmpu,tmp,GMP_RNDZ);
    mpfr_add(tmp,u,v,GMP_RNDU);
    mpfr_div_2exp(tmpv,tmp,1,GMP_RNDU);
    
    mpfr_set(u,tmpu,GMP_RNDZ);
    mpfr_set(v,tmpv,GMP_RNDU);
  }

  /* Setting of the result */
  mpfr_set(r,v,GMP_RNDN);

  /* Exactness of the result */
  mpfr_init2(prec,p);
  mpfr_sub(prec,v,u,GMP_RNDU);
  mpfr_div(tmp,prec,u,GMP_RNDU);
  mpfr_div_2exp(prec,tmp,q,GMP_RNDU);

  /* printf("entre u et v : %i ulp\n",(int) floor(mpfr_get_d(prec))); */
  


  /* Let's clean */
  mpfr_clear(u);
  mpfr_clear(v);
  mpfr_clear(tmpu);
  mpfr_clear(tmpv);
  mpfr_clear(tmp);
  mpfr_clear(zero);
  mpfr_clear(prec);

  return ;
}
