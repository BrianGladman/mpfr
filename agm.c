#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"


void mpfr_agm(mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b, unsigned char rnd_mode)
{
  int i, ntotal, p, q, go_on, no, ulps;
  double uo, vo;
  mpfr_t u, v, tmp, tmpu, tmpv;

  /* I want b>= a */
  if (mpfr_cmp(a,b) > 0) 
    return mpfr_agm(r, b, a, rnd_mode);


  /* If a or b is NaN, the result is NaN */
  if (FLAG_NAN(a) || FLAG_NAN(b)) 
    { SET_NAN(r); return; }


  /* If a or b is negative, the result is NaN */
  if ((SIGN(a)<0)||(SIGN(b)<0))
    { SET_NAN(r); return; }


  
  /* If a or b is 0, the result is 0 */
  if ((SIGN(a)==0)||(SIGN(b)==0)) 
    { SET_ZERO(r);
    return;
    }

 /* precision of the following calculus */
  q = PREC(r);
  p = q + 15;


  /* Initialisations */
  go_on=1;
  vo=mpfr_get_d(b);
  uo=mpfr_get_d(a);
  ntotal = ceil(log(p)/log(2)) +1;
  no=0;

  mpfr_init2(u,p);
  mpfr_init2(v,p);
  mpfr_init2(tmpu,p);
  mpfr_init2(tmpv,p);
  mpfr_init2(tmp,p);
  mpfr_set(u,a,GMP_RNDN);
  mpfr_set(v,b,GMP_RNDN);
 

  /* Main loop */

  while (go_on==1) {
    int can_go_on, err;

    err=ceil((3*ntotal+2)*exp(-p*log(2))+3*exp(-p*uo*log(2)/(vo-uo))/log(2));
    
    /* Calculus of un and vn */
    for(i=no;i<ntotal;i++) {    
      mpfr_mul(tmp,u,v,GMP_RNDN);
      mpfr_sqrt(tmpu,tmp,GMP_RNDN);
      mpfr_add(tmp,u,v,GMP_RNDN);
      mpfr_div_2exp(tmpv,tmp,1,GMP_RNDN);
      mpfr_set(u,tmpu,GMP_RNDN);
      mpfr_set(v,tmpv,GMP_RNDN);       
    }

    /*printf("avant can_round\n v :"); 
      mpfr_out_str(stdout,10,0,v,GMP_RNDN); printf("\n u :");
      mpfr_out_str(stdout,10,0,u,GMP_RNDN);printf("\n"); */
    

    /* Exactness of the result */

    /* Calculus of the nomber of ulps between un and vn */
    mpfr_sub(tmp,v,u,GMP_RNDN);
    mpfr_div(tmpu,tmp,u,GMP_RNDN);
    mpfr_mul_2exp(tmp,tmpu,q+1,GMP_RNDN);
    ulps = (int) floor(mpfr_get_d(tmp));
    if (ulps <0) ulps=-ulps;

    /* If there is more than 2 ulps, we have to do more iterations
       with the same precision */
    if( ulps >=2) {
      go_on=1;
      no=ntotal;
      ntotal+=5;
    }
   
    /* Else, we could have to work with more precision */
    else { 
      int round1, round2, equals;
      mpfr_t res1, res2;
      mpfr_init2(res1,q);
      mpfr_init2(res2,q);
      round1=mpfr_can_round(v,p-err-1,GMP_RNDU,rnd_mode,q);
      round2=mpfr_can_round(u,p-err-1,GMP_RNDD,rnd_mode,q);
      mpfr_set(res1, v, rnd_mode); 
      mpfr_set(res2, u, rnd_mode); 

      /* If u=v or u & v can be round to the same machine number, we
	 we have found the correct result */

      if (((ulps==0)&&round1)||((ulps==1)&&round1&&round2&&(mpfr_cmp(res1,res2)==0)))
	go_on=0;

      else {
	  go_on=1;
	  p+=5;
	  no=0;
	  ntotal+=3;
	  mpfr_clear(tmpu);
	  mpfr_clear(tmpv);
	  mpfr_clear(tmp);
	  mpfr_clear(u);
	  mpfr_clear(v);
	  mpfr_init2(u,p);
	  mpfr_init2(v,p);
	  mpfr_init2(tmpu,p);
	  mpfr_init2(tmpv,p);
	  mpfr_init2(tmp,p);
	  mpfr_set(u,a,GMP_RNDN);
	  mpfr_set(v,b,GMP_RNDN);
      }
    }
  }
  /* End of while */

  /* Setting of the result */

    mpfr_set(r,v,rnd_mode);


  /* Let's clean */
  mpfr_clear(u);
  mpfr_clear(v);
  mpfr_clear(tmpu);
  mpfr_clear(tmpv);
  mpfr_clear(tmp);
  
  return ;
}
