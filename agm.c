#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"


void mpfr_agm(mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b, unsigned char rnd_mode)
{
  int i, n, p, q, go_on, dn, dp;
  double uo, vo;
  mpfr_t u, v, tmp;

  /* I want b>= a */
  if (mpfr_cmp(a,b) > 0) 
    return mpfr_agm(r, b, a, rnd_mode);


  /* printf("debut\n");*/
  /* If a or b is NaN, let's return NaN */
  if (FLAG_NAN(a) || FLAG_NAN(b)) 
    { SET_NAN(r); return; }


  /* If a or b is negative, let's return NaN */
  if ((SIGN(a)<0)||(SIGN(b)<0))
    { SET_NAN(r); return; }


  
  /* If a or b is 0, let's return 0 as a result */
  if ((SIGN(a)==0)||(SIGN(b)==0)) 
    { SET_ZERO(r);
    return;
    }

 /* precision of the following calculus */
  q = PREC(r);
  p = q + 12;


  /* Initialisations */
  go_on=1;

  vo=mpfr_get_d(b);
  uo=mpfr_get_d(a);
  n = ceil(log(p)/log(2)) +1;
  dn=0;
  dp=0;

  /* Main loop */

  while (go_on==1) {
    int can_go_on, err;
    mpfr_t tmpu, tmpv;

    n += dn;
    p +=dp;
    err=ceil((3*n+2)*exp(-p*log(2))+3*exp(-p*uo*log(2)/(vo-uo))/log(2));
    

    /*    while (p<=err) {
      p=err+1;
      n = ceil(log(p)/log(2)) +1;
      err=ceil(log((3*n+2)*exp(-p*log(2))+3*exp(-p*uo*log(2)/(vo-uo))/log(2)));
      printf("s");
      }*/

    mpfr_init2(u,p);
    mpfr_init2(v,p);

    mpfr_set(u,a,GMP_RNDZ);
    mpfr_set(v,b,GMP_RNDU);

    mpfr_init2(tmpu,p);
    mpfr_init2(tmpv,p);
    mpfr_init2(tmp,p);

    

    /*    printf("p : %i et err : %i\n",p,err);

	  printf("internal loop\n"); */
    for(i=0;i<n;i++) { 
      /*  printf("\n v :");    
      mpfr_out_str(stdout,10,0,v,GMP_RNDN); printf("\n u :");
      mpfr_out_str(stdout,10,0,u,GMP_RNDN);printf("\n u*v :"); */

      mpfr_mul(tmp,u,v,GMP_RNDN);
      /* mpfr_out_str(stdout,10,0,tmp,GMP_RNDN);printf("\n sqrt(u*v) :");*/
      mpfr_sqrt(tmpu,tmp,GMP_RNDN);
      /*mpfr_out_str(stdout,10,0,tmpu,GMP_RNDN);printf("\n"); */
      mpfr_add(tmp,u,v,GMP_RNDN);
      mpfr_div_2exp(tmpv,tmp,1,GMP_RNDN);
      mpfr_set(u,tmpu,GMP_RNDN);
      mpfr_set(v,tmpv,GMP_RNDN);       
    }

    /* printf("avant can_round\n v :"); 
    mpfr_out_str(stdout,10,0,v,GMP_RNDN); printf("\n u :");
    mpfr_out_str(stdout,10,0,u,GMP_RNDN);printf("\n"); */
    /*   mpfr_print_raw(v);
	 printf(", %i, GMP_RNDU, %i, %i)\n",p-err,rnd_mode,q); */


    /* Exactness of the result */
    mpfr_sub(tmp,v,u,GMP_RNDU);
    /*mpfr_out_str(stdout,10,0,tmp,GMP_RNDN); printf(" as sub\n");*/
    mpfr_div(tmpu,tmp,u,GMP_RNDU);
    mpfr_mul_2exp(tmp,tmpu,q,GMP_RNDU);

    if( ((int) floor(mpfr_get_d(tmp))) >=3) {
      /*     printf("Plus d'iterations effectuer tu dois : %i ulps de diff\n", (int) floor(mpfr_get_d(tmp))); */
      go_on=1;
      dn=5;
      dp=0;
    }

    else {
      go_on=1-mpfr_can_round(v,p-err,GMP_RNDN,rnd_mode,q);
      dp=5;
      dn=0;
      /* printf("Avec plus de chiffres calculer tu dois\n");*/
    }

    if(go_on==1) {
      mpfr_clear(u);
      mpfr_clear(v);
    }

    mpfr_clear(tmpu);
    mpfr_clear(tmpv);
    mpfr_clear(tmp);

  }

  /* End of while */
  /* Setting of the result */
 
  mpfr_set(r,v,rnd_mode);


  /* Let's clean */
  mpfr_clear(u);
  mpfr_clear(v);
  mpfr_clear(tmp);

  return ;
}
