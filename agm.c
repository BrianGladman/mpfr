#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "mpfr.h"


void mpfr_agm(mpfr_ptr r, mpfr_srcptr a, mpfr_srcptr b, unsigned char rnd_mode)
{
  int i, n, p, q, go_on;
  double uo, vo, rap;
  mpfr_t u, v, tmp, zero,prec;


  /* I want b>= a */
  if (mpfr_cmp(a,b) > 0) 
    return mpfr_agm(r, b, a, rnd_mode);


  /* If a or b is NaN, let's return NaN */
  if (FLAG_NAN(a) || FLAG_NAN(b)) 
    { SET_NAN(r); return; }


  /* If a or b is negative, let's return NaN */
  if ((SIGN(a)<0)||(SIGN(b)<0))
    { SET_NAN(r); return; }

  mpfr_init2(zero,PREC(r));
  mpfr_set_si(zero,0,GMP_RNDZ);

  
  /* If a or b is 0, let's return 0 set to the precision of the result */
  if ((SIGN(a)==0)||(SIGN(b)==0)) 
    {  mpfr_set(r,zero,GMP_RNDZ);
    return;
    }

 /* precision of the following calculus */
  q = PREC(r);
  p = q + 12;


  /* Initialisations */
  go_on=1;

  vo=mpfr_get_d(b);
  uo=mpfr_get_d(a);

  rap=ceil(log((vo-uo)/uo)/log(2));

  /* Main loop */

  while (go_on==1) {
    int can_round, err;
    mpfr_t tmpu, tmpv;

    mpfr_init2(u,p);
    mpfr_init2(v,p);

    mpfr_set(u,a,GMP_RNDZ);
    mpfr_set(v,b,GMP_RNDU);

    mpfr_init2(tmpu,p);
    mpfr_init2(tmpv,p);
    mpfr_init2(tmp,p);

    n = ceil(log(p)/log(2));
    err=ceil(log((3*n+2)*exp(-p*log(2))+3*exp(p*uo*log(2)/(vo-uo))/log(2)));
    
    for(i=0;i<n;i++) {
      mpfr_mul(tmp,u,v,GMP_RNDZ);
      mpfr_sqrt(tmpu,tmp,GMP_RNDZ);
      mpfr_add(tmp,u,v,GMP_RNDU);
      mpfr_div_2exp(tmpv,tmp,1,GMP_RNDU);
    
      mpfr_set(u,tmpu,GMP_RNDZ);
      mpfr_set(v,tmpv,GMP_RNDU);
    }
    /*printf("avant can_round(");
      mpfr_out_str(stdout,10,0,v,GMP_RNDU);
      mpfr_print_raw(v);
      printf(", %i, GMP_RNDU, %i, %i)\n",p-err,rnd_mode,q);*/

    can_round=mpfr_can_round(v,p-err,GMP_RNDU,rnd_mode,q);

    /*printf("apres can_round : %i \n",can_round);*/
    if (can_round==1) {
      go_on=0;
    }
    else {
      if (rap <=1)
	printf("Erreur bizarre\n");
      else {
	p=ceil(p+rap);
	printf("Avec plus de chiffres calculer tu dois\n");
      	mpfr_clear(u);
	mpfr_clear(v);
      }
    }
    mpfr_clear(tmpu);
    mpfr_clear(tmpv);
    mpfr_clear(tmp);
  }
  /* Setting of the result */
 
  mpfr_set(r,v,rnd_mode);
 

  /* Exactness of the result */
  mpfr_init2(prec,p);
  mpfr_init2(tmp,p);
  
  mpfr_sub(prec,v,u,GMP_RNDU);
  mpfr_div(tmp,prec,u,GMP_RNDU);
  mpfr_div_2exp(prec,tmp,p,GMP_RNDU);

  if( (int) floor(mpfr_get_d(prec)) !=0)
    printf("Plus d'iterations effectuer tu dois : %i ulps de diff\n", (int) floor(mpfr_get_d(prec)));

  /*printf("entre u et v : %i ulp\n",(int) floor(mpfr_get_d(prec))); */
  


  /* Let's clean */
  mpfr_clear(u);
  mpfr_clear(v);
  mpfr_clear(tmp);
  mpfr_clear(zero);
  mpfr_clear(prec);

  return ;
}
