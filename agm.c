#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"



#define MON_INIT(xp, x, p, s) xp = (mp_ptr) TMP_ALLOC(s*BYTES_PER_MP_LIMB);    x -> _mp_prec = p; x -> _mp_d = xp; x -> _mp_size = s; x -> _mp_exp = 0; 




void 
#ifdef __STDC__
mpfr_agm(mpfr_ptr r, mpfr_srcptr op2, mpfr_srcptr op1, unsigned char rnd_mode)
#else
mpfr_agm(r, a, b, rnd_mode)
     mpfr_ptr r; 
     mpfr_srcptr a; 
     mpfr_srcptr b; 
     unsigned char rnd_mode; 
#endif
{
  int i, s, ntotal, p, q, go_on, no, ulps;
  double uo, vo;
  mp_limb_t *up, *vp, *tmpp, *tmpup, *tmpvp, *ap, *bp;
  mpfr_t u, v, tmp, tmpu, tmpv, a, b;
  TMP_DECL(marker1);
  TMP_DECL(marker2);
  TMP_DECL(marker3);


  /* If a or b is NaN, the result is NaN */
  if (FLAG_NAN(op1) || FLAG_NAN(op2)) 
    { SET_NAN(r); return; }


  /* If a or b is negative, the result is NaN */
  if ((SIGN(op1)<0)||(SIGN(op2)<0))
    { SET_NAN(r); return; }


  
  /* If a or b is 0, the result is 0 */
  if ((SIGN(op1)==0)||(SIGN(op2)==0)) 
    { SET_ZERO(r);
    return;
    }

 /* precision of the following calculus */
  q = PREC(r);
  p = q + 15;


  /* Initialisations */
  go_on=1;
  ntotal = ceil(log(p)/log(2)) +1;
  no=0;

  TMP_MARK(marker1);
  s=(p-1)/BITS_PER_MP_LIMB+1;
  MON_INIT(ap, a, p, s);  
  MON_INIT(bp, b, p, s);
  TMP_MARK(marker2);
  MON_INIT(up, u, p, s);
  MON_INIT(vp, v, p, s);   
  MON_INIT(tmpup, tmpu, p, s);  
  MON_INIT(tmpvp, tmpv, p, s);  
  MON_INIT(tmpp, tmp, p, s);  


  /* b and a will be the 2 operands but I want b>= a */
  if (mpfr_cmp(op1,op2) > 0) {
    mpfr_set(b,op1,GMP_RNDN); mpfr_set(a,op2,GMP_RNDN);  
  }
  else {
    mpfr_set(b,op2,GMP_RNDN); mpfr_set(a,op1,GMP_RNDN);  
  }
    
  vo=mpfr_get_d(b);
  uo=mpfr_get_d(a);

  mpfr_set(u,a,GMP_RNDN);
  mpfr_set(v,b,GMP_RNDN);
 
  /* Main loop */

  while (go_on==1) {
    int err;

    err=ceil((3*ntotal+2)*exp(-p*log(2))+3*exp(-p*uo*log(2)/(vo-uo))/log(2));
    
    /* Calculus of un and vn */
    for(i=no;i<ntotal;i++) {    
      mpfr_mul(tmp,u,v,GMP_RNDN);
      /*  printf("dont la racine est\n");   */
      mpfr_sqrt(tmpu,tmp,GMP_RNDN); 
      /*mpfr_print_raw(tmpu);  printf("\n");*/
      mpfr_add(tmp,u,v,GMP_RNDN);
      mpfr_div_2exp(tmpv,tmp,1,GMP_RNDN);
      mpfr_set(u,tmpu,GMP_RNDN);
      mpfr_set(v,tmpv,GMP_RNDN); 
    }

    /*    printf("avant can_round\n v : ");  */
    /* mpfr_out_str(stdout,10,0,v,GMP_RNDN); printf("\n u :");
      mpfr_out_str(stdout,10,0,u,GMP_RNDN);printf("\n"); */
    

    /* Exactness of the result */

    /* Calculus of the nomber of ulps between un and vn */
    /*printf("avant sub entre\n");
     mpfr_print_raw(v);  printf("    -\n");
    mpfr_print_raw(u);  printf("   , en RNDN\n");*/
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
      int round1, round2; 
      mpfr_t res1, res2;
      mp_limb_t *res1p, *res2p;

      TMP_MARK(marker3);
      s=(q-1)/BITS_PER_MP_LIMB+1;
      MON_INIT(res1p, res1, q, s);  
      MON_INIT(res2p, res2, q, s);  
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
	  TMP_FREE(marker2); 
	  TMP_MARK(marker2);
	  s=(p-1)/BITS_PER_MP_LIMB+1;
	  MON_INIT(up, u, p, s);
	  MON_INIT(vp, v, p, s);   
	  MON_INIT(tmpup, tmpu, p, s);  
	  MON_INIT(tmpvp, tmpv, p, s);  
	  MON_INIT(tmpp, tmp, p, s);
	  mpfr_set(u,a,GMP_RNDN);
	  mpfr_set(v,b,GMP_RNDN);
      }
      TMP_FREE(marker3); 
    }
  }
  /* End of while */

  /* Setting of the result */

    mpfr_set(r,v,rnd_mode);


  /* Let's clean */
    TMP_FREE(marker1); 
  
  return ;
}

