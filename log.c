#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "mpfr.h"


void mpfr_log(mpfr_ptr r, mpfr_srcptr a, unsigned char rnd_mode) {
  int p, m, q, bool;
  mpfr_t cst, rapport, agm, tmp1, tmp2, s, mm;
  double x, ref;

  /* If a is NaN, let's return NaN */
  if (FLAG_NAN(a)) 
    { SET_NAN(r); return; }

  /* If a is negative or null, let's return NaN */
  if (SIGN(a)<=0)
    { SET_NAN(r); return; }

  q=PREC(r);
  x=mpfr_get_d(a);
  
  p=q+5;
  bool=1;

  while (bool==1) {
    ref=exp(((double) p) *log(2)/2);
    m=q/2 - EXP(a) - 2;
    while (x<=ref) {
      m++;
      x*=2;
    }

  /*  printf("m : %i\n",m); */
  mpfr_init2(cst,p);
  mpfr_init2(rapport,p);
  mpfr_init2(agm,p);
  mpfr_init2(tmp1,p);
  mpfr_init2(tmp2,p);
  mpfr_init2(s,p);
  mpfr_init2(mm,p);
  
  mpfr_set_si(mm,-m,rnd_mode);

  mpfr_set_si(tmp1,1,rnd_mode);
  /*printf("1 : \t");
    mpfr_out_str(stdout,10,0,tmp1,GMP_RNDN);printf("\n"); */

  mpfr_set_si(tmp2,4,rnd_mode);
  /* printf("4 : \t");
     mpfr_out_str(stdout,10,0,tmp2,GMP_RNDN);printf("\n"); */

  mpfr_mul_2exp(s,a,m,rnd_mode);
  /*printf("s : \t");
    mpfr_out_str(stdout,10,0,s,GMP_RNDN);printf("\n");*/

  /*printf("div entre ");mpfr_out_str(stdout,10,0,tmp2,GMP_RNDN);printf(" et ");
  mpfr_out_str(stdout,10,0,s,GMP_RNDN);printf("\n"); 
  printf("div entre \n");
  mpfr_print_raw(tmp2);printf("    et \n");
  mpfr_print_raw(s);  printf(" en arrondi \n"); */

  mpfr_div(rapport,tmp2,s,rnd_mode);
  /*  printf("4/s : \t");
      mpfr_out_str(stdout,10,0,rapport,GMP_RNDN);printf("\n"); */
 
  mpfr_agm(agm,tmp1,rapport,rnd_mode);
  /*printf("AG : \t");
    mpfr_out_str(stdout,10,0,agm,GMP_RNDN);printf("\n");*/
 
  mpfr_mul_2exp(tmp1,agm,1,rnd_mode);
  /*printf("2AG : \t");
    mpfr_out_str(stdout,10,0,tmp1,GMP_RNDN);printf("\n");*/
 
  mpfr_pi(cst, rnd_mode);
  /*printf("pi : \t");
    mpfr_out_str(stdout,10,0,cst,GMP_RNDN);printf("\n");*/
  
  mpfr_div(tmp2,cst,tmp1,rnd_mode);
  /*printf("pi/2AG : \t");
    mpfr_out_str(stdout,10,0,tmp2,GMP_RNDN);printf("\n");*/
 
  mpfr_log2(cst,rnd_mode);
  /*printf("log2 : \t");
    mpfr_out_str(stdout,10,0,cst,GMP_RNDN);printf("\n");*/
 
  mpfr_mul(tmp1,cst,mm,rnd_mode);
  /*printf("-mlog2 : \t");
    mpfr_out_str(stdout,10,0,tmp1,GMP_RNDN);printf("\n");*/
 
  mpfr_add(cst,tmp1,tmp2,rnd_mode);
  /*printf("res : \t");
    mpfr_out_str(stdout,10,0,cst,GMP_RNDN);printf("\n");*/
 
 
  if(mpfr_can_round(cst,p-4,rnd_mode,rnd_mode,q)==1) {
    mpfr_set(r,cst,rnd_mode);
    bool=0;
  }
  else {
    printf("Avec plus de precisions calculer tu dois !!\n");
    p+=5;
  }

  mpfr_clear(rapport);
  mpfr_clear(agm);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  mpfr_clear(s);
  mpfr_clear(mm);
  mpfr_clear(cst);
  }
 
}


