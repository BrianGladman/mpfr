
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"


#undef A
#undef B
#define C
#define C1  3
#define C2  2
#define GENERIC mpfr_sin_aux 
#include "generic.c" 
#undef C
#undef C1
#undef C2
#undef GENERIC


#undef A
#undef B
#define C
#define C1  1
#define C2  2
#define GENERIC mpfr_cos_aux 
#include "generic.c" 



int mpfr_extract(mpz_ptr , mpfr_srcptr , int );

#define shift (BITS_PER_MP_LIMB / 2)

int
#if __STDC__
mpfr_sin_cos(mpfr_ptr sinus, mpfr_ptr cosinus, mpfr_srcptr x, mp_rnd_t rnd_mode)
#else
mpfr_sin_cos(sinus, cosinus,x,rnd_mode)
mpfr_ptr sinus;
mpfr_ptr cosinus;
mpfr_srcptr x; 
mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t_sin, t_cos;
  mpfr_t x_copy;
  int i,k;
  mpz_t uk;
  mpz_t square;
  mpfr_t tmp_sin, tmp_cos;
  mpfr_t tmp;
  mpfr_t inter;
  int ttt;
  int twopoweri;
  int Prec;
  int loop;
  int prec_x;
  int shift_x = 0;
  int good = 0;
  int realprec = 0;
  int iter;
  int factor = 2;
  int logn;
  int tmp_factor;
  int tmpi;
  if (FLAG_NAN(x)) { SET_NAN(sinus);  SET_NAN(cosinus); return 1; }
  if (!NOTZERO(x)) { 
    mpfr_set_ui(sinus, 0, GMP_RNDN); 
    mpfr_set_ui(cosinus, 1, GMP_RNDN); 
    return 0;
 }

  prec_x = (int) ceil(log
		      ((double) (PREC(x)) / (double) BITS_PER_MP_LIMB)
		      /log(2.0));  
  logn =  (int) ceil(log
		      ((double) prec_x)
		      /log(2.0));  
  if (logn < 2) logn = 2;
  ttt = EXP(x);
  mpfr_init2(x_copy,PREC(x));
  mpfr_set(x_copy,x,GMP_RNDD);
  mpz_init(square);
  /* on fait le shift pour que le nombre soit inferieur a 1 */
  if (ttt > 0) 
    {
      shift_x = ttt;
      mpfr_mul_2exp(x_copy,x,-ttt, GMP_RNDN); 
      ttt = EXP(x_copy);
    }
  realprec = PREC(sinus)+logn;
  while (!good){
  try_again:
    Prec = realprec + 2*shift + 2 + shift_x + factor;
    k = (int) ceil(log
		   ((double) (Prec) / (double) BITS_PER_MP_LIMB)
		   /log(2.0));
    /* Maintenant, il faut extraire : */
    mpfr_init2(t_sin, Prec);
    mpfr_init2(t_cos, Prec);
    mpfr_init2(tmp, Prec);
    mpfr_init2(tmp_sin, Prec);
    mpfr_init2(tmp_cos, Prec);
    mpfr_init2(inter, Prec);
    mpfr_set_ui(tmp_sin,0,GMP_RNDN);
    mpfr_set_ui(tmp_cos,1,GMP_RNDN);
    twopoweri = BITS_PER_MP_LIMB;   
    if (k <= prec_x) iter = k; else iter= prec_x;
    for(i = 0; i <= iter; i++){
      mpfr_extract(uk,x_copy,i);
	if (i)
	  {
	    mpz_set(square,uk);
	    mpz_mul(square, square, square);
	    mpz_neg(square, square);
	    mpfr_sin_aux(t_sin,square, 2*(twopoweri - ttt) + 2, k - i + 1);
	    mpfr_cos_aux(t_cos,square, 2*(twopoweri - ttt) + 2, k - i + 1);

	    mpfr_set_z(tmp,uk,GMP_RNDD);
	    mpfr_mul(t_sin, t_sin, tmp,GMP_RNDD);
	    mpfr_div_2exp(t_sin,t_sin, twopoweri - ttt, GMP_RNDD);
	  }
	else
	  {
	    /* cas particulier : on est oblige de faire les calculs avec x/2^. 
	       puis elever au carre (plus rapide) */    
	    mpz_set(square,uk);
	    mpz_mul(square, square, square);
	    mpz_neg(square, square);
	    /* pour l'instant, shift = 0 ... */
	    /* ATTENTION, IL FAUT AUSSI MULTIPLIER LE DENOMINATEUR */
	    mpfr_sin_aux(t_sin,square, 2*(shift + twopoweri - ttt) + 2, k+1);
	    mpfr_cos_aux(t_cos,square, 2*(shift + twopoweri - ttt) + 2, k+1);
	    mpfr_set_z(tmp,uk,GMP_RNDD);
	    mpfr_mul(t_sin, t_sin, tmp,GMP_RNDD);
	    /* LA AUSSI, IL FAUT PENSER A DECALER DE twopoweri - ttt) */
	    mpfr_div_2exp(t_sin,t_sin, twopoweri - ttt + shift, GMP_RNDD);
       	    for (loop= 0 ; loop < shift; loop++){
	      /* t_sin = sin(a)
		 t_cos = cos(a) */
	      /* on veut t_sin = 2 sin a cos a
		 et t_cos = 2 * cos^2 - 1 */
	      mpfr_mul(t_sin, t_sin, t_cos, GMP_RNDD);
	      mpfr_mul_2exp(t_sin, t_sin, 1, GMP_RNDD);
	      
	      mpfr_mul(t_cos, t_cos, t_cos, GMP_RNDD);
	      mpfr_mul_2exp(t_cos, t_cos, 1, GMP_RNDD);
      	      mpfr_sub_ui(t_cos, t_cos, 1,  GMP_RNDD);
	    }
	  }
	mpz_clear(uk);
	/* on utilise cos(a+b) = cos a cos b - sin a sin b 
	   sin(a+b) = sin a cos b + cos a sin b */
	/* Donnees : 
	   tmp = cos(a) 
	   tmp_sin = sin(a)
	   t_sin = sin(b)
	   t_cos = cos(b) */	  
	mpfr_set(tmp, tmp_cos,GMP_RNDD);
	/* inter = sin a sin b */
	mpfr_mul(inter, tmp_sin, t_sin, GMP_RNDD);
	/* tmp_cos = cos a cos b */
	mpfr_mul(tmp_cos, tmp_cos, t_cos, GMP_RNDD);
	/* tmp_cos = cos (a+b) */
	mpfr_sub(tmp_cos, tmp_cos, inter, GMP_RNDD);
	/* inter = cos a sin b */
	mpfr_mul(inter, tmp, t_sin, GMP_RNDD);
	/* tmp_sin = sin a cos b */
	mpfr_mul(tmp_sin, tmp_sin, t_cos, GMP_RNDD);
	/* tmp_sin = sin (a+b) */
	mpfr_add(tmp_sin, tmp_sin, inter, GMP_RNDD);
	twopoweri <<= 1;
    }
    tmp_factor= factor;
    for (loop= 0 ; loop < shift_x; loop++){
      mpfr_mul(tmp_sin, tmp_sin, tmp_cos, GMP_RNDD);
      mpfr_mul_2exp(tmp_sin, tmp_sin, 1, GMP_RNDD);
      mpfr_mul(tmp_cos, tmp_cos, tmp_cos, GMP_RNDD);
      mpfr_mul_2exp(tmp_cos, tmp_cos, 1, GMP_RNDD);
      mpfr_set_ui(tmp, 1,GMP_RNDN);      
      tmpi = -EXP(tmp_cos);
      mpfr_sub(tmp_cos, tmp_cos, tmp,  GMP_RNDD);
      /* rep\'erer si le nombre de chiffres obtenu est suffisant pour 
	 avoir la bonne pr\'ecision. Le probl\`eme : comment faire ? 
         la pr\'ecision s'obtient en comparant 
	 (Prec-factor) a la pr\'ecision obtenue r\'eellement, celle-ci
	 \'etant donn\'ee par Prec + EXP(tmp_cos). 
	 il faut donc comparer EXP(tmp_cos) a factor */
      tmp_factor -= -EXP(tmp_cos) + tmpi;    
      if (tmp_factor <= 0)
	{
	  factor += -tmp_factor  + 5;
	  mpfr_clear(t_sin);
	  mpfr_clear(t_cos);
	  mpfr_clear(tmp);
	  mpfr_clear(tmp_sin);
	  mpfr_clear(inter);
	  goto try_again;
	}      
    }
    if (mpfr_can_round(tmp_sin, realprec, GMP_RNDD, rnd_mode, PREC(sinus))){
	mpfr_set(sinus,tmp_sin,rnd_mode);
	mpfr_set(cosinus,tmp_cos,rnd_mode);
	good = 1;
    } else {
      mpfr_clear(t_sin);
      mpfr_clear(t_cos);
      mpfr_clear(tmp);
      mpfr_clear(tmp_sin);
      mpfr_clear(inter);
      realprec += 3*logn;
    }
    }
  mpz_clear(square);
  return 0;
} 



