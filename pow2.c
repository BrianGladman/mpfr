/* mpfr_pow -- power function x^y 

Copyright (C) 2001 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

 /* The computation of y=pow(x,z) is done by

    y=exp(z*log(x))=x^z
 */

/* check if the mpfr input is an integer : 1 input is an integer, -1 not*/
int mpfr_pow_si _PROTO ((mpfr_ptr, mpfr_srcptr, long int, mp_rnd_t)); 

int
#if __STDC__
mpfr_pow_si (mpfr_ptr y, mpfr_srcptr x, long int n, mp_rnd_t rnd_mode)
#else
mpfr_pow_si (y, x, n, rnd)
     mpfr_ptr y;
     mpfr_srcptr x;
     long int n;
     mp_rnd_t rnd_mode;
#endif
{

  if (n>0)
      return mpfr_pow_ui(y,x,(unsigned long int)n,rnd_mode);
  else
    {

      int inexact = 0;

      n=-n;

      /* x is NaN*/
      if (MPFR_IS_NAN(x)) 
        {
          MPFR_SET_NAN(y); 
          return 1; 
        }
      MPFR_CLEAR_NAN(y);

      /* n=0 */
      if(n==0)
        return mpfr_set_ui(y,1,GMP_RNDN);;

      /* case x is INF */
      if(MPFR_IS_INF(x))
        {
          if(MPFR_SIGN(x)>0) /* +Inf */
            {
              MPFR_SET_ZERO(y);
              if(MPFR_SIGN(y) < 0)
                MPFR_CHANGE_SIGN(y);
              return 0;
            }
          else
            {
              MPFR_SET_ZERO(y); /* -Inf */
              if(!(n%2))        /* n is odd */
                {
                  if(MPFR_SIGN(y) > 0)
                    MPFR_CHANGE_SIGN(y);
                }
              else              /* n is not odd */
                {
                  if(MPFR_SIGN(y) < 0)
                    MPFR_CHANGE_SIGN(y);
                }
              return 0;
            }
        }

      /* case x=0 */
      if(mpfr_cmp_ui(x,0) == 0)
        {
          if(!(n%2))             /* n is odd */
           {
             MPFR_SET_INF(y);
             MPFR_SET_SAME_SIGN(y,x);
             DIVIDE_BY_ZERO;    /* Execption GMP*/
             return 0;
           }
         else                   /* n is not odd */
           {
             MPFR_SET_INF(y);
             if(MPFR_SIGN(y) < 0)
               MPFR_CHANGE_SIGN(y);
             DIVIDE_BY_ZERO;    /* Execption GMP*/
             return 0;
           }
        }

      MPFR_CLEAR_INF(y);

      /* General case */
      {
        /* Declaration of the intermediary variable */
        mpfr_t t, ti;

        /* Declaration of the size variable */
        mp_prec_t Nx = MPFR_PREC(x);   /* Precision of input variable */
        mp_prec_t Ny = MPFR_PREC(y);   /* Precision of input variable */

        mp_prec_t Nt;   /* Precision of the intermediary variable */
        long int err;  /* Precision of error */
                
        /* compute the precision of intermediary variable */
        Nt=MAX(Nx,Ny);
        /* the optimal number of bits : see algorithms.ps */
        Nt=Nt+3+_mpfr_ceil_log2(Nt);

        /* initialise of intermediary	variable */
        mpfr_init(t);
        mpfr_init(ti);

        do {

          /* reactualisation of the precision */
          mpfr_set_prec(t,Nt);                    
          mpfr_set_prec(ti,Nt);             

          /* compute 1/(x^n) n>0*/
          mpfr_pow_ui(ti,y,(unsigned long int)(n),GMP_RNDN);
          mpfr_ui_div(t,1,ti,GMP_RNDN);

          /* estimation of the error -- see pow function in algorithms.ps*/
          err = Nt - 3;

          /* actualisation of the precision */
          Nt += 10;

        } while (err<0 || !mpfr_can_round(t,err,GMP_RNDN,rnd_mode,Ny));
 
        inexact = mpfr_set(y,t,rnd_mode);
        mpfr_clear(t);
        mpfr_clear(ti);
      }
      return inexact;
    }
}

int mpfr_isinteger _PROTO((mpfr_srcptr));

int 
#if __STDC__
mpfr_isinteger(mpfr_srcptr x)
#else
mpfr_isinteger(x)
     mpfr_srcptr x;
#endif
{

  mpfr_t u;
  int expo;
  mp_prec_t prec;

  expo=(int)MPFR_EXP(x);
  prec=MPFR_PREC(x);

  if (expo<=0) 
    return 0;

  if (expo>=prec) 
    return 1;

  mpfr_init2(u,prec);
  mpfr_trunc(u,x);

  if(mpfr_cmp(x,u)==0) return 1;
  else return 0;
}

int mpfr_pow _PROTO ((mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mp_rnd_t));

int
#if __STDC__
mpfr_pow (mpfr_ptr z, mpfr_srcptr x ,mpfr_srcptr y , mp_rnd_t rnd_mode) 
#else
mpfr_pow (z, x, y, rnd_mode)
     mpfr_ptr z;
     mpfr_srcptr x;
     mpfr_srcptr y;
     mp_rnd_t rnd_mode;
#endif
{
  int inexact = 0;
 
  if (MPFR_IS_NAN(x) || MPFR_IS_NAN(y) ) 
    {
      MPFR_SET_NAN(z); 
      return 1; 
    }

  if (MPFR_IS_INF(y))
    {
      mpfr_t px;
      mpfr_init2(px,MPFR_PREC(x));
      mpfr_abs(px,x,GMP_RNDN);
      if(MPFR_SIGN(y)>0)
        {
          if(MPFR_IS_INF(x))
          {
            if(MPFR_SIGN(x)>0)
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_INF(z);
              if(MPFR_SIGN(z) <0)
                MPFR_CHANGE_SIGN(z);
              mpfr_clear(px);
              return 0;
            }
            else
            {
              MPFR_CLEAR_FLAGS(z);  
              MPFR_SET_ZERO(z);
              if(MPFR_SIGN(z) <0)
                MPFR_CHANGE_SIGN(z);
              mpfr_clear(px);
              return 0;
            }
          }
          if(mpfr_cmp_ui(px,1) > 0)
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_INF(z);
              if(MPFR_SIGN(z) <0)
                MPFR_CHANGE_SIGN(z);
              mpfr_clear(px);
              return 0;
            }
          if(mpfr_cmp_ui(px,1) < 0)
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_ZERO(z);
              if(MPFR_SIGN(z) <0)
                MPFR_CHANGE_SIGN(z);
              mpfr_clear(px);
              return 0;
            }
          if(mpfr_cmp_ui(px,1)==0)
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_NAN(z);
              mpfr_clear(px);
              return 1;
            }
        }
      else
        {
          if(MPFR_IS_INF(x))
          {
            if(MPFR_SIGN(x)>0)
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_ZERO(z);
              if(MPFR_SIGN(z) <0)
                MPFR_CHANGE_SIGN(z);
              mpfr_clear(px);
              return 0;
            }
            else
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_INF(z);
              if(MPFR_SIGN(z) <0)
                MPFR_CHANGE_SIGN(z);
              mpfr_clear(px);
              return 0;
            }
          }
          if(mpfr_cmp_ui(px,1) > 0)
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_ZERO(z);
              if(MPFR_SIGN(z) <0)
                MPFR_CHANGE_SIGN(z);
              mpfr_clear(px);
              return 0;
            }
          if(mpfr_cmp_ui(px,1) < 0)
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_INF(z);
              if(MPFR_SIGN(z) <0)
                MPFR_CHANGE_SIGN(z);
              mpfr_clear(px);
              return 0;
            }
          if(mpfr_cmp_ui(px,1)==0)
            {
              MPFR_CLEAR_FLAGS(z);
              MPFR_SET_NAN(z);
              mpfr_clear(px);
              return 1;
            }
        }
    }

  if(MPFR_IS_ZERO(y))
    {
      return mpfr_set_ui(z,1,GMP_RNDN);
    }

  if(mpfr_isinteger(y))
    {
      mpz_t zi;
      long int zii;
      int exptol;
    
      mpz_init(zi);  
      exptol=mpz_set_fr(zi,y);     
        
      if (exptol>0)
        mpz_mul_2exp(zi, zi, exptol);
      else
        mpz_tdiv_q_2exp(zi, zi, (unsigned long int) (-exptol));

      zii=mpz_get_ui(zi);
        
      mpz_clear(zi);
      return mpfr_pow_si(z,x,zii,rnd_mode); 
    }
  if (MPFR_IS_INF(x))
    {
      if (MPFR_SIGN(x) > 0)
        {
        if (MPFR_SIGN(y) >0)
          {
            MPFR_CLEAR_FLAGS(z);
            MPFR_SET_INF(z);
            if(MPFR_SIGN(z) <0)
              MPFR_CHANGE_SIGN(z);
            return 0;
          }
        else
          {
            MPFR_CLEAR_FLAGS(z);
            MPFR_SET_ZERO(z);
            if(MPFR_SIGN(z) <0)
              MPFR_CHANGE_SIGN(z);
            return 0;
          }
        }
      else
        {
          MPFR_CLEAR_FLAGS(z);
          MPFR_SET_NAN(z); 
          return 1; 
        }
    }       
    
  MPFR_CLEAR_INF(z);
  if(MPFR_SIGN(x) < 0)
    {
      MPFR_CLEAR_FLAGS(z);
      MPFR_SET_NAN(z); 
      return 1; 
    }
  MPFR_CLEAR_NAN(z);

  if(mpfr_cmp_ui(x,0) == 0)
    {
      MPFR_CLEAR_FLAGS(z);
      MPFR_SET_ZERO(z);
      return 0;
    }
  /* General case */
  {
    /* Declaration of the intermediary variable */
      mpfr_t t, te, ti;

      /* Declaration of the size variable */
      mp_prec_t Nx = MPFR_PREC(x);   /* Precision of input variable */
      mp_prec_t Ny = MPFR_PREC(y);   /* Precision of input variable */

      mp_prec_t Nt;   /* Precision of the intermediary variable */
      long int err;  /* Precision of error */
                
      /* compute the precision of intermediary variable */
      Nt=MAX(Nx,Ny);
      /* the optimal number of bits : see algorithms.ps */
      Nt=Nt+5+_mpfr_ceil_log2(Nt);

      /* initialise of intermediary	variable */
      mpfr_init(t);
      mpfr_init(ti);
      mpfr_init(te);             

      do {

	/* reactualisation of the precision */
	mpfr_set_prec(t,Nt);                    
	mpfr_set_prec(ti,Nt);                    
   	mpfr_set_prec(te,Nt);             

	/* compute   exp(y*ln(x))*/
        mpfr_log(ti,x,GMP_RNDU);         /* ln(n) */
        mpfr_mul(te,y,ti,GMP_RNDU);       /* y*ln(n) */
        mpfr_exp(t,te,GMP_RNDN);         /* exp(x*ln(n))*/

	/* estimation of the error -- see pow function in algorithms.ps*/
	err = Nt - _mpfr_ceil_log2(1+pow(2,MPFR_EXP(te)+2));

	/* actualisation of the precision */
	Nt += 10;

      } while (err<0 || !mpfr_can_round(t,err,GMP_RNDN,rnd_mode,Ny));
      inexact = mpfr_set(z,t,rnd_mode);
      mpfr_clear(t);
      mpfr_clear(ti);
      mpfr_clear(te);
    }
    return inexact;
}







