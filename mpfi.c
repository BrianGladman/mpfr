/* mpfi.c -- Implementation of mpfi.

Copyright (C) 1999 Free Software Foundation.

This file is part of the MPFR Library.
 
The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.
 
The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.
 
You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include "mpfi.h"

/* *********************************************************** */
/*                   MPFI                                      */
/* *********************************************************** */

#define MPFR_SIGN_PART(x) mpfr_cmp_ui_2exp(x,0,0)

#define MPFI_ISPOS(x) ((MPFR_SIGN_PART((&(x->left)))>=0) && (MPFR_SIGN_PART((&(x->right)))>0))
#define MPFI_MPFR_ISNEG(x) ((MPFR_SIGN_PART((&(x->left)))<0) && (MPFR_SIGN_PART((&(x->right)))<=0))
#define MPFI_ISNULL(x) ((MPFR_SIGN_PART((&(x->left)))==0) && (MPFR_SIGN_PART((&(x->right)))==0))
#define MPFI_HASZERO(x) ((MPFR_SIGN_PART((&(x->left)))<0) && (MPFR_SIGN_PART((&(x->right)))>0))


void mpfi_set_prec(const unsigned int p)
{
  mpfr_set_default_prec(p);
}

int mpfi_cmp_default (mpfi_srcptr a,mpfi_srcptr b)
{
  return((mpfr_cmp(&(a->right),&(b->left))<0) ? -1 :(mpfr_cmp(&(b->right),&(a->left))>0));
}

int (*mpfi_cmp)(mpfi_srcptr,mpfi_srcptr)=mpfi_cmp_default;

int mpfi_comp_d (mpfi_srcptr a,const double b)
{
  mpfi_t tmp;
  mpfi_set_d(tmp,b);
  return(mpfi_cmp(a,tmp));
}

int mpfi_comp_ui (mpfi_srcptr a,const unsigned long b)
{
  mpfi_t tmp;
  mpfi_set_ui(tmp,b);
  return(mpfi_cmp(a,tmp));
}

int mpfi_comp_si (mpfi_srcptr a,const unsigned long b)
{
  mpfi_t tmp;
  mpfi_set_si(tmp,b);
  return(mpfi_cmp(a,tmp));
}

int mpfi_is_pos_default(mpfi_srcptr a)
{
  return((MPFR_SIGN_PART(&(a->left))>0));
}

int mpfi_is_neg_default(mpfi_srcptr a)
{
  return((MPFR_SIGN_PART(&(a->right))<0));
}

int mpfi_is_zero_default(mpfi_srcptr a)
{
  return((MPFR_SIGN_PART(&(a->right))==0) && 
	 (MPFR_SIGN_PART(&(a->left))==0));
}

int (*mpfi_is_pos)  (mpfi_srcptr)=mpfi_is_pos_default;
int (*mpfi_is_neg)  (mpfi_srcptr)=mpfi_is_neg_default;
int (*mpfi_is_zero)  (mpfi_srcptr)=mpfi_is_zero_default;


/*  interval manipulation */

static int mpfi_error=0;

#define MPFI_ERROR(s) \
{\
if(!mpfi_error) mpfi_error=1;\
}

#define MPFI_CHECK_ERROR_INT(x,s) \
{\
if ((MPFI_HASZERO(x))) MPFI_ERROR(s);\
}

#ifdef MPFI_USE_CHECK_ERROR
#define MPFI_CHECK_ERROR(x,s) MPFI_CHECK_ERROR_INT(x,s)
#else 
#define MPFI_CHECK_ERROR(x,s) {}
#endif

#ifdef MPFI_USE_CHECK_ERROR
void mpfi_reset_error() 
{
  mpfi_error=0;
}

int mpfi_is_error()
{
  return(mpfi_error==1);
}
#endif

void mpfi_init (mpfi_t x)
{
  mpfr_init(&(x->left));
  mpfr_init(&(x->right));
}

void mpfi_init2 (mpfi_t x, mp_prec_t p)
{
  mpfr_init2 (&(x->left), p);
  mpfr_init2 (&(x->right), p);
}

void mpfi_set(mpfi_ptr a, mpfi_srcptr b)
{
  mpfr_set(&(a->left),&(b->left),MPFI_RNDD);
  mpfr_set(&(a->right),&(b->right),MPFI_RNDU);
  MPFI_CHECK_ERROR(a,"mpfi_set");
}

void mpfi_set_si(mpfi_ptr a,const long b)
{
  mpfr_set_si(&(a->left),b,MPFI_RNDD);
  mpfr_set_si(&(a->right),b,MPFI_RNDU);
  MPFI_CHECK_ERROR(a,"mpfi_set_si");
}

void mpfi_set_ui(mpfi_ptr a,const unsigned long b)
{
  mpfr_set_ui(&(a->left),b,MPFI_RNDD);
  mpfr_set_ui(&(a->right),b,MPFI_RNDU);
  MPFI_CHECK_ERROR(a,"mpfi_set_ui");
}

void mpfi_set_d(mpfi_ptr a, const double b)
{
  mpfr_set_d(&(a->left),b,MPFI_RNDD);
  mpfr_set_d(&(a->right),b,MPFI_RNDU);
  MPFI_CHECK_ERROR(a,"mpfi_set_ui");
}

void mpfi_add (mpfi_ptr a, mpfi_srcptr b, mpfi_srcptr c)
{
  if (MPFI_ISNULL(c)) {
    mpfi_set (a, b);
  }
  else {
    if (MPFI_ISNULL(b)) {
      mpfi_set (a, c);
    }
    else {
      mpfr_add (&(a->left), &(b->left), &(c->left), MPFI_RNDD);
      mpfr_add (&(a->right), &(b->right), &(c->right), MPFI_RNDU);
    }
  }
  MPFI_CHECK_ERROR(a,"mpfi_add");
}

void mpfi_sub (mpfi_ptr a, mpfi_srcptr b, mpfi_srcptr c)
{
  /* if using temporary variables: check that their precision agrees
     with that of input! */
  mpfr_sub (&(a->left), &(b->left), &(c->left), MPFI_RNDD);
  mpfr_sub (&(a->right), &(b->right), &(c->right), MPFI_RNDU);
  MPFI_CHECK_ERROR(a,"mpfi_sub");
}

void mpfi_mul (mpfi_ptr a, mpfi_srcptr b, mpfi_srcptr c)
{
  mpfr_t u, v;
  int in_place;

  if (MPFI_ISNULL(b) || MPFI_ISNULL(c)) {
    mpfi_set_ui (a, 0);
  }
  else {
    if (MPFI_ISPOS(c) || MPFI_ISPOS(b)) {
      /* works here even if a = b or a = c */
      mpfr_mul(&(a->left), &(b->left), &(c->left), MPFI_RNDD);
      mpfr_mul(&(a->right), &(b->right), &(c->right), MPFI_RNDU);
    }
    else {
      if (MPFI_MPFR_ISNEG(c)) {
	in_place = (b->right)._mp_d == (a->right)._mp_d;
	if (!in_place) u[0] = b->right;
	else {
	  mpfr_init2 (u, MPFR_PREC(&(b->right)));
	  mpfr_set (u, &(b->right), GMP_RNDD);
	}
	mpfr_mul (&(a->right), &(b->left), &(c->right), MPFI_RNDU);
	mpfr_mul (&(a->left), u, &(c->left), MPFI_RNDD);    
	if (in_place) mpfr_clear (u);
      }
      else {
	fprintf (stderr, "mpfi_mul: not yet implemented\n");
	exit (1);
      }  
    }
  }
  MPFI_CHECK_ERROR(a,"mpfi_mul");
}

void mpfi_div (mpfi_ptr a, mpfi_srcptr u, mpfi_srcptr c)
{
  mpfi_t b;
  int in_place;

  in_place = (u->left)._mp_d == (a->left)._mp_d;

  if (!in_place) {
    b[0]=(*u);
  }
  else {
    mpfi_init(b);
    mpfi_set(b,u);
  }
  if (!(MPFI_ISNULL(b))) {
    if (MPFI_ISPOS(c)) {
      mpfr_div(&(a->left),&(b->left),&(c->right),MPFI_RNDD);
      mpfr_div(&(a->right),&(b->right),&(c->left),MPFI_RNDU);
    }
    else {
      if (MPFI_MPFR_ISNEG(c)) {
	mpfr_div(&(a->right),&(b->left),&(c->left),MPFI_RNDU);
	mpfr_div(&(a->left),&(b->right),&(c->right),MPFI_RNDD);    
      }
      else {
	/* zero dans l'intervalle .... */
	MPFI_ERROR("mpfi_div : division by zero");
      }  
    }
  }
  MPFI_CHECK_ERROR(a,"mpfi_div");
  if (in_place) mpfi_clear (b);
}

void   mpfi_add_d(mpfi_ptr a, mpfi_srcptr b, const double c)
{
  mpfi_t tmp;
  mpfi_set_d(tmp,c);
  mpfi_add(a,b,tmp);
}

void   mpfi_sub_d(mpfi_ptr a, mpfi_srcptr b, const double c)
{
  mpfi_t tmp;
  mpfi_set_d(tmp,c);
  mpfi_sub(a,b,tmp);
}

void   mpfi_mul_d(mpfi_ptr a, mpfi_srcptr b, const double c)
{
  mpfi_t tmp;
  mpfi_set_d(tmp,c);
  mpfi_mul(a,b,tmp);
}

void mpfi_mul_ui (mpfi_ptr a, mpfi_srcptr b, unsigned int c)
{
  mpfr_mul_ui (&(a->left), &(b->left), c, GMP_RNDD);
  mpfr_mul_ui (&(a->right), &(b->right), c, GMP_RNDU);
}

void mpfi_sub_ui (mpfi_ptr a, mpfi_srcptr b, unsigned int c)
{
  mpfr_sub_ui (&(a->left), &(b->left), c, GMP_RNDD);
  mpfr_sub_ui (&(a->right), &(b->right), c, GMP_RNDU);
}

void mpfi_ui_div (mpfi_ptr a, unsigned int b, mpfi_srcptr c)
{
  if (MPFI_ISPOS(c) || MPFI_MPFR_ISNEG(c)) {
    mpfr_t tmp;
    int in_place = (a->left)._mp_d == (c->left)._mp_d;
    if (in_place) {
      mpfr_init2 (tmp, MPFR_PREC(&(a->left)));
      mpfr_set (tmp, &(a->left), GMP_RNDN);
    }
    else tmp[0] = a->left;
    mpfr_ui_div (&(a->left), b, &(c->right), GMP_RNDD);
    mpfr_ui_div (&(a->right), b, tmp, GMP_RNDU);
    if (in_place) mpfr_clear (tmp);
  }
  else
    MPFI_ERROR("mpfi_ui_div : division by zero");
}

void mpfi_div_d (mpfi_ptr a, mpfi_srcptr b, const double c)
{
  mpfi_t tmp;
  mpfi_set_d(tmp,c);
  mpfi_div(a,b,tmp);
}

void     mpfi_dadd(mpfi_ptr a,mpfi_srcptr b)
{
  mpfi_add(a,a,b);
}

void     mpfi_dsub(mpfi_ptr a,mpfi_srcptr b)
{
  mpfi_sub(a,a,b);
}

void     mpfi_dmul(mpfi_ptr a,mpfi_srcptr b)
{
  mpfi_mul(a,a,b);
}

void     mpfi_ddiv(mpfi_ptr a,mpfi_srcptr b)
{
  mpfi_div(a,a,b);
}

void mpfi_mul_2exp(mpfi_ptr a, mpfi_srcptr b,unsigned long c)
{
  mpfr_mul_2exp(&(a->left),&(b->left),c,MPFI_RNDD);
  mpfr_mul_2exp(&(a->right),&(b->right),c,MPFI_RNDU);
  MPFI_CHECK_ERROR(a,"mpfi_mul_2exp");
}

void mpfi_div_2exp(mpfi_ptr a, mpfi_srcptr b,unsigned long c)
{
  mpfr_div_2exp(&(a->left),&(b->left),c,MPFI_RNDD);
  mpfr_div_2exp(&(a->right),&(b->right),c,MPFI_RNDU);
  MPFI_CHECK_ERROR(a,"mpfi_mul_2exp");
}

void mpfi_neg(mpfi_ptr a, mpfi_srcptr b)
{
  mpfr_t tmp;

  mpfr_init2 (tmp, MPFR_PREC(&(b->left)));
  mpfr_set (tmp, &(b->left), MPFI_RNDD);
  mpfr_neg (&(a->left), &(b->right), MPFI_RNDD);
  mpfr_neg (&(a->right), tmp, MPFI_RNDU);
  MPFI_CHECK_ERROR (a,"mpfi_neg");
  mpfr_clear (tmp);
}

void mpfi_inv(mpfi_ptr a, mpfi_srcptr b)
{
  mpfi_ui_div(a, 1, b);
}

void   mpfi_set_left (mpfi_ptr a, mpfr_srcptr b)
{
  mpfr_set(&(a->left),b,MPFI_RNDD);
}

void   mpfi_set_right (mpfi_ptr a, mpfr_srcptr b)
{
  mpfr_set(&(a->right),b,MPFI_RNDU);
}

void   mpfi_get_left (mpfi_srcptr a, mpfr_ptr b)
{
  mpfr_set(b,&(a->left),MPFI_RNDD);  
}

void   mpfi_get_right (mpfi_srcptr a, mpfr_ptr b)
{
  mpfr_set(b,&(a->right),MPFI_RNDU);  
}

void mpfi_clear(mpfi_ptr a)
{
  mpfr_clear(&(a->right));  
  mpfr_clear(&(a->left));  
  
}

void mpfi_out_str (FILE *stream, int base, size_t n_digits, mpfi_srcptr op)
{
  fprintf (stream, "[");
  mpfr_out_str (stream, base, n_digits, &(op->left), GMP_RNDD);
  fprintf (stream, ",");
  mpfr_out_str (stream, base, n_digits, &(op->right), GMP_RNDU);
  fprintf (stream, "]");
}
