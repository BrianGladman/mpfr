/* mpfr_cos -- cosine of a floating-point number

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
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

int 
#if __STDC__
mpfr_mul_mpz (mpfr_ptr y, mpfr_srcptr x, mpz_srcptr z,mp_rnd_t rnd_mode) 
#else
mpfr_mul_mpz (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mpz_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t;
  int res;
  mpfr_init(t);
  mpfr_set_z(t,z,rnd_mode);
  res=mpfr_mul(y,x,t,rnd_mode);
  mpfr_clear(t);
  return(res);
}

int 
#if __STDC__
mpfr_div_mpz (mpfr_ptr y, mpfr_srcptr x, mpz_srcptr z, mp_rnd_t rnd_mode) 
#else
mpfr_div_mpz (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mpz_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t;
  int res;
  mpfr_init(t);
  mpfr_set_z(t,z,rnd_mode);
  res=mpfr_div(y,x,t,rnd_mode);
  mpfr_clear(t);
  return(res);
}

int 
#if __STDC__
mpfr_add_mpz (mpfr_ptr y, mpfr_srcptr x, mpz_srcptr z, mp_rnd_t rnd_mode) 
#else
mpfr_add_mpz (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mpz_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t;
  int res;
  mpfr_init(t);
  mpfr_set_z(t,z,rnd_mode);
  res=mpfr_add(y,x,t,rnd_mode);
  mpfr_clear(t);
  return(res);
}

int 
#if __STDC__
mpfr_sub_mpz (mpfr_ptr y, mpfr_srcptr x, mpz_srcptr z,mp_rnd_t rnd_mode) 
#else
mpfr_sub_mpz (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mpz_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t;
  int res;
  mpfr_init(t);
  mpfr_set_z(t,z,rnd_mode);
  res=mpfr_sub(y,x,t,rnd_mode);
  mpfr_clear(t);
  return(res);
}

int 
#if __STDC__
mpfr_mul_mpq (mpfr_ptr y, mpfr_srcptr x, mpq_srcptr z,mp_rnd_t rnd_mode) 
#else
mpfr_mul_mpq (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mpq_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t;
  int res;
  mpfr_init(t);
  mpfr_set_q(t,z,rnd_mode);
  res=mpfr_mul(y,x,t,rnd_mode);
  mpfr_clear(t);
  return(res);
}

int 
#if __STDC__
mpfr_div_mpq (mpfr_ptr y, mpfr_srcptr x, mpq_srcptr z, mp_rnd_t rnd_mode) 
#else
mpfr_div_mpq (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mpq_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t;
  int res;
  mpfr_init(t);
  mpfr_set_q(t,z,rnd_mode);
  res=mpfr_div(y,x,t,rnd_mode);
  mpfr_clear(t);
  return(res);
}

int 
#if __STDC__
mpfr_add_mpq (mpfr_ptr y, mpfr_srcptr x, mpq_srcptr z, mp_rnd_t rnd_mode) 
#else
mpfr_add_mpq (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mpq_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t;
  int res;
  mpfr_init(t);
  mpfr_set_q(t,z,rnd_mode);
  res=mpfr_add(y,x,t,rnd_mode);
  mpfr_clear(t);
  return(res);
}

int 
#if __STDC__
mpfr_sub_mpq (mpfr_ptr y, mpfr_srcptr x, mpq_srcptr z,mp_rnd_t rnd_mode) 
#else
mpfr_sub_mpq (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mpq_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
  mpfr_t t;
  int res;
  mpfr_init(t);
  mpfr_set_q(t,z,rnd_mode);
  res=mpfr_sub(y,x,t,rnd_mode);
  mpfr_clear(t);
  return(res);
}
