/* mpfr.h -- Include file for mpfr.

Copyright (C) 1999 Free Software Foundation.

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

#ifndef __MPFR_H
#define __MPFR_H
#include <stdio.h>

/* Definition of rounding modes */

#define GMP_RNDN 0
#define GMP_RNDZ 1
#define GMP_RNDU 2
#define GMP_RNDD 3

/* Definition of exponent limits */

#define MPFR_EMAX_DEFAULT ((mp_exp_t) ((1UL<<31)-1))
#define MPFR_EMIN_DEFAULT (-(MPFR_EMAX_DEFAULT))

#define MPFR_EMIN_MIN MPFR_EMIN_DEFAULT
#define MPFR_EMIN_MAX MPFR_EMAX_DEFAULT
#define MPFR_EMAX_MIN MPFR_EMIN_DEFAULT
#define MPFR_EMAX_MAX MPFR_EMAX_DEFAULT

/* Flags */

#define MPFR_FLAGS_UNDERFLOW 1
#define MPFR_FLAGS_OVERFLOW 2
#define MPFR_FLAGS_NAN 4
#define MPFR_FLAGS_INEXACT 8
#define MPFR_FLAGS_ALL 15

/* Definitions of types and their semantics */

typedef unsigned long int mp_prec_t; /* easy to change if necessary */
typedef int mp_rnd_t;                /* preferred to char */

typedef struct {  
  mp_prec_t _mpfr_prec; /* WARNING : for the mpfr type, the precision */
                              /* should be understood as the number of BITS,*/
			      /* not the number of mp_limb_t's. This means  */
			      /* that the corresponding number of allocated
				 limbs is >= ceil(_mp_prec/BITS_PER_MP_LIMB) */
  mp_size_t _mpfr_size;         /* MPFR_ABSSIZE(.) is the number of allocated 
				 limbs the field _mp_d points to.
				 The sign is that of _mpfr_size.
				 The number 0 is such that _mp_d[k-1]=0
				 where k = ceil(_mp_prec/BITS_PER_MP_LIMB) */
  mp_exp_t _mpfr_exp; 
  mp_limb_t *_mpfr_d;
}
__mpfr_struct; 

/*
   The number represented is

    sign(_mpfr_size)*(_mpfr_d[k-1]/B+_mpfr_d[k-2]/B^2+...+_mpfr_d[0]/B^k)*2^_mpfr_exp

   where k=ceil(_mp_prec/BITS_PER_MP_LIMB) and B=2^BITS_PER_MP_LIMB.

   For the msb (most significant bit) normalized representation, we must have
   _mpfr_d[k-1]>=B/2, unless the number is zero (in that case its sign is still
   given by sign(_mpfr_size)).

   We must also have the last k*BITS_PER_MP_LIMB-_mp_prec bits set to zero.
*/

typedef __mpfr_struct mpfr_t[1]; 
typedef __mpfr_struct *mpfr_ptr; 
typedef __gmp_const __mpfr_struct *mpfr_srcptr;

#define MPFR_SIGN(x) (((x)->_mpfr_size >> 31) ? -1 : 1)

/* Prototypes */

#ifndef _PROTO
#if defined (__STDC__) || defined (__cplusplus)
#define _PROTO(x) x
#else
#define _PROTO(x) ()
#endif
#endif

#if defined (__cplusplus)
extern "C" {
#endif  

extern unsigned int __mpfr_flags;
extern mp_exp_t __mpfr_emin;
extern mp_exp_t __mpfr_emax;
mp_exp_t mpfr_get_emin _PROTO ((void));
int mpfr_set_emin _PROTO ((mp_exp_t));
mp_exp_t mpfr_get_emax _PROTO ((void));
int mpfr_set_emax _PROTO ((mp_exp_t));
void mpfr_clear_flags _PROTO ((void));
void mpfr_clear_underflow _PROTO ((void));
void mpfr_clear_overflow _PROTO ((void));
void mpfr_clear_nanflag _PROTO ((void));
void mpfr_clear_inexflag _PROTO ((void));
int mpfr_check_range _PROTO ((mpfr_ptr, mp_rnd_t));
int mpfr_underflow_p _PROTO ((void));
int mpfr_overflow_p _PROTO ((void));
int mpfr_nanflag_p _PROTO ((void));
int mpfr_inexflag_p _PROTO ((void));

void mpfr_init2 _PROTO ((mpfr_ptr, mp_prec_t));
void mpfr_init _PROTO ((mpfr_ptr));
int mpfr_round _PROTO ((mpfr_ptr, mp_rnd_t, mp_prec_t)); 
int mpfr_can_round _PROTO ((mpfr_ptr, mp_prec_t, mp_rnd_t, mp_rnd_t,
			    mp_prec_t));
int mpfr_set_d _PROTO ((mpfr_ptr, double, mp_rnd_t)); 
int mpfr_set_z _PROTO ((mpfr_ptr, mpz_srcptr, mp_rnd_t)); 
mp_exp_t mpz_set_fr _PROTO ((mpz_ptr, mpfr_srcptr)); 
void mpfr_set_q _PROTO ((mpfr_ptr, mpq_srcptr, mp_rnd_t)); 
double mpfr_get_d _PROTO ((mpfr_srcptr)); 
int mpfr_set_f _PROTO ((mpfr_ptr, mpf_srcptr, mp_rnd_t));
int mpfr_set_si _PROTO ((mpfr_ptr, long, mp_rnd_t));
int mpfr_set_ui _PROTO ((mpfr_ptr, unsigned long, mp_rnd_t));
void mpfr_print_raw _PROTO ((mpfr_srcptr)); 
void mpfr_random _PROTO ((mpfr_ptr));
void mpfr_random2 _PROTO ((mpfr_ptr, mp_size_t, mp_exp_t)); 
void mpfr_urandomb _PROTO ((mpfr_ptr, gmp_randstate_t)); 
void mpfr_clear _PROTO ((mpfr_ptr)); 
void mpfr_set_str_raw _PROTO ((mpfr_ptr, char *));
int mpfr_set_str _PROTO ((mpfr_ptr, __gmp_const char *, int, mp_rnd_t));
int mpfr_init_set_str _PROTO ((mpfr_ptr, char *, int, mp_rnd_t));
size_t mpfr_inp_str _PROTO ((mpfr_ptr, FILE *, int, mp_rnd_t));
char* mpfr_get_str _PROTO ((char *, mp_exp_t *, int, size_t, mpfr_srcptr, mp_rnd_t));
size_t mpfr_out_str _PROTO ((FILE *, int, size_t, mpfr_srcptr, mp_rnd_t));
int mpfr_mul _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
int mpfr_pow_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
int mpfr_ui_pow_ui _PROTO ((mpfr_ptr, unsigned long int, unsigned long int,
			     mp_rnd_t));
void mpfr_div _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
void mpfr_agm _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
int mpfr_sqrt _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_sqrt_ui _PROTO ((mpfr_ptr, unsigned long, mp_rnd_t));  
void mpfr_add _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
void mpfr_add_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long, mp_rnd_t));
int mpfr_sub_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long, mp_rnd_t));
void mpfr_add_one_ulp _PROTO ((mpfr_ptr));
int mpfr_sub _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
int mpfr_ui_sub _PROTO ((mpfr_ptr, unsigned long, mpfr_srcptr, mp_rnd_t));
void mpfr_reldiff _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
void mpfr_const_pi _PROTO ((mpfr_ptr, mp_rnd_t));
void mpfr_const_log2 _PROTO ((mpfr_ptr, mp_rnd_t));
int mpfr_const_euler _PROTO ((mpfr_ptr, mp_rnd_t));
int mpfr_log _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_exp _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_exp2 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_sin _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_cos _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_tan _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_mul_ui _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
void mpfr_set_machine_rnd_mode _PROTO ((mp_rnd_t));
int mpfr_cmp_ui_2exp _PROTO ((mpfr_srcptr, unsigned long int, int));
int mpfr_cmp_si_2exp _PROTO ((mpfr_srcptr, long int, int));
int mpfr_mul_2exp _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int,mp_rnd_t));
int mpfr_div_2exp _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int,mp_rnd_t));
int mpfr_set_prec _PROTO((mpfr_ptr, mp_prec_t));
void mpfr_set_prec_raw _PROTO((mpfr_ptr, mp_prec_t));
void mpfr_set_default_prec _PROTO((mp_prec_t));
mp_prec_t mpfr_get_default_prec _PROTO((void));
extern mp_prec_t __mpfr_default_fp_bit_precision;
extern mp_rnd_t __gmp_default_rounding_mode;
char * mpfr_print_rnd_mode _PROTO((mp_rnd_t)); 
int mpfr_neg _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 
void mpfr_sub_one_ulp _PROTO((mpfr_ptr)); 
int mpfr_div_ui _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t)); 
void mpfr_ui_div _PROTO((mpfr_ptr, unsigned long int, mpfr_srcptr, mp_rnd_t)); 
mp_prec_t mpfr_get_prec _PROTO((mpfr_srcptr));
void mpfr_set_default_rounding_mode _PROTO((mp_rnd_t));
int mpfr_eq _PROTO((mpfr_srcptr, mpfr_srcptr, unsigned long));
void mpfr_floor _PROTO((mpfr_ptr, mpfr_srcptr));
void mpfr_trunc _PROTO((mpfr_ptr, mpfr_srcptr));
void mpfr_ceil _PROTO((mpfr_ptr, mpfr_srcptr));
void mpfr_extract _PROTO((mpz_ptr, mpfr_srcptr, unsigned int));
void mpfr_swap _PROTO((mpfr_ptr, mpfr_ptr));
void mpfr_dump _PROTO((mpfr_srcptr, mp_rnd_t));
int mpfr_set4 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t, int));
int mpfr_cmp3 _PROTO ((mpfr_srcptr, mpfr_srcptr, int));
int mpfr_nan_p _PROTO((mpfr_srcptr));
int mpfr_inf_p _PROTO((mpfr_srcptr));
int mpfr_number_p _PROTO((mpfr_srcptr));
int mpfr_atanh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 
int mpfr_acosh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 
int mpfr_asinh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 
int mpfr_cosh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 
int mpfr_sinh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 
int mpfr_tanh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_asin _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_atan _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_factorial _PROTO ((mpfr_ptr, unsigned long int, mp_rnd_t));
int mpfr_fma _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
int mpfr_hypot _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
  /*En cours ...*/
int mpfr_pow _PROTO ((mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mp_rnd_t)); 
int mpfr_dim _PROTO ((mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mp_rnd_t)); 
int mpfr_max _PROTO ((mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mp_rnd_t)); 
int mpfr_min _PROTO ((mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mp_rnd_t)); 
int mpfr_log2 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 
int mpfr_log10 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 

#if defined (__cplusplus)
}
#endif  

#define mpfr_cmp_ui(b,i) mpfr_cmp_ui_2exp((b),(i),0)
#define mpfr_cmp_si(b,i) mpfr_cmp_si_2exp((b),(i),0)
#define mpfr_set(a,b,r) mpfr_set4(a,b,r,MPFR_SIGN(b))
#define mpfr_abs(a,b,r) mpfr_set4(a,b,r,1)
#define mpfr_cmp(b, c) mpfr_cmp3(b, c, 1)
#define mpfr_sgn(x) mpfr_cmp_ui(x,0)

#define mpfr_init_set_si(x, i, rnd) \
 ( mpfr_init(x), mpfr_set_si((x), (i), (rnd)) )
#define mpfr_init_set_ui(x, i, rnd) \
 ( mpfr_init(x), mpfr_set_ui((x), (i), (rnd)) )
#define mpfr_init_set_d(x, d, rnd) \
 ( mpfr_init(x), mpfr_set_d((x), (d), (rnd)) )
#define mpfr_init_set(x, y, rnd) \
 ( mpfr_init(x), mpfr_set((x), (y), (rnd)) )
#define mpfr_init_set_f(x, y, rnd) \
 ( mpfr_init(x), mpfr_set_f((x), (y), (rnd)) )

#endif
