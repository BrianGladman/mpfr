/* mpfr.h -- Include file for mpfr.

Copyright (C) 1999 PolKA project, Inria Lorraine and Loria

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

/* Cygnus does not know about *rand48 functions */
#ifdef __CYGWIN32__
#define mrand48 rand
#define drand48 rand
#define lrand48 rand
#define srand48 srand
#endif

/* Definition of rounding modes */

#define GMP_RNDN 0
#define GMP_RNDZ 1
#define GMP_RNDU 2
#define GMP_RNDD 3

/* Definitions of types and their semantics */

typedef unsigned long int mp_prec_t; /* easy to change if necessary */
typedef int mp_rnd_t;                /* preferred to char */

typedef struct {  
  mp_prec_t _mp_prec; /* WARNING : for the mpfr type, the precision */
                              /* should be understood as the number of BITS,*/
			      /* not the number of mp_limb_t's. This means  */
			      /* that the corresponding number of allocated
				 limbs is >= ceil(_mp_prec/BITS_PER_MP_LIMB) */
  mp_size_t _mp_size;         /* ABSSIZE(.) is the number of allocated 
				 limbs the field _mp_d points to.
				 The sign is that of _mp_size.
				 The number 0 is such that _mp_d[k-1]=0
				 where k = ceil(_mp_prec/BITS_PER_MP_LIMB) */
  mp_exp_t _mp_exp; 
  mp_limb_t *_mp_d;
}
__mpfr_struct; 

/*
   The number represented is

    sign(_mp_size)*(_mp_d[k-1]/B+_mp_d[k-2]/B^2+...+_mp_d[0]/B^k)*2^_mp_exp

   where k=ceil(_mp_prec/BITS_PER_MP_LIMB) and B=2^BITS_PER_MP_LIMB.

   For the msb (most significant bit) normalized representation, we must have
   _mp_d[k-1]>=B/2, unless the number is zero (in that case its sign is still
   given by sign(_mp_size)).

   We must also have the last k*BITS_PER_MP_LIMB-_mp_prec bits set to zero.
*/

typedef __mpfr_struct mpfr_t[1]; 
typedef __mpfr_struct *mpfr_ptr; 
typedef __gmp_const __mpfr_struct *mpfr_srcptr;


/* Prototypes */

#ifndef _PROTO
#if defined (__STDC__) || defined (__cplusplus)
#define _PROTO(x) x
#else
#define _PROTO(x) ()
#endif
#endif

/* bit 31 of _mp_size is used for sign,
   bit 30 of _mp_size is used for Nan flag,
   remaining bits are used to store the number of allocated limbs */
#define FLAG_NAN(x) (((x)->_mp_size >> 30)&1)
#define SET_NAN(x) ((x)->_mp_size |= (1<<30))
#define ABSSIZE(x) ((x)->_mp_size & ((1<<30)-1))
#define SIZE(x) ((x)->_mp_size)
#define EXP(x) ((x)->_mp_exp)
#define MANT(x) ((x)->_mp_d)
#define MPFR_SIGN(x) (((x)->_mp_size >> 31) ? -1 : 1)
#define ISNONNEG(x) (MPFR_SIGN(x)>=0)
#define ISNEG(x) (MPFR_SIGN(x)==-1)
#define CHANGE_SIGN(x) (SIZE(x) = SIZE(x) ^ (1<<31))
#define PREC(x) ((x)->_mp_prec)
#define NOTZERO(x) (MANT(x)[(PREC(x)-1)/BITS_PER_MP_LIMB])
#define SET_ZERO(x) (MANT(x)[(PREC(x)-1)/BITS_PER_MP_LIMB] = 0)
#define mpfr_sgn(x) ((NOTZERO(x)) ? MPFR_SIGN(x) : 0)

/* reallocates the mantissa of x to q bits and sets the precision to q */
#define _mpfr_realloc(x, q) { \
    (x)->_mp_d = (mp_ptr) (*_mp_reallocate_func) \
       ((x)->_mp_d, (PREC(x)-1)>>3, (q+7)>>3); \
    PREC(x) = q; }

void mpfr_init2 _PROTO ((mpfr_ptr, mp_prec_t));
void mpfr_init _PROTO ((mpfr_ptr));
int mpfr_round_raw _PROTO ((mp_limb_t *, mp_limb_t *, mp_prec_t, int,
			     mp_prec_t, mp_rnd_t));
int mpfr_round_raw2 _PROTO((mp_limb_t *, mp_prec_t, int, mp_rnd_t, mp_prec_t));
void mpfr_round _PROTO ((mpfr_ptr, mp_rnd_t, mp_prec_t)); 
int mpfr_can_round _PROTO ((mpfr_ptr, mp_prec_t, mp_rnd_t, mp_rnd_t,
			    mp_prec_t));
int mpfr_can_round_raw _PROTO ((mp_limb_t *, mp_prec_t, int, mp_prec_t, 
				mp_rnd_t, mp_rnd_t, mp_prec_t));
void mpfr_set_d _PROTO ((mpfr_ptr, double, mp_rnd_t)); 
int mpfr_set_z _PROTO ((mpfr_ptr, mpz_srcptr, mp_rnd_t)); 
mp_exp_t mpz_set_fr _PROTO ((mpz_ptr, mpfr_srcptr)); 
void mpfr_set_q _PROTO ((mpfr_ptr, mpq_srcptr, mp_rnd_t)); 
double mpfr_get_d _PROTO ((mpfr_srcptr)); 
double mpfr_get_d2 _PROTO ((mpfr_srcptr, long)); 
void mpfr_set_f _PROTO ((mpfr_ptr, mpf_srcptr, mp_rnd_t));
void mpfr_set_si _PROTO ((mpfr_ptr, long, mp_rnd_t)); 
void mpfr_set_ui _PROTO ((mpfr_ptr, unsigned long, mp_rnd_t)); 
void mpfr_print_raw _PROTO ((mpfr_srcptr)); 
void mpfr_random _PROTO ((mpfr_ptr)); 
void mpfr_random2 _PROTO ((mpfr_ptr, mp_size_t, mp_exp_t)); 
void mpfr_urandomb _PROTO ((mpfr_ptr, gmp_randstate_t)); 
void mpfr_clear _PROTO ((mpfr_ptr)); 
void mpfr_set_str_raw _PROTO ((mpfr_ptr, char *));
int mpfr_set_str _PROTO ((mpfr_ptr, char *, int, mp_rnd_t));
int mpfr_init_set_str _PROTO ((mpfr_ptr, char *, int, mp_rnd_t));
size_t mpfr_inp_str _PROTO ((mpfr_ptr, FILE *, int, mp_rnd_t));
void mpfr_get_str_raw _PROTO ((char *, mpfr_srcptr));
char* mpfr_get_str _PROTO ((char *, mp_exp_t *, int, size_t, mpfr_srcptr, mp_rnd_t));
size_t mpfr_out_str _PROTO ((FILE *, int, size_t, mpfr_srcptr, mp_rnd_t));
void mpfr_mul _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
int mpfr_pow_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
int mpfr_ui_pow_ui _PROTO ((mpfr_ptr, unsigned long int, unsigned long int,
			     mp_rnd_t));
mp_limb_t mpn_divrem_n _PROTO ((mp_limb_t *, mp_limb_t *, mp_limb_t *, mp_size_t));
mp_size_t kara_sqrtrem _PROTO ((mp_limb_t *, mp_limb_t *, mp_limb_t *, mp_size_t));
void mpfr_div _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
void mpfr_agm _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
int mpfr_sqrt _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_sqrt_ui _PROTO ((mpfr_ptr, unsigned long, mp_rnd_t));  
void mpfr_add _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
void mpfr_add_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long, mp_rnd_t));
void mpfr_sub_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long, mp_rnd_t));
int mpfr_add_one_ulp _PROTO ((mpfr_ptr));
void mpfr_sub _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
void mpfr_ui_sub _PROTO ((mpfr_ptr, unsigned long, mpfr_srcptr, mp_rnd_t));
void mpfr_reldiff _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
void mpfr_set4 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t, int));
void mpfr_const_pi _PROTO ((mpfr_ptr, mp_rnd_t));
void mpfr_const_log2 _PROTO ((mpfr_ptr, mp_rnd_t));
int mpfr_log _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_exp _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_zeta _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
void mpfr_mul_ui _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
void mpfr_set_machine_rnd_mode _PROTO ((mp_rnd_t));
int mpfr_cmp3 _PROTO ((mpfr_srcptr, mpfr_srcptr, long int));
int mpfr_cmp_ui_2exp _PROTO ((mpfr_srcptr, unsigned long int, int));
int mpfr_cmp_si_2exp _PROTO ((mpfr_srcptr, long int, int));
int mpfr_cmp2 _PROTO ((mpfr_srcptr, mpfr_srcptr));
void mpfr_mul_2exp _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int,mp_rnd_t));
void mpfr_div_2exp _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int,mp_rnd_t));
void mpfr_set_prec _PROTO((mpfr_ptr, mp_prec_t));
void mpfr_set_prec_raw _PROTO((mpfr_ptr, mp_prec_t));
void mpfr_set_default_prec _PROTO((mp_prec_t));
extern mp_prec_t __gmp_default_fp_bit_precision;
extern mp_rnd_t __gmp_default_rounding_mode;
char * mpfr_print_rnd_mode _PROTO((mp_rnd_t)); 
void mpfr_neg _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t)); 
int mpfr_sub_one_ulp _PROTO((mpfr_ptr)); 
int mpfr_div_ui _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t)); 
void mpfr_ui_div _PROTO((mpfr_ptr, unsigned long int, mpfr_srcptr, mp_rnd_t)); 
mp_prec_t mpfr_get_prec _PROTO((mpfr_t));
void mpfr_set_default_rounding_mode _PROTO((mp_rnd_t));
int mpfr_eq(mpfr_srcptr, mpfr_srcptr, unsigned long); 
void mpfr_floor _PROTO((mpfr_t, mpfr_srcptr));
void mpfr_trunc _PROTO((mpfr_t, mpfr_srcptr));
void mpfr_ceil _PROTO((mpfr_t, mpfr_srcptr));

#define mpfr_cmp_ui(b,i) mpfr_cmp_ui_2exp(b,i,0)
#define mpfr_cmp_si(b,i) mpfr_cmp_si_2exp(b,i,0)
#define mpfr_set(a,b,r) mpfr_set4(a,b,r,MPFR_SIGN(b))
#define mpfr_abs(a,b,r) mpfr_set4(a,b,r,1)
#define mpfr_cmp(b,c) mpfr_cmp3(b,c,1)

#if (BITS_PER_MP_LIMB==32)
#define MPFR_LIMBS_PER_DOUBLE 2
#elif (BITS_PER_MP_LIMB==64)
#define MPFR_LIMBS_PER_DOUBLE 1
#endif

/* gmp-2.0.2 had only one threshold for both multiplication and squaring */
#ifndef KARATSUBA_MUL_THRESHOLD 
#ifdef KARATSUBA_THRESHOLD
#define KARATSUBA_MUL_THRESHOLD KARATSUBA_THRESHOLD
#else
#define KARATSUBA_MUL_THRESHOLD 16
#endif
#endif

#define mpfr_init_set_si(x, i, rnd) \
 do { mpfr_init(x); mpfr_set_si((x), (i), (rnd)); } while (0)
#define mpfr_init_set_ui(x, i, rnd) \
 do { mpfr_init(x); mpfr_set_ui((x), (i), (rnd)); } while (0)
#define mpfr_init_set_d(x, d, rnd) \
 do { mpfr_init(x); mpfr_set_d((x), (d), (rnd)); } while (0)
#define mpfr_init_set(x, y, rnd) \
 do { mpfr_init(x); mpfr_set((x), (y), (rnd)); } while (0)
#define mpfr_init_set_f(x, y, rnd) \
 do { mpfr_init(x); mpfr_set_f((x), (y), (rnd)); } while (0)
#define mpfr_init_set_str_raw(x, y, rnd) \
 do { mpfr_init(x); mpfr_set_str_raw((x), (y), (rnd)); } while (0)

