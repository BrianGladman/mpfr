/* Test file for mpfr_add_[q,z], mpfr_sub_[q,z], mpfr_div_[q,z], mpfr_mul_[q,z]

Copyright 2004 Free Software Foundation.

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

#include <stdio.h>
#include <stdlib.h>
#include "mpfr-test.h"


static void
test_specialz (int (*mpfr_func)(mpfr_ptr, mpfr_srcptr, mpz_srcptr, mp_rnd_t),
	       void (*mpz_func)(mpz_ptr, mpz_srcptr, mpz_srcptr),
	       const char *op)
{
  mpfr_t x1, x2;
  mpz_t  z1, z2;
  int res;

  mpfr_inits2 (128, x1, x2, NULL);
  mpz_init (z1); mpz_init(z2);
  mpz_fac_ui (z1, 19); /* 19!+1 fits perfectly in a 128 bits mantissa */
  mpz_add_ui (z1, z1, 1);
  mpz_fac_ui (z2, 20); /* 20!+1 fits perfectly in a 128 bits mantissa */
  mpz_add_ui (z2, z2, 1);
  res = mpfr_set_z(x1, z1, GMP_RNDN);
  if (res)
    {
      printf("Specialz %s: set_z1 error\n", op);
      abort();
    }
  /* (19!+1) * (20!+1) fits in a 128 bits number */
  res = mpfr_func(x1, x1, z2, GMP_RNDN);
  if (res)
    {
      printf("Specialz %s: wrong inexact flag.\n", op);
      abort();
    }
  mpz_func(z1, z1, z2);
  res = mpfr_set_z (x2, z1, GMP_RNDN);
  if (res)
    {
      printf("Specialz %s: set_z2 error\n", op);
      abort();
    }
  if (mpfr_cmp(x1, x2))
    {
      printf("Specialz %s: results differ.\nx1=", op);
      mpfr_print_binary(x1);
      printf("\nx2=");
      mpfr_print_binary(x2);
      putchar('\n');
      abort();
    }
  
  mpz_clear (z1); mpz_clear(z2);
  mpfr_clears(x1, x2, NULL);
}

static void
test_genericz (mp_prec_t p0, mp_prec_t p1, unsigned int N,
	       int (*func)(mpfr_ptr, mpfr_srcptr, mpz_srcptr, mp_rnd_t),
	       const char *op)
{
  mp_prec_t prec;
  mpfr_t arg1, dst_big, dst_small, tmp;
  mpz_t  arg2;
  mp_rnd_t rnd;
  int inexact, compare, compare2;
  unsigned int n;

  mpfr_inits (arg1, dst_big, dst_small, tmp, NULL);
  mpz_init (arg2);

  for (prec = p0; prec <= p1; prec++)
    {
      mpfr_set_prec (arg1, prec);
      mpfr_set_prec (tmp, prec);
      mpfr_set_prec (dst_small, prec);

      for (n=0; n<N; n++)
        {
          mpfr_urandomb (arg1, RANDS);
	  mpz_urandomb (arg2, RANDS, 1024);
          rnd = randlimb () % 4;
          mpfr_set_prec (dst_big, 2*prec);
          compare = func(dst_big, arg1, arg2, rnd);
          if (mpfr_can_round (dst_big, 2*prec, rnd, rnd, prec))
            {
              mpfr_set (tmp, dst_big, rnd);
              inexact = func(dst_small, arg1, arg2, rnd);
              if (mpfr_cmp (tmp, dst_small))
                {
		  printf ("Results differ for prec=%u rnd_mode=%s and %s \n"
			  "arg1=",
                          (unsigned) prec, mpfr_print_rnd_mode (rnd), op);
		  mpfr_print_binary (arg1);
		  printf("\narg2=");  
		  mpz_out_str(stdout, 2, arg2);
		  printf ("\ngot      ");
		  mpfr_print_binary (dst_small);
		  printf ("\nexpected ");
		  mpfr_print_binary (tmp);
		  printf ("\napprox  ");
		  mpfr_print_binary (dst_big);
		  exit (1);
		}
	      compare2 = mpfr_cmp (tmp, dst_big);
	      /* if rounding to nearest, cannot know the sign of t - f(x)
		 because of composed rounding: y = o(f(x)) and t = o(y) */
	      if (compare * compare2 >= 0)
		compare = compare + compare2;
	      else
		compare = inexact; /* cannot determine sign(t-f(x)) */
	      if (((inexact == 0) && (compare != 0)) ||
		  ((inexact > 0) && (compare <= 0)) ||
		  ((inexact < 0) && (compare >= 0)))
		{
		  printf ("Wrong inexact flag for rnd=%s and %s:\n"
			  "expected %d, got %d\n", 
			  mpfr_print_rnd_mode (rnd), op, compare, inexact);
		  printf ("\narg1="); mpfr_print_binary (arg1);
		  printf ("\narg2="); mpz_out_str(stdout, 2, arg2);
		  printf ("\ndstl="); mpfr_print_binary (dst_big); 
		  printf ("\ndsts="); mpfr_print_binary (dst_small); 
                  printf ("\ntmp ="); mpfr_print_binary (tmp);
		  exit (1);
		}
	    }
	}
    }

  mpz_clear (arg2);
  mpfr_clears (arg1, dst_big, dst_small, tmp, NULL);
}

static void
test_genericq (mp_prec_t p0, mp_prec_t p1, unsigned int N,
               int (*func)(mpfr_ptr, mpfr_srcptr, mpq_srcptr, mp_rnd_t),
	       const char *op)
{
  mp_prec_t prec;
  mpfr_t arg1, dst_big, dst_small, tmp;
  mpq_t  arg2;
  mp_rnd_t rnd;
  int inexact, compare, compare2;
  unsigned int n;

  mpfr_inits (arg1, dst_big, dst_small, tmp, NULL);
  mpq_init (arg2);

  for (prec = p0; prec <= p1; prec++)
    {
      mpfr_set_prec (arg1, prec);
      mpfr_set_prec (tmp, prec);
      mpfr_set_prec (dst_small, prec);

      for (n=0; n<N; n++)
        {
          mpfr_urandomb (arg1, RANDS);
          mpq_set_ui (arg2, randlimb (), randlimb() );
          rnd = randlimb () % 4;
          mpfr_set_prec (dst_big, prec+10);
          compare = func(dst_big, arg1, arg2, rnd);
          if (mpfr_can_round (dst_big, prec+10, rnd, rnd, prec))
            {
              mpfr_set (tmp, dst_big, rnd);
              inexact = func(dst_small, arg1, arg2, rnd);
              if (mpfr_cmp (tmp, dst_small))
                {
                  printf ("Results differ for prec=%u rnd_mode=%s and %s\n"
			  "arg1=",
                          (unsigned) prec, mpfr_print_rnd_mode (rnd), op);
                  mpfr_print_binary (arg1);
                  printf("\narg2=");
                  mpq_out_str(stdout, 2, arg2);
                  printf ("\ngot      ");
                  mpfr_print_binary (dst_small);
                  printf ("\nexpected ");
                  mpfr_print_binary (tmp);
                  printf ("\napprox  ");
                  mpfr_print_binary (dst_big);
                  exit (1);
                }
              compare2 = mpfr_cmp (tmp, dst_big);
              /* if rounding to nearest, cannot know the sign of t - f(x)
                 because of composed rounding: y = o(f(x)) and t = o(y) */
              if (compare * compare2 >= 0)
                compare = compare + compare2;
              else
                compare = inexact; /* cannot determine sign(t-f(x)) */
              if (((inexact == 0) && (compare != 0)) ||
                  ((inexact > 0) && (compare <= 0)) ||
                  ((inexact < 0) && (compare >= 0)))
                {
                  printf ("Wrong inexact flag for rnd=%s and %s:\n"
			  "expected %d, got %d", 
			  mpfr_print_rnd_mode (rnd), op, compare, inexact);
                  printf ("\narg1="); mpfr_print_binary (arg1);
                  printf ("\narg2="); mpq_out_str(stdout, 2, arg2);
                  printf ("\ndstl="); mpfr_print_binary (dst_big);
                  printf ("\ndsts="); mpfr_print_binary (dst_small);
                  printf ("\ntmp ="); mpfr_print_binary (tmp);
		  putchar('\n');
                  exit (1);
                }
            }
        }
    }

  mpq_clear (arg2);
  mpfr_clears (arg1, dst_big, dst_small, tmp, NULL);
}

static void
test_specialq (mp_prec_t p0, mp_prec_t p1, unsigned int N,
	       int (*mpfr_func)(mpfr_ptr, mpfr_srcptr, mpq_srcptr, mp_rnd_t),
               void (*mpq_func)(mpq_ptr, mpq_srcptr, mpq_srcptr),
	       const char *op)
{
  mpfr_t fra, frb, frq;
  mpq_t  q1, q2, qr;
  unsigned int n;
  mp_prec_t prec;

  for(prec = p0 ; prec < p1 ; prec++)
    {
      mpfr_inits2 (prec, fra, frb, frq, NULL);
      mpq_init (q1); mpq_init(q2); mpq_init (qr);
      
      for( n = 0 ; n < N ; n++)
	{
	  mpq_set_ui(q1, randlimb(), randlimb() );
	  mpq_set_ui(q2, randlimb(), randlimb() );
	  mpq_func (qr, q1, q2);
	  mpfr_set_q (fra, q1, GMP_RNDD);
	  mpfr_func (fra, fra, q2, GMP_RNDD);
	  mpfr_set_q (frb, q1, GMP_RNDU);
	  mpfr_func (frb, frb, q2, GMP_RNDU);
	  mpfr_set_q (frq, qr, GMP_RNDN);
	  /* We should have fra <= qr <= frb */
	  if ( (mpfr_cmp(fra, frq) > 0) || (mpfr_cmp (frq, frb) > 0))
	    {
	      printf("Range error for prec=%lu and %s", prec, op);
	      printf ("\nq1="); mpq_out_str(stdout, 2, q1);
	      printf ("\nq2="); mpq_out_str(stdout, 2, q2);
	      printf ("\nfr_dn="); mpfr_print_binary (fra);
	      printf ("\nfr_q ="); mpfr_print_binary (frq);
	      printf ("\nfr_up="); mpfr_print_binary (frb);
	      putchar('\n');
	      exit (1);
	    }
	}
      
      mpq_clear (q1); mpq_clear (q2); mpq_clear (qr);
      mpfr_clears (fra, frb, frq, NULL);
    }
}


int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  test_genericz (2, 150, 100, mpfr_add_z, "add");
  test_genericz (2, 150, 100, mpfr_sub_z, "sub");
  test_genericz (2, 150, 100, mpfr_mul_z, "mul");
  test_genericz (2, 150, 100, mpfr_div_z, "div");
  test_specialz (mpfr_add_z, mpz_add, "add");
  test_specialz (mpfr_sub_z, mpz_sub, "sub");
  test_specialz (mpfr_mul_z, mpz_mul, "mul");

  test_genericq (2, 150, 100, mpfr_add_q, "add");
  test_genericq (2, 150, 100, mpfr_sub_q, "sub");
  test_genericq (2, 150, 100, mpfr_mul_q, "mul");
  test_genericq (2, 150, 100, mpfr_div_q, "div");
  test_specialq (2, 150, 100, mpfr_mul_q, mpq_mul, "mul");
  test_specialq (2, 150, 100, mpfr_div_q, mpq_div, "div");
  test_specialq (2, 150, 100, mpfr_add_q, mpq_add, "add");
  test_specialq (2, 150, 100, mpfr_sub_q, mpq_sub, "sub");

  tests_end_mpfr ();
  return 0;
}

