/* Test file for mpfr_mul.

Copyright 1999, 2001, 2002, 2003 Free Software Foundation, Inc.

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
#include <time.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#include "mpfr-test.h"

/* Workaround for sparc gcc 2.95.x bug, see notes in tadd.c. */
#define check(x,y,rnd_mode,px,py,pz,res)  _check(x,y,res,rnd_mode,px,py,pz)

/* checks that x*y gives the same results in double
   and with mpfr with 53 bits of precision */
static void
_check (double x, double y, double res, mp_rnd_t rnd_mode, unsigned int px,
        unsigned int py, unsigned int pz)
{
  double z2;
  mpfr_t xx, yy, zz;

  mpfr_init2 (xx, px);
  mpfr_init2 (yy, py);
  mpfr_init2 (zz, pz);
  mpfr_set_d(xx, x, rnd_mode);
  mpfr_set_d(yy, y, rnd_mode);
  mpfr_mul(zz, xx, yy, rnd_mode);
  z2 = mpfr_get_d1 (zz);
  if (res != z2)
    {
      printf ("mpfr_mul failed for x=%1.20e y=%1.20e with rnd_mode=%s\n",
              x, y, mpfr_print_rnd_mode (rnd_mode));
      printf ("correct is %1.20e, mpfr_mul gives %1.20e\n", res, z2);
      if (res != 0.0)
        exit (1);
    }
  mpfr_clear(xx); mpfr_clear(yy); mpfr_clear(zz);
}

static void
check53 (double x, double y, mp_rnd_t rnd_mode, double z1)
{
  double z2;
  mpfr_t xx, yy, zz;

  mpfr_init2 (xx, 53);
  mpfr_init2 (yy, 53);
  mpfr_init2 (zz, 53);
  mpfr_set_d (xx, x, rnd_mode);
  mpfr_set_d (yy, y, rnd_mode);
  mpfr_mul (zz, xx, yy, rnd_mode);
  z2 = mpfr_get_d1 (zz);
  if (z1 != z2 && !(Isnan(z1) && Isnan(z2)))
    {
      printf ("mpfr_mul failed for x=%1.20e y=%1.20e with rnd_mode=%s\n",
              x, y, mpfr_print_rnd_mode(rnd_mode));
      printf ("correct result is %1.20e, mpfr_mul gives %1.20e\n", z1, z2);
      exit (1);
    }
  mpfr_clear (xx);
  mpfr_clear (yy);
  mpfr_clear (zz);
}

/* checks that x*y gives the same results in double
   and with mpfr with 24 bits of precision */
static void
check24 (float x, float y, mp_rnd_t rnd_mode, float z1)
{
  float z2;
  mpfr_t xx, yy, zz;

  mpfr_init2 (xx, 24);
  mpfr_init2 (yy, 24);
  mpfr_init2 (zz, 24);
  mpfr_set_d (xx, x, rnd_mode);
  mpfr_set_d (yy, y, rnd_mode);
  mpfr_mul (zz, xx, yy, rnd_mode);
  z2 = (float) mpfr_get_d1 (zz);
  if (z1 != z2)
    {
      printf ("mpfr_mul failed for x=%1.0f y=%1.0f with prec=24 and"
	      "rnd=%s\n", x, y, mpfr_print_rnd_mode(rnd_mode));
      printf ("correct result is gives %.10e, mpfr_mul gives %.10e\n", z1, z2);
      exit (1);
    }
  mpfr_clear(xx);
  mpfr_clear(yy);
  mpfr_clear(zz);
}

/* the following examples come from the paper "Number-theoretic Test 
   Generation for Directed Rounding" from Michael Parks, Table 1 */
static void
check_float (void)
{
  float z;

  check24(8388609.0,  8388609.0, GMP_RNDN, 70368760954880.0);
  check24(16777213.0, 8388609.0, GMP_RNDN, 140737479966720.0);
  check24(8388611.0,  8388609.0, GMP_RNDN, 70368777732096.0);
  check24(12582911.0, 8388610.0, GMP_RNDN, 105553133043712.0);
  check24(12582914.0, 8388610.0, GMP_RNDN, 105553158209536.0);
  check24(13981013.0, 8388611.0, GMP_RNDN, 117281279442944.0);
  check24(11184811.0, 8388611.0, GMP_RNDN, 93825028587520.0);
  check24(11184810.0, 8388611.0, GMP_RNDN, 93825020198912.0);
  check24(13981014.0, 8388611.0, GMP_RNDN, 117281287831552.0);

  check24(8388609.0,  8388609.0, GMP_RNDZ, 70368760954880.0);
  check24(16777213.0, 8388609.0, GMP_RNDZ, 140737471578112.0);
  z = 70368777732096.0;
  check24(8388611.0,  8388609.0, GMP_RNDZ, z);
  check24(12582911.0, 8388610.0, GMP_RNDZ, 105553124655104.0);
  check24(12582914.0, 8388610.0, GMP_RNDZ, 105553158209536.0);
  check24(13981013.0, 8388611.0, GMP_RNDZ, 117281271054336.0);
  check24(11184811.0, 8388611.0, GMP_RNDZ, 93825028587520.0);
  check24(11184810.0, 8388611.0, GMP_RNDZ, 93825011810304.0);
  check24(13981014.0, 8388611.0, GMP_RNDZ, 117281287831552.0);

  check24(8388609.0,  8388609.0, GMP_RNDU, 70368769343488.0);
  check24(16777213.0, 8388609.0, GMP_RNDU, 140737479966720.0);
  check24(8388611.0,  8388609.0, GMP_RNDU, 70368786120704.0);
  check24(12582911.0, 8388610.0, GMP_RNDU, 105553133043712.0);
  check24(12582914.0, 8388610.0, GMP_RNDU, 105553166598144.0);
  check24(13981013.0, 8388611.0, GMP_RNDU, 117281279442944.0);
  check24(11184811.0, 8388611.0, GMP_RNDU, 93825036976128.0);
  check24(11184810.0, 8388611.0, GMP_RNDU, 93825020198912.0);
  check24(13981014.0, 8388611.0, GMP_RNDU, 117281296220160.0);

  check24(8388609.0,  8388609.0, GMP_RNDD, 70368760954880.0);
  check24(16777213.0, 8388609.0, GMP_RNDD, 140737471578112.0);
  check24(8388611.0,  8388609.0, GMP_RNDD, 70368777732096.0);
  check24(12582911.0, 8388610.0, GMP_RNDD, 105553124655104.0);
  check24(12582914.0, 8388610.0, GMP_RNDD, 105553158209536.0);
  check24(13981013.0, 8388611.0, GMP_RNDD, 117281271054336.0);
  check24(11184811.0, 8388611.0, GMP_RNDD, 93825028587520.0);
  check24(11184810.0, 8388611.0, GMP_RNDD, 93825011810304.0);
  check24(13981014.0, 8388611.0, GMP_RNDD, 117281287831552.0);
}

/* check sign of result */
static void
check_sign (void)
{
  mpfr_t a, b;

  mpfr_init2 (a, 53);
  mpfr_init2 (b, 53);
  mpfr_set_si (a, -1, GMP_RNDN);
  mpfr_set_ui (b, 2, GMP_RNDN);
  mpfr_mul(a, b, b, GMP_RNDN);
  if (mpfr_cmp_ui (a, 4) )
    {
      printf ("2.0*2.0 gives %1.20e\n", mpfr_get_d1 (a));
      exit (1);
    }
  mpfr_clear(a); mpfr_clear(b);
}

/* checks that the inexact return value is correct */
static void
check_exact (void)
{
  mpfr_t a, b, c, d;
  mp_prec_t prec;
  int i, inexact;
  mp_rnd_t rnd;

  mpfr_init (a);
  mpfr_init (b);
  mpfr_init (c);
  mpfr_init (d);

  mpfr_set_prec (a, 17);
  mpfr_set_prec (b, 17);
  mpfr_set_prec (c, 32);
  mpfr_set_str_binary (a, "1.1000111011000100e-1");
  mpfr_set_str_binary (b, "1.0010001111100111e-1");
  if (mpfr_mul (c, a, b, GMP_RNDZ))
    {
      printf ("wrong return value (1)\n");
      exit (1);
    }

  for (prec = 2; prec < 100; prec++)
    {
      mpfr_set_prec (a, prec);
      mpfr_set_prec (b, prec);
      mpfr_set_prec (c, 2 * prec - 2);
      mpfr_set_prec (d, 2 * prec);
      for (i = 0; i < 1000; i++)
        {
          mpfr_random (a);
          mpfr_random (b);
          rnd = randlimb () % 4;
          inexact = mpfr_mul (c, a, b, rnd);
          if (mpfr_mul (d, a, b, rnd)) /* should be always exact */
            {
              printf ("unexpected inexact return value\n");
              exit (1);
            }
          if ((inexact == 0) && mpfr_cmp (c, d))
            {
              printf ("inexact=0 but results differ\n");
              exit (1);
            }
          else if (inexact && (mpfr_cmp (c, d) == 0))
            {
              printf ("inexact!=0 but results agree\n");
              printf ("prec=%u rnd=%s a=", (unsigned int) prec,
                      mpfr_print_rnd_mode (rnd));
              mpfr_out_str (stdout, 2, 0, a, rnd);
              printf ("\nb=");
              mpfr_out_str (stdout, 2, 0, b, rnd);
              printf ("\nc=");
              mpfr_out_str (stdout, 2, 0, c, rnd);
              printf ("\nd=");
              mpfr_out_str (stdout, 2, 0, d, rnd);
              printf ("\n");
              exit (1);
            }
        }
    }

  mpfr_clear (a);
  mpfr_clear (b);
  mpfr_clear (c);
  mpfr_clear (d);
}

static void
check_max(void)
{
  mpfr_t xx, yy, zz;

  mpfr_init2(xx, 4);
  mpfr_init2(yy, 4);
  mpfr_init2(zz, 4);
  mpfr_set_str (xx, /*11.0/16*/ "0.68750", 10, GMP_RNDN);
  mpfr_mul_2si(xx, xx, MPFR_EMAX_DEFAULT/2, GMP_RNDN);
  mpfr_set_str (yy, /*11.0/16*/ "0.68750", 10, GMP_RNDN);
  mpfr_mul_2si(yy, yy, MPFR_EMAX_DEFAULT - MPFR_EMAX_DEFAULT/2 + 1, GMP_RNDN);
  mpfr_clear_flags();
  mpfr_mul(zz, xx, yy, GMP_RNDU);
  if (!(mpfr_overflow_p() && MPFR_IS_INF(zz)))
    {
      printf("check_max failed (should be an overflow)\n");
      exit(1);
    }

  mpfr_clear_flags();
  mpfr_mul(zz, xx, yy, GMP_RNDD);
  if (mpfr_overflow_p() || MPFR_IS_INF(zz))
    {
      printf("check_max failed (should NOT be an overflow)\n");
      exit(1);
    }
  mpfr_set_str (xx, /*15.0/16*/ "0.93750", 10, GMP_RNDN);
  mpfr_mul_2si(xx, xx, MPFR_EMAX_DEFAULT, GMP_RNDN);
  if (!(MPFR_IS_FP(xx) && MPFR_IS_FP(zz)))
    {
      printf("check_max failed (internal error)\n");
      exit(1);
    }
  if (mpfr_cmp(xx, zz) != 0)
    {
      printf("check_max failed: got ");
      mpfr_out_str(stdout, 2, 0, zz, GMP_RNDZ);
      printf(" instead of ");
      mpfr_out_str(stdout, 2, 0, xx, GMP_RNDZ);
      printf("\n");
      exit(1);
    }

  mpfr_clear(xx);
  mpfr_clear(yy);
  mpfr_clear(zz);
}

static void
check_min(void)
{
  mpfr_t xx, yy, zz;

  mpfr_init2(xx, 4);
  mpfr_init2(yy, 4);
  mpfr_init2(zz, 3);
  mpfr_set_str(xx, /*15.0/16*/ "0.9375", 10, GMP_RNDN);
  mpfr_mul_2si(xx, xx, MPFR_EMIN_DEFAULT/2, GMP_RNDN);
  mpfr_set_str(yy, /*15.0/16*/ "0.9375", 10, GMP_RNDN);
  mpfr_mul_2si(yy, yy, MPFR_EMIN_DEFAULT - MPFR_EMIN_DEFAULT/2 - 1, GMP_RNDN);
  mpfr_mul(zz, xx, yy, GMP_RNDD);
  if (mpfr_sgn(zz) != 0)
    {
      printf("check_min failed: got ");
      mpfr_out_str(stdout, 2, 0, zz, GMP_RNDZ);
      printf(" instead of 0\n");
      exit(1);
    }

  mpfr_mul(zz, xx, yy, GMP_RNDU);
  mpfr_set_str (xx, "0.5", 10, GMP_RNDN);
  mpfr_mul_2si(xx, xx, MPFR_EMIN_DEFAULT, GMP_RNDN);
  if (mpfr_sgn(xx) <= 0)
    {
      printf("check_min failed (internal error)\n");
      exit(1);
    }
  if (mpfr_cmp(xx, zz) != 0)
    {
      printf("check_min failed: got ");
      mpfr_out_str(stdout, 2, 0, zz, GMP_RNDZ);
      printf(" instead of ");
      mpfr_out_str(stdout, 2, 0, xx, GMP_RNDZ);
      printf("\n");
      exit(1);
    }

  mpfr_clear(xx);
  mpfr_clear(yy);
  mpfr_clear(zz);
}

static void
check_nans (void)
{
  mpfr_t  p, x, y;

  mpfr_init2 (x, 123L);
  mpfr_init2 (y, 123L);
  mpfr_init2 (p, 123L);

  /* nan * 0 == nan */
  mpfr_set_nan (x);
  mpfr_set_ui (y, 0L, GMP_RNDN);
  mpfr_mul (p, x, y, GMP_RNDN);
  ASSERT_ALWAYS (mpfr_nan_p (p));

  /* 1 * nan == nan */
  mpfr_set_ui (x, 1L, GMP_RNDN);
  mpfr_set_nan (y);
  mpfr_mul (p, x, y, GMP_RNDN);
  ASSERT_ALWAYS (mpfr_nan_p (p));

  /* 0 * +inf == nan */
  mpfr_set_ui (x, 0L, GMP_RNDN);
  mpfr_set_nan (y);
  mpfr_mul (p, x, y, GMP_RNDN);
  ASSERT_ALWAYS (mpfr_nan_p (p));

  /* +1 * +inf == +inf */
  mpfr_set_ui (x, 1L, GMP_RNDN);
  mpfr_set_inf (y, 1);
  mpfr_mul (p, x, y, GMP_RNDN);
  ASSERT_ALWAYS (mpfr_inf_p (p));
  ASSERT_ALWAYS (mpfr_sgn (p) > 0);

  /* -1 * +inf == -inf */
  mpfr_set_si (x, -1L, GMP_RNDN);
  mpfr_set_inf (y, 1);
  mpfr_mul (p, x, y, GMP_RNDN);
  ASSERT_ALWAYS (mpfr_inf_p (p));
  ASSERT_ALWAYS (mpfr_sgn (p) < 0);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (p);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  check_nans ();
  check_exact ();
  check_float ();

  check53(6.9314718055994530941514e-1, 0.0, GMP_RNDZ, 0.0);
  check53(0.0, 6.9314718055994530941514e-1, GMP_RNDZ, 0.0);
  check_sign();
  check53(-4.165000000e4, -0.00004801920768307322868063274915, GMP_RNDN, 2.0); 
  check53(2.71331408349172961467e-08, -6.72658901114033715233e-165, GMP_RNDZ,
	  -1.8251348697787782844e-172);
  check53(0.31869277231188065, 0.88642843322303122, GMP_RNDZ,
	  2.8249833483992453642e-1);
  check(8.47622108205396074254e-01, 3.24039313247872939883e-01, GMP_RNDU,
	28, 45, 2, 0.375);
  check(2.63978122803639081440e-01, 6.8378615379333496093e-1, GMP_RNDN,
	34, 23, 31, 0.180504585267044603);
  check(1.0, 0.11835170935876249132, GMP_RNDU, 6, 41, 36, 0.1183517093595583);
  check53(67108865.0, 134217729.0, GMP_RNDN, 9.007199456067584e15);
  check(1.37399642157394197284e-01, 2.28877275604219221350e-01, GMP_RNDN,
	49, 15, 32, 0.0314472340833162888);
  check(4.03160720978664954828e-01, 5.85483042917246621073e-01, GMP_RNDZ,
	51, 22, 32, 0.2360436821472831);
  check(3.90798504668055102229e-14, 9.85394674650308388664e-04, GMP_RNDN,
	46, 22, 12, 0.385027296503914762e-16);
  check(4.58687081072827851358e-01, 2.20543551472118792844e-01, GMP_RNDN,
	49, 3, 2, 0.09375);
  check_max();
  check_min();

  tests_end_mpfr ();
  return 0;
}
