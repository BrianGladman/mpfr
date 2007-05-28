/* Test file for mpfr_exp2.

Copyright 2001, 2002, 2003, 2004, 2006, 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_exp2
#include "tgeneric.c"

static void
special_overflow (void)
{
  mpfr_t x, y;
  int inex;
  mp_exp_t emin, emax;

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();

  set_emin (-125);
  set_emax (128);

  mpfr_init2 (x, 24);
  mpfr_init2 (y, 24);

  mpfr_set_str_binary (x, "0.101100100000000000110100E15");
  inex = mpfr_exp2 (y, x, GMP_RNDN);
  if (!mpfr_inf_p (y) || inex <= 0)
    {
      printf ("Overflow error.\n");
      mpfr_dump (y);
      printf ("inex = %d\n", inex);
      exit (1);
    }

  mpfr_clear (y);
  mpfr_clear (x);
  set_emin (emin);
  set_emax (emax);
}

static void
emax_m_eps (void)
{
  if (mpfr_get_emax () <= LONG_MAX)
    {
      mpfr_t x, y;
      int inex, ov;

      mpfr_init2 (x, sizeof(mp_exp_t) * CHAR_BIT * 4);
      mpfr_init2 (y, 8);
      mpfr_set_si (x, mpfr_get_emax (), GMP_RNDN);

      mpfr_clear_flags ();
      inex = mpfr_exp2 (y, x, GMP_RNDN);
      ov = mpfr_overflow_p ();
      if (!ov || !mpfr_inf_p (y) || inex <= 0)
        {
          printf ("Overflow error for x = emax, GMP_RNDN.\n");
          mpfr_dump (y);
          printf ("inex = %d, %soverflow\n", inex, ov ? "" : "no ");
          exit (1);
        }

      mpfr_nextbelow (x);

      mpfr_clear_flags ();
      inex = mpfr_exp2 (y, x, GMP_RNDN);
      ov = mpfr_overflow_p ();
      if (!ov || !mpfr_inf_p (y) || inex <= 0)
        {
          printf ("Overflow error for x = emax - eps, GMP_RNDN.\n");
          mpfr_dump (y);
          printf ("inex = %d, %soverflow\n", inex, ov ? "" : "no ");
          exit (1);
        }

      mpfr_clear_flags ();
      inex = mpfr_exp2 (y, x, GMP_RNDD);
      ov = mpfr_overflow_p ();
      if (ov || mpfr_inf_p (y) || inex >= 0 ||
          (mpfr_nextabove (y), !mpfr_inf_p (y)))
        {
          printf ("Overflow error for x = emax - eps, GMP_RNDD.\n");
          mpfr_dump (y);
          printf ("inex = %d, %soverflow\n", inex, ov ? "" : "no ");
          exit (1);
        }

      mpfr_clear (x);
      mpfr_clear (y);
    }
}

static void
exp_range (void)
{
  mpfr_t x;

  set_emin (3);
  mpfr_init2 (x, 8);
  mpfr_set_ui (x, 5, GMP_RNDN);
  mpfr_exp2 (x, x, GMP_RNDN);
  set_emin (MPFR_EMIN_MIN);
  if (mpfr_nan_p (x) || mpfr_cmp_ui (x, 32) != 0)
    {
      printf ("Error in mpfr_exp2 for x = 5, with emin = 3\n");
      printf ("Expected 32, got ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }
  mpfr_clear (x);
}

int
main (int argc, char *argv[])
{
  mpfr_t x, y;
  mp_exp_t emin, emax;

  tests_start_mpfr ();

  special_overflow ();
  emax_m_eps ();
  exp_range ();

  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_ui (x, 4, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDD);
  mpfr_exp2 (y, x, GMP_RNDU);
  if (mpfr_cmp_ui (y, 16) != 0)
    {
      printf ("Error for 2^4\n");
      exit (1);
    }

  mpfr_set_si (x, -4, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDD);
  mpfr_exp2 (y, x, GMP_RNDU);
  if (mpfr_cmp_ui_2exp (y, 1, -4) != 0)
    {
      printf ("Error for 2^(-4)\n");
      exit (1);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);
  mpfr_set_str (x, /*-1683977482443233.0 / 2199023255552.0*/
                "-7.6578429909351734750089235603809357e2", 10,GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDN);
  if (mpfr_cmp_str1 (y, "2.991959870867646566478e-231"))
    {
      printf ("Error for x=-1683977482443233/2^41\n");
      exit (1);
    }

  mpfr_set_prec (x, 10);
  mpfr_set_prec (y, 10);
  /* save emin */
  emin = mpfr_get_emin ();
  set_emin (-10);
  mpfr_set_si (x, -12, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 0) || mpfr_sgn (y) < 0)
    {
      printf ("Error for x=emin-2, RNDN\n");
      printf ("Expected +0\n");
      printf ("Got      "); mpfr_print_binary (y); puts ("");
      exit (1);
    }
  /* restore emin */
  set_emin (emin);

  /* save emax */
  emax = mpfr_get_emax ();
  set_emax (10);
  mpfr_set_ui (x, 11, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) < 0)
    {
      printf ("Error for x=emax+1, RNDN\n");
      exit (1);
    }
  /* restore emax */
  set_emax (emax);

  MPFR_SET_INF(x);
  MPFR_SET_POS(x);
  mpfr_exp2 (y, x, GMP_RNDN);
  if(!MPFR_IS_INF(y))
    {
      printf ("evaluation of function in INF does not return INF\n");
      exit (1);
    }

  MPFR_CHANGE_SIGN(x);
  mpfr_exp2 (y, x, GMP_RNDN);
  if(!MPFR_IS_ZERO(y))
    {
      printf ("evaluation of function in -INF does not return 0\n");
      exit (1);
    }

  MPFR_SET_NAN(x);
  mpfr_exp2 (y, x, GMP_RNDN);
  if(!MPFR_IS_NAN(y))
    {
      printf ("evaluation of function in NaN does not return NaN\n");
      exit (1);
    }

  if ((mp_exp_unsigned_t) 8 << 31 != 0 ||
      mpfr_get_emax () <= (mp_exp_unsigned_t) 100000 * 100000)
    {
      /* emax <= 10000000000 */
      mpfr_set_prec (x, 40);
      mpfr_set_prec (y, 40);
      mpfr_set_str (x, "10000000000.5", 10, GMP_RNDN);
      mpfr_clear_flags ();
      mpfr_exp2 (y, x, GMP_RNDN);
      if (!(MPFR_IS_INF (y) && MPFR_IS_POS (y) && mpfr_overflow_p ()))
        {
          printf ("exp2(10000000000.5) should overflow.\n");
          exit (1);
        }
    }

  mpfr_set_prec (x, 2);
  mpfr_set_prec (y, 2);
  mpfr_set_str_binary (x, "-1.0E-26");
  mpfr_exp2 (y, x, GMP_RNDD);
  mpfr_set_str_binary (x, "1.1E-1");
  if (mpfr_cmp (x, y))
    {
      printf ("Error for exp(-2^(-26)) for prec=2\n");
      exit (1);
    }

  test_generic (2, 100, 100);

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();
  return 0;
}
