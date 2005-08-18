/* Test file for mpfr_exp.

Copyright 1999, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

#ifdef CHECK_EXTERNAL
static int
test_exp (mpfr_ptr a, mpfr_srcptr b, mp_rnd_t rnd_mode)
{
  int res;
  int ok = rnd_mode == GMP_RNDN && mpfr_number_p (b) && mpfr_get_prec (a)>=53;
  if (ok)
    {
      mpfr_print_raw (b);
    }
  res = mpfr_exp (a, b, rnd_mode);
  if (ok)
    {
      printf (" ");
      mpfr_print_raw (a);
      printf ("\n");
    }
  return res;
}
#else
#define test_exp mpfr_exp
#endif

/* returns the number of ulp of error */
static void
check3 (const char *op, mp_rnd_t rnd, const char *res)
{
  mpfr_t x, y;

  mpfr_inits2 (53, x, y, NULL);
  /* y negative. If we forget to set the sign in mpfr_exp, we'll see it. */
  mpfr_set_si (y, -1, GMP_RNDN);
  mpfr_set_str1 (x, op);
  test_exp (y, x, rnd);
  if (mpfr_cmp_str1 (y, res) )
    {
      printf ("mpfr_exp failed for x=%s, rnd=%s\n",
              op, mpfr_print_rnd_mode (rnd));
      printf ("expected result is %s, got ", res);
      mpfr_out_str (stdout, 10, 0, y, GMP_RNDN);
      putchar('\n');
      exit (1);
    }
  mpfr_clears (x, y, NULL);
}

/* expx is the value of exp(X) rounded towards -infinity */
static void
check_worst_case (const char *Xs, const char *expxs)
{
  mpfr_t x, y;

  mpfr_inits2(53, x, y, NULL);
  mpfr_set_str1(x, Xs);
  test_exp(y, x, GMP_RNDD);
  if (mpfr_cmp_str1 (y, expxs))
    {
      printf ("exp(x) rounded towards -infinity is wrong\n");
      exit(1);
    }
  mpfr_set_str1(x, Xs);
  test_exp(x, x, GMP_RNDU);
  mpfr_nexttoinf (y);
  if (mpfr_cmp(x,y))
    {
      printf ("exp(x) rounded towards +infinity is wrong\n");
      exit(1);
    }
  mpfr_clears(x,y,NULL);
}

/* worst cases communicated by Jean-Michel Muller and Vincent Lefevre */
static int
check_worst_cases (void)
{
  mpfr_t x; mpfr_t y;

  mpfr_init(x);
  mpfr_set_prec (x, 53);

  check_worst_case("4.44089209850062517562e-16", "1.00000000000000022204");
  check_worst_case("6.39488462184069720009e-14", "1.00000000000006372680");
  check_worst_case("1.84741111297455401935e-12", "1.00000000000184718907");
  check_worst_case("1.76177628026265550074e-10", "1.00000000017617751702");
  check3("1.76177628026265550074e-10", GMP_RNDN, "1.00000000017617773906");
  check_worst_case("7.54175277499595900852e-10", "1.00000000075417516676");
  check3("7.54175277499595900852e-10", GMP_RNDN, "1.00000000075417538881");
  /* bug found by Vincent Lefe`vre on December 8, 1999 */
  check3("-5.42410311287441459172e+02", GMP_RNDN, "2.7176584868845723e-236");
  /* further cases communicated by Vincent Lefe`vre on January 27, 2000 */
  check3("-1.32920285897904911589e-10", GMP_RNDN, "0.999999999867079769622");
  check3("-1.44037948245738330735e-10", GMP_RNDN, "0.9999999998559621072757");
  check3("-1.66795910430705305937e-10", GMP_RNDZ, "0.9999999998332040895832");
  check3("-1.64310953745426656203e-10", GMP_RNDN, "0.9999999998356891017792");
  check3("-1.38323574826034659172e-10", GMP_RNDZ, "0.9999999998616764251835");
  check3("-1.23621668465115401498e-10", GMP_RNDZ, "0.9999999998763783315425");

  mpfr_set_prec (x, 601);
  mpfr_set_str (x, "0.88b6ba510e10450edc258748bc9dfdd466f21b47ed264cdf24aa8f64af1f3fad9ec2301d43c0743f534b5aa20091ff6d352df458ef1ba519811ef6f5b11853534fd8fa32764a0a6d2d0dd20@0", 16, GMP_RNDZ);
  mpfr_init2 (y, 601);
  mpfr_exp_2 (y, x, GMP_RNDD);
  mpfr_exp_3 (x, x, GMP_RNDD);
  if (mpfr_cmp (x, y))
    {
      printf ("mpfr_exp_2 and mpfr_exp_3 differ for prec=601\n");
      printf ("mpfr_exp_2 gives ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\nmpfr_exp_3 gives ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_prec (x, 13001);
  mpfr_set_prec (y, 13001);
  mpfr_random (x);
  mpfr_exp_3 (y, x, GMP_RNDN);
  mpfr_exp_2 (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("mpfr_exp_2 and mpfr_exp_3 differ for prec=13001\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
  return 0;
}

static void
compare_exp2_exp3 (int n)
{
  mpfr_t x, y, z; int prec; mp_rnd_t rnd;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);
  for (prec = 20; prec <= n; prec++)
    {
      mpfr_set_prec (x, prec);
      mpfr_set_prec (y, prec);
      mpfr_set_prec (z, prec);
      mpfr_random (x);
      rnd = (mp_rnd_t) RND_RAND();
      mpfr_exp_2 (y, x, rnd);
      mpfr_exp_3 (z, x, rnd);
      if (mpfr_cmp (y,z))
        {
          printf ("mpfr_exp_2 and mpfr_exp_3 disagree for rnd=%s and\nx=",
                  mpfr_print_rnd_mode (rnd));
          mpfr_print_binary (x);
          puts ("");
          printf ("mpfr_exp_2 gives  ");
          mpfr_print_binary (y);
          puts ("");
          printf ("mpfr_exp_3 gives ");
          mpfr_print_binary (z);
          puts ("");
          exit (1);
        }
  }

  /* bug found by Patrick Pe'lissier on 7 Jun 2004 */
  prec = 203780;
  mpfr_set_prec (x, prec);
  mpfr_set_prec (z, prec);
  mpfr_set_ui (x, 3, GMP_RNDN);
  mpfr_sqrt (x, x, GMP_RNDN);
  mpfr_sub_ui (x, x, 1, GMP_RNDN);
  mpfr_exp_3 (z, x, GMP_RNDN);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

#define TEST_FUNCTION test_exp
#include "tgeneric.c"

static void
check_special ()
{
  mpfr_t x, y, z;
  mp_exp_t emin, emax;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);

  /* check exp(NaN) = NaN */
  mpfr_set_nan (x);
  test_exp (y, x, GMP_RNDN);
  if (!mpfr_nan_p (y))
    {
      printf ("Error for exp(NaN)\n");
      exit (1);
    }

  /* check exp(+inf) = +inf */
  mpfr_set_inf (x, 1);
  test_exp (y, x, GMP_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) < 0)
    {
      printf ("Error for exp(+inf)\n");
      exit (1);
    }

  /* check exp(-inf) = +0 */
  mpfr_set_inf (x, -1);
  test_exp (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 0) || mpfr_sgn (y) < 0)
    {
      printf ("Error for exp(-inf)\n");
      exit (1);
    }

  /* Check overflow. Corner case of mpfr_exp_2 */
  mpfr_set_prec (x, 64);
  mpfr_set_emax (MPFR_EMAX_DEFAULT);
  mpfr_set_emin (MPFR_EMIN_DEFAULT);
  mpfr_set_str (x, 
    "0.1011000101110010000101111111010100001100000001110001100111001101E30",
                2, GMP_RNDN);
  mpfr_exp (x, x, GMP_RNDD);
  if (mpfr_cmp_str (x, 
".1111111111111111111111111111111111111111111111111111111111111111E1073741823",
                    2, GMP_RNDN) != 0)
    {
      printf ("Wrong overflow detection in mpfr_exp\n");
      mpfr_dump (x);
      exit (1);
    }
  /* Check underflow. Corner case of mpfr_exp_2 */
  mpfr_set_str (x,
"-0.1011000101110010000101111111011111010001110011110111100110101100E30",
                2, GMP_RNDN);
  mpfr_exp (x, x, GMP_RNDN);
  if (mpfr_cmp_str (x, "0.1E-1073741823", 2, GMP_RNDN) != 0)
    {
      printf ("Wrong underflow (1) detection in mpfr_exp\n");
      mpfr_dump (x);
      exit (1);
    }
  mpfr_set_str (x,
"-0.1011001101110010000101111111011111010001110011110111100110111101E30",
                2, GMP_RNDN);
  mpfr_exp (x, x, GMP_RNDN);
  if (mpfr_cmp_ui (x, 0) != 0)
    {
      printf ("Wrong underflow (2) detection in mpfr_exp\n");
      mpfr_dump (x);
      exit (1);
    }
  /* Check overflow. Corner case of mpfr_exp_3 */
  if (MPFR_PREC_MAX > MPFR_EXP_THRESHOLD+10) {
    mpfr_set_prec (x, MPFR_EXP_THRESHOLD+10);
    mpfr_set_str (x, 
     "0.1011000101110010000101111111010100001100000001110001100111001101E30",
                  2, GMP_RNDN);
    mpfr_clear_overflow ();
    mpfr_exp (x, x, GMP_RNDD);
    if (!mpfr_overflow_p ())
      {
        printf ("Wrong overflow detection in mpfr_exp_3\n");
        mpfr_dump (x);
        exit (1);
      }
    /* Check underflow. Corner case of mpfr_exp_3 */
    mpfr_set_str (x,
    "-0.1011000101110010000101111111011111010001110011110111100110101100E30",
                  2, GMP_RNDN);
    mpfr_clear_underflow ();
    mpfr_exp (x, x, GMP_RNDN);
    if (!mpfr_underflow_p ())
      {
        printf ("Wrong underflow detection in mpfr_exp_3\n");
        mpfr_dump (x);
        exit (1);
      }
    mpfr_set_prec (x, 53);
  }

  /* check overflow */
  emax = mpfr_get_emax ();
  set_emax (10);
  mpfr_set_ui (x, 7, GMP_RNDN);
  test_exp (y, x, GMP_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) < 0)
    {
      printf ("Error for exp(7) for emax=10\n");
      exit (1);
    }
  set_emax (emax);

  /* check underflow */
  emin = mpfr_get_emin ();
  set_emin (-10);
  mpfr_set_si (x, -9, GMP_RNDN);
  test_exp (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 0) || mpfr_sgn (y) < 0)
    {
      printf ("Error for exp(-9) for emin=-10\n");
      printf ("Expected +0\n");
      printf ("Got      "); mpfr_print_binary (y); puts ("");
      exit (1);
    }
  set_emin (emin);

  /* check case EXP(x) < -precy */
  mpfr_set_prec (y, 2);
  mpfr_set_str_binary (x, "-0.1E-3");
  test_exp (y, x, GMP_RNDD);
  if (mpfr_cmp_ui_2exp (y, 3, -2))
    {
      printf ("Error for exp(-1/16), prec=2, RNDD\n");
      exit (1);
    }
  test_exp (y, x, GMP_RNDZ);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error for exp(-1/16), prec=2, RNDZ\n");
      exit (1);
    }
  mpfr_set_str_binary (x, "0.1E-3");
  test_exp (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error for exp(1/16), prec=2, RNDN\n");
      exit (1);
    }
  test_exp (y, x, GMP_RNDU);
  if (mpfr_cmp_ui_2exp (y, 3, -1))
    {
      printf ("Error for exp(1/16), prec=2, RNDU\n");
      exit (1);
    }

  /* bug reported by Franky Backeljauw, 28 Mar 2003 */
  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);
  mpfr_set_str_binary (x, "1.1101011000111101011110000111010010101001101001110111e28");
  test_exp (y, x, GMP_RNDN);

  mpfr_set_prec (x, 153);
  mpfr_set_prec (z, 153);
  mpfr_set_str_binary (x, "1.1101011000111101011110000111010010101001101001110111e28");
  test_exp (z, x, GMP_RNDN);
  mpfr_prec_round (z, 53, GMP_RNDN);

  if (mpfr_cmp (y, z))
    {
      printf ("Error in mpfr_exp for large argument\n");
      exit (1);
    }

  /* corner cases in mpfr_exp_3 */
  mpfr_set_prec (x, 2);
  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_set_prec (y, 2);
  mpfr_exp_3 (y, x, GMP_RNDN);

  /* Check some little things about overflow detection */
  set_emin (-125);
  set_emax (128);
  mpfr_set_prec (x, 107);
  mpfr_set_prec (y, 107);
  mpfr_set_str_binary (x, "0.11110000000000000000000000000000000000000000000"
                       "0000000000000000000000000000000000000000000000000000"
                       "00000000E4");
  test_exp (y, x, GMP_RNDN);
  if (mpfr_cmp_str (y, "0.11000111100001100110010101111101011010010101010000"
                    "1101110111100010111001011111111000110111001011001101010"
                    "01E22", 2, GMP_RNDN))
    {
      printf ("Special overflow error (1)\n");
      mpfr_dump (y);
      exit (1);
    }

  set_emin (MPFR_EMIN_MIN);
  set_emax (MPFR_EMAX_MAX);

  /* Check for overflow producing a segfault with HUGE exponent */
  mpfr_set_ui  (x, 3, GMP_RNDN);
  mpfr_mul_2ui (x, x, 32, GMP_RNDN);
  test_exp (y, x, GMP_RNDN); /* Can't test return value: May overflow or not*/

  /* Bug due to wrong approximation of (x)/log2 */
  mpfr_set_prec (x, 163);

  mpfr_set_str (x, "-4.28ac8fceeadcda06bb56359017b1c81b85b392e7", 16,
                GMP_RNDN);
  mpfr_exp (x, x, GMP_RNDN);
  if (mpfr_cmp_str (x, "3.fffffffffffffffffffffffffffffffffffffffe8@-2",
                    16, GMP_RNDN))
    {
      printf ("Error for x= -4.28ac8fceeadcda06bb56359017b1c81b85b392e7");
      printf ("expected  3.fffffffffffffffffffffffffffffffffffffffe8@-2");
      printf ("Got       ");
      mpfr_out_str (stdout, 16, 0, x, GMP_RNDN); putchar ('\n');
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

/* check sign of inexact flag */
static void
check_inexact (void)
{
  mpfr_t x, y;
  int inexact;

  mpfr_init2 (x, 53);
  mpfr_init2 (y, 53);

  mpfr_set_str_binary (x,
        "1.0000000000001001000110100100101000001101101011100101e2");
  inexact = test_exp (y, x, GMP_RNDN);
  if (inexact <= 0)
    {
      printf ("Wrong inexact flag (Got %d instead of 1)\n", inexact);
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
check_exp10(void)
{
  mpfr_t x;
  int inexact;

  mpfr_init2 (x, 200);
  mpfr_set_ui(x, 4, GMP_RNDN);

  inexact = mpfr_exp10 (x, x, GMP_RNDN);
  if (mpfr_cmp_ui(x, 10*10*10*10))
    {
      printf ("exp10: Wrong returned value\n");
      exit (1);
    }
  if (inexact != 0)
    {
      printf ("exp10: Wrong inexact flag\n");
      exit (1);
    }

  mpfr_clear (x);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  check_inexact ();
  check_special ();

  test_generic (2, 100, 100);

  compare_exp2_exp3(500);
  check_worst_cases();
  check3("0.0", GMP_RNDU, "1.0");
  check3("-1e-170", GMP_RNDU, "1.0");
  check3("-1e-170", GMP_RNDN, "1.0");
  check3("-8.88024741073346941839e-17", GMP_RNDU, "1.0");
  check3("8.70772839244701057915e-01", GMP_RNDN, "2.38875626491680437269");
  check3("1.0", GMP_RNDN, "2.71828182845904509080");
  check3("-3.42135637628104173534e-07", GMP_RNDZ, "0.999999657864420798958");
  /* worst case for argument reduction, very near from 5*log(2),
     thanks to Jean-Michel Muller  */
  check3("3.4657359027997265421", GMP_RNDN, "32.0");
  check3("3.4657359027997265421", GMP_RNDU, "32.0");
  check3("3.4657359027997265421", GMP_RNDD, "31.999999999999996447");
  check3("2.26523754332090625496e+01", GMP_RNDD, "6.8833785261699581146e9");
  check3("1.31478962104089092122e+01", GMP_RNDZ, "5.12930793917860137299e+05");
  check3("4.25637507920002378103e-01", GMP_RNDU, "1.53056585656161181497e+00");
  check3("6.26551618962329307459e-16", GMP_RNDU, "1.00000000000000066613e+00");
  check3("-3.35589513871216568383e-03",GMP_RNDD, "9.96649729583626853291e-01");
  check3("1.95151388850007272424e+01", GMP_RNDU, "2.98756340674767792225e+08");
  check3("2.45045953503350730784e+01", GMP_RNDN, "4.38743344916128387451e+10");
  check3("2.58165606081678085104e+01", GMP_RNDD, "1.62925781879432281494e+11");
  check3("-2.36539020084338638128e+01",GMP_RNDZ, "5.33630792749924762447e-11");
  check3("2.39211946135858077866e+01", GMP_RNDU, "2.44817704330214385986e+10");
  check3("-2.78190533055889162029e+01",GMP_RNDZ, "8.2858803483596879512e-13");
  check3("2.64028186174889789584e+01", GMP_RNDD, "2.9281844652878973388e11");
  check3("2.92086338843268329413e+01", GMP_RNDZ, "4.8433797301907177734e12");
  check3("-2.46355324071459982349e+01",GMP_RNDZ, "1.9995129297760994791e-11");
  check3("-2.23509444608605427618e+01",GMP_RNDZ, "1.9638492867489702307e-10");
  check3("-2.41175390197331687148e+01",GMP_RNDD, "3.3564940885530624592e-11");
  check3("2.46363885231578088053e+01", GMP_RNDU, "5.0055014282693267822e10");
  check3("111.1263531080090984914932050742208957672119140625", GMP_RNDN, "1.8262572323517295459e48");
  check3("-3.56196340354684821250e+02",GMP_RNDN, "2.0225297096141478156e-155");
  check3("6.59678273772710895173e+02", GMP_RNDU, "3.1234469273830195529e286");
  check3("5.13772529701934331570e+02", GMP_RNDD, "1.3445427121297197752e223");
  check3("3.57430211008718345056e+02", GMP_RNDD, "1.6981197246857298443e155");
  check3("3.82001814471465536371e+02", GMP_RNDU, "7.9667300591087367805e165");
  check3("5.92396038219384422518e+02", GMP_RNDD, "1.880747529554661989e257");
  check3("-5.02678550462488090034e+02",GMP_RNDU, "4.8919201895446217839e-219");
  check3("5.30015757134837031117e+02", GMP_RNDD, "1.5237672861171573939e230");
  check3("5.16239362447650933063e+02", GMP_RNDZ, "1.5845518406744492105e224");
  check3("6.00812634798592370977e-01", GMP_RNDN, "1.823600119339019443");
  check_exp10 ();

  tests_end_mpfr ();
  return 0;
}
