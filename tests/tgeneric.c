/* Generic test file for functions with one or two arguments (the second being
   either mpfr_t or double).

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 Free Software Foundation, Inc.
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

/* define TWO_ARGS for two-argument functions like mpfr_pow
   define DOUBLE_ARG1 or DOUBLE_ARG2 for function with a double operand in
   first or second place like sub_d or d_sub */

#ifndef TEST_RANDOM_POS
/* For the random function: one number on two is negative. */
#define TEST_RANDOM_POS 256
#endif

#ifndef TEST_RANDOM_POS2
/* For the random function: one number on two is negative. */
#define TEST_RANDOM_POS2 256
#endif

#ifndef TEST_RANDOM_EMIN
#define TEST_RANDOM_EMIN -256
#endif

#ifndef TEST_RANDOM_EMAX
#define TEST_RANDOM_EMAX 255
#endif

#define TGENERIC_FAIL(S, X, U)                                          \
  do                                                                    \
    {                                                                   \
      printf ("%s\nx = ", (S));                                         \
      mpfr_out_str (stdout, 2, 0, (X), GMP_RNDN);                       \
      printf ("\n");                                                    \
      if (U)                                                            \
        {                                                               \
          printf ("u = ");                                              \
          mpfr_out_str (stdout, 2, 0, (U), GMP_RNDN);                   \
          printf ("\n");                                                \
        }                                                               \
      printf ("yprec = %u, rnd_mode = %s, inexact = %d, flags = %u\n",  \
              (unsigned int) yprec, mpfr_print_rnd_mode (rnd), compare, \
              (unsigned int) __gmpfr_flags);                            \
      exit (1);                                                         \
    }                                                                   \
  while (0)

#undef TGENERIC_CHECK
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
#define TGENERIC_CHECK(S, EXPR) \
  do if (!(EXPR)) TGENERIC_FAIL (S, x, u); while (0)
#else
#define TGENERIC_CHECK(S, EXPR) \
  do if (!(EXPR)) TGENERIC_FAIL (S, x, 0); while (0)
#endif

#if DEBUG_TGENERIC
#define STR(F) #F
#define TGENERIC_IAUX(F,P,X,U)                                          \
  do                                                                    \
    {                                                                   \
      printf ("tgeneric: testing function " STR(F)                      \
              ", %s, target prec = %lu\nx = ",                          \
              mpfr_print_rnd_mode (rnd), (unsigned long) (P));          \
      mpfr_out_str (stdout, 2, 0, (X), GMP_RNDN);                       \
      printf ("\n");                                                    \
      if (U)                                                            \
        {                                                               \
          printf ("u = ");                                              \
          mpfr_out_str (stdout, 2, 0, (U), GMP_RNDN);                   \
          printf ("\n");                                                \
        }                                                               \
    }                                                                   \
  while (0)
#undef TGENERIC_INFO
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
#define TGENERIC_INFO(F,P) TGENERIC_IAUX(F,P,x,u)
#else
#define TGENERIC_INFO(F,P) TGENERIC_IAUX(F,P,x,0)
#endif
#endif

/* For some functions (for example cos), the argument reduction is too
   expensive when using mpfr_get_emax(). Then simply define REDUCE_EMAX
   to some reasonable value before including tgeneric.c. */
#ifndef REDUCE_EMAX
#define REDUCE_EMAX mpfr_get_emax ()
#endif

static void
test_generic (mp_prec_t p0, mp_prec_t p1, unsigned int N)
{
  mp_prec_t prec, xprec, yprec;
  mpfr_t x, y, z, t;
#ifdef TWO_ARGS
  mpfr_t u;
#elif defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
  mpfr_t u;
  double d;
#endif
  mp_rnd_t rnd;
  int inexact, compare, compare2;
  unsigned int n;
  unsigned long ctrt = 0, ctrn = 0;
  mp_exp_t old_emin, old_emax;

  old_emin = mpfr_get_emin ();
  old_emax = mpfr_get_emax ();

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);
  mpfr_init (t);
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
  mpfr_init (u);
#endif

  /* generic test */
  for (prec = p0; prec <= p1; prec++)
    {
      mpfr_set_prec (z, prec);
      mpfr_set_prec (t, prec);
      yprec = prec + 10;
      mpfr_set_prec (y, yprec);

      /* Note: in precision p1, we test 4 special cases. */
      for (n = 0; n < (prec == p1 ? N + 4 : N); n++)
        {
          xprec = prec;
          if (randlimb () & 1)
            {
              xprec *= (double) randlimb () / MP_LIMB_T_MAX;
              if (xprec < MPFR_PREC_MIN)
                xprec = MPFR_PREC_MIN;
            }
          mpfr_set_prec (x, xprec);
#ifdef TWO_ARGS
          mpfr_set_prec (u, xprec);
#elif defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
          mpfr_set_prec (u, IEEE_DBL_MANT_DIG);
#endif

          if (n > 3 || prec < p1)
            {
#if defined(RAND_FUNCTION)
              RAND_FUNCTION (x);
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
              RAND_FUNCTION (u);
#endif
#else
              tests_default_random (x, TEST_RANDOM_POS,
                                    TEST_RANDOM_EMIN, TEST_RANDOM_EMAX);
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
              tests_default_random (u, TEST_RANDOM_POS2,
                                    TEST_RANDOM_EMIN, TEST_RANDOM_EMAX);
#endif
#endif
            }
          else if (n <= 1)
            {
              /* Special cases tested in precision p1 if n <= 1. */
              mpfr_set_si (x, n == 0 ? 1 : -1, GMP_RNDN);
              mpfr_set_exp (x, mpfr_get_emin ());
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
              mpfr_set_si (u, randlimb () % 2 == 0 ? 1 : -1, GMP_RNDN);
              mpfr_set_exp (u, mpfr_get_emin ());
#endif
            }
          else  /* 2 <= n <= 3 */
            {
              /* Special cases tested in precision p1 if 2 <= n <= 3. */
              if (getenv ("MPFR_CHECK_MAX") == NULL)
                goto next_n;
              mpfr_set_si (x, n == 0 ? 1 : -1, GMP_RNDN);
              mpfr_setmax (x, REDUCE_EMAX);
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
              mpfr_set_si (u, randlimb () % 2 == 0 ? 1 : -1, GMP_RNDN);
              mpfr_setmax (u, mpfr_get_emax ());
#endif
            }

          rnd = (mp_rnd_t) RND_RAND ();
          mpfr_clear_flags ();
#ifdef DEBUG_TGENERIC
          TGENERIC_INFO (TEST_FUNCTION, MPFR_PREC (y));
#endif
#if defined(TWO_ARGS)
          compare = TEST_FUNCTION (y, x, u, rnd);
#elif defined(DOUBLE_ARG1)
          d = mpfr_get_d (u, rnd);
          compare = TEST_FUNCTION (y, d, x, rnd);
#elif defined(DOUBLE_ARG2)
          d = mpfr_get_d (u, rnd);
          compare = TEST_FUNCTION (y, x, d, rnd);
#else
          compare = TEST_FUNCTION (y, x, rnd);
#endif
          TGENERIC_CHECK ("Bad inexact flag",
                          (compare != 0) ^ (mpfr_inexflag_p () == 0));
          ctrt++;
          if (MPFR_IS_SINGULAR (y))
            {
              if (MPFR_IS_NAN (y) || mpfr_nanflag_p ())
                TGENERIC_CHECK ("Bad NaN flag",
                                MPFR_IS_NAN (y) && mpfr_nanflag_p ());
              else if (MPFR_IS_INF (y))
                TGENERIC_CHECK ("Bad overflow flag",
                                (compare != 0) ^ (mpfr_overflow_p () == 0));
              else if (MPFR_IS_ZERO (y))
                TGENERIC_CHECK ("Bad underflow flag",
                                (compare != 0) ^ (mpfr_underflow_p () == 0));
            }
          else if (mpfr_overflow_p ())
            {
              TGENERIC_CHECK ("Bad compare value (overflow)", compare != 0);
              mpfr_nexttoinf (y);
              TGENERIC_CHECK ("Should have been max MPFR number",
                              MPFR_IS_INF (y));
            }
          else if (mpfr_underflow_p ())
            {
              TGENERIC_CHECK ("Bad compare value (underflow)", compare != 0);
              mpfr_nexttozero (y);
              TGENERIC_CHECK ("Should have been min MPFR number",
                              MPFR_IS_ZERO (y));
            }
          else if (mpfr_can_round (y, yprec, rnd, rnd, prec))
            {
              ctrn++;
              mpfr_set (t, y, rnd);
              /* Risk of failures are known when some flags are already set
                 before the function call. Do not set the erange flag, as
                 it will remain set after the function call and no checks
                 are performed in such a case (see the mpfr_erangeflag_p
                 test below). */
              if (randlimb () & 1)
                __gmpfr_flags = MPFR_FLAGS_ALL ^ MPFR_FLAGS_ERANGE;
#ifdef DEBUG_TGENERIC
              TGENERIC_INFO (TEST_FUNCTION, MPFR_PREC (z));
#endif
#if defined(TWO_ARGS)
              inexact = TEST_FUNCTION (z, x, u, rnd);
#elif defined(DOUBLE_ARG1)
              inexact = TEST_FUNCTION (z, d, x, rnd);
#elif defined(DOUBLE_ARG2)
              inexact = TEST_FUNCTION (z, x, d, rnd);
#else
              inexact = TEST_FUNCTION (z, x, rnd);
#endif
              if (mpfr_erangeflag_p ())
                goto next_n;
              if (mpfr_nan_p (z) || mpfr_cmp (t, z) != 0)
                {
                  printf ("results differ for x=");
                  mpfr_out_str (stdout, 2, xprec, x, GMP_RNDN);
#ifdef TWO_ARGS
                  printf ("\nu=");
                  mpfr_out_str (stdout, 2, xprec, u, GMP_RNDN);
#elif defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
                  printf ("\nu=");
                  mpfr_out_str (stdout, 2, IEEE_DBL_MANT_DIG, u, GMP_RNDN);
#endif
                  printf (" prec=%u rnd_mode=%s\n", (unsigned) prec,
                          mpfr_print_rnd_mode (rnd));
                  printf ("got      ");
                  mpfr_out_str (stdout, 2, prec, z, GMP_RNDN);
                  puts ("");
                  printf ("expected ");
                  mpfr_out_str (stdout, 2, prec, t, GMP_RNDN);
                  puts ("");
                  printf ("approx   ");
                  mpfr_print_binary (y);
                  puts ("");
                  exit (1);
                }
              compare2 = mpfr_cmp (t, y);
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
                  printf ("Wrong inexact flag for rnd=%s: expected %d, got %d"
                          "\n", mpfr_print_rnd_mode (rnd), compare, inexact);
                  printf ("x="); mpfr_print_binary (x); puts ("");
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
                  printf ("u="); mpfr_print_binary (u); puts ("");
#endif
                  printf ("y="); mpfr_print_binary (y); puts ("");
                  printf ("t="); mpfr_print_binary (t); puts ("");
                  exit (1);
                }
            }

        next_n:
          /* In case the exponent range has been changed by
             tests_default_random()... */
          mpfr_set_emin (old_emin);
          mpfr_set_emax (old_emax);
        }
    }

#ifndef TGENERIC_NOWARNING
  if (3 * ctrn < 2 * ctrt)
    printf ("Warning! Too few normal cases in generic tests (%lu / %lu)\n",
            ctrn, ctrt);
#endif

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
#if defined(TWO_ARGS) || defined(DOUBLE_ARG1) || defined(DOUBLE_ARG2)
  mpfr_clear (u);
#endif
}

#undef TEST_RANDOM_POS
#undef TEST_RANDOM_POS2
#undef TEST_RANDOM_EMIN
#undef TEST_RANDOM_EMAX
#undef RAND_FUNCTION
#undef TWO_ARGS
#undef TWO_ARGS_UI
#undef TEST_FUNCTION
#undef test_generic
