/* Test file for mpfr_compound.

Copyright 2021-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

/* TODO: Generate and test cases that are very close to the overflow
   and underflow thresholds (for all rounding modes), both below and
   above the thresholds.
   Note: This is the goal of the ofuf_thresholds tests in tests/tests.c
   (work in progress); hardcoded tests, e.g. from ofuf_thresholds failures,
   could also be added as ofuf_thresholds may be modified. */

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_compound
#define TWO_ARGS
#define TEST_RANDOM_POS 16 /* the 2nd argument is negative with prob. 16/512 */
#include "tgeneric.c"

static void
check_misc (void)
{
  mpfr_t x, y, z;
  long i;
  const char *sy[] = { "-1.5", "-1", "-0.5", "-0", "+0", "0.5", "1", "1.5",
                       "-Inf", "+Inf", "NaN" };
  int j;
  mpfr_prec_t prec = 2; /* we need at least 2 so that 3/4 is exact */

  mpfr_init2 (x, prec);
  mpfr_init2 (y, prec);
  mpfr_init2 (z, prec);

  /* compound(x,y) = NaN for x < -1, and set invalid exception */
  for (i = 0; i < numberof(sy); i++)
    for (j = 0; j < 2; j++)
      {
        const char *s;
        int inex;
        mpfr_flags_t ex_flags = MPFR_FLAGS_NAN, flags;

        if (j == 0)
          {
            mpfr_set_si (x, -2, MPFR_RNDN);
            s = "-2";
          }
        else
          {
            mpfr_set_inf (x, -1);
            s = "-Inf";
          }
        mpfr_set_str (y, sy[i], 10, MPFR_RNDN);
        mpfr_clear_flags ();
        inex = mpfr_compound (z, x, y, MPFR_RNDN);
        flags = __gmpfr_flags;
        if (! MPFR_IS_NAN (z))
          {
            printf ("Error, compound(%s,%s) should give NaN.\n", s, sy[i]);
            printf ("Got "); mpfr_dump (z);
            exit (1);
          }
        if (inex != 0)
          {
            printf ("Error, compound(%s,%s) should be exact.\n", s, sy[i]);
            printf ("Got inex = %d instead of 0.\n", inex);
            exit (1);
          }
        if (flags != ex_flags)
          {
            printf ("Bad flags for compound(%s,%s).\n", s, sy[i]);
            printf ("Expected flags:");
            flags_out (ex_flags);
            printf ("Got flags:     ");
            flags_out (flags);
            exit (1);
          }
      }

  /* compound(x,0) = 1 for x >= -1 or x = NaN */
  for (i = -2; i <= 2; i++)
    {
      int inex;
      mpfr_flags_t ex_flags = 0, flags;

      if (i == -2)
        mpfr_set_nan (x);
      else if (i == 2)
        mpfr_set_inf (x, 1);
      else
        mpfr_set_si (x, i, MPFR_RNDN);
      mpfr_set_si (y, 0, MPFR_RNDN);
      mpfr_clear_flags ();
      inex = mpfr_compound (z, x, y, MPFR_RNDN);
      flags = __gmpfr_flags;
      if (mpfr_cmp_ui0 (z, 1) != 0)
        {
          printf ("Error, compound(x,0) should give 1 on x = ");
          mpfr_dump (x);
          printf ("got "); mpfr_dump (z);
          exit (1);
        }
      if (inex != 0)
        {
          printf ("Error, compound(x,0) should be exact on x = ");
          mpfr_dump (x);
          printf ("got inex = %d instead of 0.\n", inex);
          exit (1);
        }
      if (flags != ex_flags)
        {
          printf ("Bad flags for compound(x,0) on x = ");
          mpfr_dump (x);
          printf ("Expected flags:");
          flags_out (ex_flags);
          printf ("Got flags:     ");
          flags_out (flags);
          exit (1);
        }
    }

  /* compound(-1,NaN) gives NaN.
     compound(-1,+/-0) = 1.
     compound(-1,y) = +Inf for y < 0, and if y is finite,
     raise the divide-by-zero flag.
     compound(-1,y) = +0 for y > 0.
  */
  mpfr_set_si (x, -1, MPFR_RNDN);
  for (i = 0; i < numberof (sy); i++)
    {
      int inex;
      mpfr_flags_t ex_flags, flags;

      mpfr_set_str (y, sy[i], 10, MPFR_RNDN);
      mpfr_clear_flags ();
      inex = mpfr_compound (z, x, y, MPFR_RNDN);
      flags = __gmpfr_flags;
      if (inex != 0)
        {
          printf ("Error, compound(-1,y) should be exact.\n");
          printf ("got inex = %d instead of 0.\n", inex);
          exit (1);
        }
      if (MPFR_IS_NAN (y))
        {
          if (!MPFR_IS_NAN (z))
            {
              printf ("Error, compound(-1,NaN) should give NaN.\n");
              printf ("Got "); mpfr_dump (z);
              exit (1);
            }
          ex_flags = MPFR_FLAGS_NAN;
        }
      else if (MPFR_IS_ZERO (y))
        {
          if (mpfr_cmp_ui0 (z, 1) != 0)
            {
              printf ("Error, compound(-1,0) should give 1.\n");
              printf ("Got "); mpfr_dump (z);
              exit (1);
            }
          ex_flags = 0;
        }
      else if (MPFR_SIGN (y) < 0)
        {
          if (!MPFR_IS_INF (z) || MPFR_IS_NEG (z))
            {
              printf ("Error, compound(-1,y) should give +Inf with y = ");
              mpfr_dump (y);
              printf ("Got "); mpfr_dump (z);
              exit (1);
            }
          ex_flags = mpfr_inf_p (y) ? 0 : MPFR_FLAGS_DIVBY0;
        }
      else
        {
          MPFR_ASSERTN (MPFR_SIGN (y) > 0);
          if (!MPFR_IS_ZERO (z) || MPFR_IS_NEG (z))
            {
              printf ("Error, compound(-1,y) should give +0 with y = ");
              mpfr_dump (y);
              printf ("Got "); mpfr_dump (z);
              exit (1);
            }
          ex_flags = 0;
        }
      if (flags != ex_flags)
        {
          printf ("Bad flags for compound(-1,y) with y = ");
          mpfr_dump (y);
          printf ("Expected flags:");
          flags_out (ex_flags);
          printf ("Got flags:     ");
          flags_out (flags);
          exit (1);
        }
    }

  /* compound(+/-0,y) = 1 */
  for (i = 0; i < numberof (sy); i++)
    {
      mpfr_set_str (y, sy[i], 10, MPFR_RNDN);
      for (j = 0; j <= 1; j++)
        {
          int inex;
          mpfr_flags_t ex_flags = 0, flags;

          mpfr_set_zero (x, j ? -1 : +1);
          mpfr_clear_flags ();
          inex = mpfr_compound (z, x, y, MPFR_RNDN);
          flags = __gmpfr_flags;
          if (mpfr_cmp_ui0 (z, 1) != 0)
            {
              printf ("Error, compound(x,y) should give 1 on\n  x = ");
              mpfr_dump (x);
              printf ("  y = ");
              mpfr_dump (y);
              printf ("got ");
              mpfr_dump (z);
              exit (1);
            }
          if (inex != 0)
            {
              printf ("Error, compound(x,y) should be exact on\n  x = ");
              mpfr_dump (x);
              printf ("  y = ");
              mpfr_dump (y);
              printf ("got inex = %d instead of 0.\n", inex);
              exit (1);
            }
          if (flags != ex_flags)
            {
              printf ("Bad flags for compound(x,y) on\n  x = ");
              mpfr_dump (x);
              printf ("  y = ");
              mpfr_dump (y);
              printf ("Expected flags:");
              flags_out (ex_flags);
              printf ("Got flags:     ");
              flags_out (flags);
              exit (1);
            }
        }
    }

  /* compound(+Inf,y) = +0 for y < 0
   * compound(+Inf,y) = +Inf for y > 0
   */
  mpfr_set_inf (x, 1);
  for (i = 0; i < numberof (sy); i++)
    {
      int inex;
      mpfr_flags_t ex_flags = 0, flags;
      const char *s;
      int bad_abs_val;

      mpfr_set_str (y, sy[i], 10, MPFR_RNDN);
      if (MPFR_IS_NAN (y) || MPFR_IS_ZERO (y))
        continue;
      /* Now, y < 0 or y > 0. */

      mpfr_clear_flags ();
      inex = mpfr_compound (z, x, y, MPFR_RNDN);
      flags = __gmpfr_flags;

      if (MPFR_IS_NEG (y))
        {
          s = "0";
          bad_abs_val = ! MPFR_IS_ZERO (z);
        }
      else
        {
          s = "Inf";
          bad_abs_val = ! MPFR_IS_INF (z);
        }

      if (bad_abs_val || MPFR_IS_NEG (z))
        {
          printf ("Error, compound(+Inf,%s) should give +%s.\n", sy[i], s);
          printf ("Got "); mpfr_dump (z);
          exit (1);
        }
      if (inex != 0)
        {
          printf ("Error, compound(+Inf,%s) should be exact.\n", sy[i]);
          printf ("Got inex = %d instead of 0.\n", inex);
          exit (1);
        }
      if (flags != ex_flags)
        {
          printf ("Bad flags for compound(+Inf,%s).\n", sy[i]);
          printf ("Expected flags:");
          flags_out (ex_flags);
          printf ("Got flags:     ");
          flags_out (flags);
          exit (1);
        }
    }

  /* compound(NaN,y) gives NaN for y <> 0 */
  mpfr_set_nan (x);
  for (i = 0; i < numberof (sy); i++)
    {
      int inex;
      mpfr_flags_t ex_flags = MPFR_FLAGS_NAN, flags;

      mpfr_set_str (y, sy[i], 10, MPFR_RNDN);
      if (MPFR_IS_ZERO (y))
        continue;

      mpfr_clear_flags ();
      inex = mpfr_compound (z, x, y, MPFR_RNDN);
      flags = __gmpfr_flags;
      if (! MPFR_IS_NAN (z))
        {
          printf ("Error, compound(NaN,%s) should give NaN.\n", sy[i]);
          printf ("Got "); mpfr_dump (z);
          exit (1);
        }
      if (inex != 0)
        {
          printf ("Error, compound(NaN,%s) should be exact.\n", sy[i]);
          printf ("Got inex = %d instead of 0.\n", inex);
          exit (1);
        }
      if (flags != ex_flags)
        {
          printf ("Bad flags for compound(NaN,%s).\n", sy[i]);
          printf ("Expected flags:");
          flags_out (ex_flags);
          printf ("Got flags:     ");
          flags_out (flags);
          exit (1);
        }
    }

  /* hard-coded test: x is the 32-bit nearest approximation of 17/42 */
  mpfr_set_prec (x, 32);
  mpfr_set_prec (z, 32);
  mpfr_set_ui_2exp (x, 3476878287UL, -33, MPFR_RNDN);
  mpfr_set_si (y, 12, MPFR_RNDN);
  mpfr_compound (z, x, y, MPFR_RNDN);
  mpfr_set_ui_2exp (x, 1981447393UL, -25, MPFR_RNDN);
  if (!mpfr_equal_p (z, x))
    {
      printf ("Error for compound(3476878287/2^33,12)\n");
      printf ("expected "); mpfr_dump (x);
      printf ("got      "); mpfr_dump (z);
      exit (1);
    }

  /* test for negative integer y */
  i = -1;
  mpfr_set_prec (y, 32);
  while (1)
    {
      /* i has the form -(2^k-1) */
      mpfr_set_si_2exp (x, -1, -1, MPFR_RNDN); /* x = -0.5 */
      mpfr_set_si (y, i, MPFR_RNDN);
      mpfr_compound (z, x, y, MPFR_RNDN);
      mpfr_set_ui_2exp (x, 1, -i, MPFR_RNDN);
      if (!mpfr_equal_p (z, x))
        {
          printf ("Error for compound(-0.5,%ld)\n", i);
          printf ("expected "); mpfr_dump (x);
          printf ("got      "); mpfr_dump (z);
          exit (1);
        }
      if (i == -2147483647) /* largest possible value on 32-bit machine */
        break;
      i = 2 * i - 1;
    }

  /* The "#if" makes sure that 64-bit constants are supported, avoiding
     a compilation failure. The "if" makes sure that the constant is
     representable in a long (this would not be the case with 32-bit
     unsigned long and 64-bit limb). */
#if GMP_NUMB_BITS >= 64 || MPFR_PREC_BITS >= 64
  /* another test for large negative integer y */
  if (4994322635099777669 <= LONG_MAX)
    {
      i = -4994322635099777669;
      mpfr_set_prec (y, 63);
      mpfr_set_ui (x, 1, MPFR_RNDN);
      mpfr_set_si (y, i, MPFR_RNDN);
      mpfr_compound (z, x, y, MPFR_RNDN);
      mpfr_set_si (x, 1, MPFR_RNDN);
      mpfr_mul_2si (x, x, i, MPFR_RNDN);
      if (!mpfr_equal_p (z, x))
        {
          printf ("Error for compound(1,%ld)\n", i);
          printf ("expected "); mpfr_dump (x);
          printf ("got      "); mpfr_dump (z);
          exit (1);
        }
    }
#endif

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

/* Failure with mpfr_compound_si from 2021-02-15 (e.g. in MPFR 4.2.0)
   due to incorrect underflow detection. */
static void
bug_20230206 (void)
{
  if (MPFR_PREC_MIN == 1)
    {
      mpfr_t x, y1, y2, z;
      int inex1, inex2;
      mpfr_flags_t flags1, flags2;
#if MPFR_PREC_BITS >= 64
      mpfr_exp_t emin;
#endif

      mpfr_inits2 (1, x, y1, y2, (mpfr_ptr) 0);
      mpfr_init2 (z, 63);
      mpfr_set_ui_2exp (x, 1, -1, MPFR_RNDN);  /* x = 1/2 */

      /* This first test is useful mainly for a 32-bit mpfr_exp_t type
         (no failure with a 64-bit mpfr_exp_t type since the underflow
         threshold in the extended exponent range is much lower). */

      mpfr_set_ui_2exp (y1, 1, -1072124363, MPFR_RNDN);
      inex1 = -1;
      flags1 = MPFR_FLAGS_INEXACT;
      mpfr_clear_flags ();
      /* -1832808704 ~= -2^30 / log2(3/2) */
      mpfr_set_si (z, -1832808704, MPFR_RNDN);
      inex2 = mpfr_compound (y2, x, z, MPFR_RNDN);
      flags2 = __gmpfr_flags;
      if (!(mpfr_equal_p (y1, y2) &&
            SAME_SIGN (inex1, inex2) &&
            flags1 == flags2))
        {
          printf ("Error in bug_20230206 (1):\n");
          printf ("Expected ");
          mpfr_dump (y1);
          printf ("  with inex = %d, flags =", inex1);
          flags_out (flags1);
          printf ("Got      ");
          mpfr_dump (y2);
          printf ("  with inex = %d, flags =", inex2);
          flags_out (flags2);
          exit (1);
        }

      /* This second test is for a 64-bit mpfr_exp_t type
         (it is disabled with a 32-bit mpfr_exp_t type). */

      /* The "#if" makes sure that 64-bit constants are supported, avoiding
         a compilation failure. The "if" makes sure that the constant is
         representable in a long (this would not be the case with 32-bit
         unsigned long and 64-bit limb). It also ensures that mpfr_exp_t
         has at least 64 bits. */
#if MPFR_PREC_BITS >= 64
      emin = mpfr_get_emin ();
      set_emin (MPFR_EMIN_MIN);
      mpfr_set_ui_2exp (y1, 1, -4611686018427366846, MPFR_RNDN);
      inex1 = 1;
      flags1 = MPFR_FLAGS_INEXACT;
      mpfr_clear_flags ();
      /* -7883729320669216768 ~= -2^62 / log2(3/2) */
      mpfr_set_si (z, -7883729320669216768, MPFR_RNDN);
      inex2 = mpfr_compound (y2, x, z, MPFR_RNDN);
      flags2 = __gmpfr_flags;
      if (!(mpfr_equal_p (y1, y2) &&
            SAME_SIGN (inex1, inex2) &&
            flags1 == flags2))
        {
          printf ("Error in bug_20230206 (2):\n");
          printf ("Expected ");
          mpfr_dump (y1);
          printf ("  with inex = %d, flags =", inex1);
          flags_out (flags1);
          printf ("Got      ");
          mpfr_dump (y2);
          printf ("  with inex = %d, flags =", inex2);
          flags_out (flags2);
          exit (1);
        }
      set_emin (emin);
#endif

      mpfr_clears (x, y1, y2, z, (mpfr_ptr) 0);
    }
}

/* Reported by Patrick Pelissier on 2023-02-11
   (tgeneric_ui.c with GMP_CHECK_RANDOMIZE=1412991715).
   On a 32-bit host, one gets Inf (overflow) instead of 0.1E1071805703.
*/
static void
bug_20230211 (void)
{
  mpfr_t x, y1, y2, z;
  int inex1, inex2;
  mpfr_flags_t flags1, flags2;

  mpfr_inits2 (1, x, y1, y2, (mpfr_ptr) 0);
  mpfr_init2 (z, 31);
  mpfr_set_ui_2exp (x, 1, -1, MPFR_RNDN);  /* x = 1/2 */
  mpfr_set_ui_2exp (y1, 1, 1071805702, MPFR_RNDN);
  inex1 = 1;
  flags1 = MPFR_FLAGS_INEXACT;
  mpfr_clear_flags ();
  mpfr_set_ui (z, 1832263949, MPFR_RNDN);
  inex2 = mpfr_compound (y2, x, z, MPFR_RNDN);
  flags2 = __gmpfr_flags;
  if (!(mpfr_equal_p (y1, y2) &&
        SAME_SIGN (inex1, inex2) &&
        flags1 == flags2))
    {
      printf ("Error in bug_20230211:\n");
      printf ("Expected ");
      mpfr_dump (y1);
      printf ("  with inex = %d, flags =", inex1);
      flags_out (flags1);
      printf ("Got      ");
      mpfr_dump (y2);
      printf ("  with inex = %d, flags =", inex2);
      flags_out (flags2);
      exit (1);
    }
  mpfr_clears (x, y1, y2, z, (mpfr_ptr) 0);
}

/* Integer overflow with compound.c d04caeae04c6a83276916c4fbac1fe9b0cec3c8b
   (2023-02-23) on "n * (kx - 1) + 1". Note: if the only effect is just a
   random value, this probably doesn't affect the result (one might enter
   the "if" while one shouldn't, but the real check is done inside the "if").
   This test fails if -fsanitize=undefined -fno-sanitize-recover is used or
   if the processor emits a signal in case of integer overflow.
   This test has been made obsolete by the "kx < ex" condition
   in 2cb3123891dd46fe0258d4aec7f8655b8ec69aaf. */
static void
bug_20230517 (void)
{
  mpfr_exp_t old_emax;
  mpfr_t x, z;

  old_emax = mpfr_get_emax ();
  set_emax (MPFR_EMAX_MAX);

  mpfr_init2 (x, 123456);
  mpfr_init2 (z, 64);
  mpfr_set_ui (x, 65536, MPFR_RNDN);
  mpfr_nextabove (x);
  mpfr_set_ui (z, LONG_MAX >> 16, MPFR_RNDN);
  mpfr_compound (x, x, z, MPFR_RNDN);
  mpfr_clear (x);
  mpfr_clear (z);

  set_emax (old_emax);
}

int
main (void)
{
  tests_start_mpfr ();

  check_misc ();
  bug_20230206 ();
  bug_20230211 ();
  bug_20230517 ();

  test_generic (MPFR_PREC_MIN, 100, 100);

  tests_end_mpfr ();
  return 0;
}
