/* Test file for mpfr_set_float128 and mpfr_get_float128.

Copyright 2012-2016 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef WITH_FPU_CONTROL
#include <fpu_control.h>
#endif

#ifdef MPFR_WANT_FLOAT128

#include "mpfr-test.h"

static void
check_special (void)
{
  __float128 f;
  mpfr_t x;

  mpfr_init2 (x, 113);

  /* check NaN */
  f = 0.0 / 0.0;
  mpfr_set_float128 (x, f, MPFR_RNDN);
  if (mpfr_nan_p (x) == 0)
    {
      printf ("Error in mpfr_set_float128(x, NaN)\n");
      exit (1);
    }
  f = mpfr_get_float128 (x, MPFR_RNDN);
  if (! DOUBLE_ISNAN (f))
    {
      printf ("Error in mpfr_get_float128(NaN)\n");
      printf ("got %f\n", (double) f);
      exit (1);
    }

  /* check +Inf */
  f = 1.0 / 0.0;
  mpfr_set_float128 (x, f, MPFR_RNDN);
  if (mpfr_inf_p (x) == 0 || mpfr_sgn (x) < 0)
    {
      printf ("Error in mpfr_set_float128(x, +Inf)\n");
      exit (1);
    }
  f = mpfr_get_float128 (x, MPFR_RNDN);
  if (f != (1.0 / 0.0))
    {
      printf ("Error in mpfr_get_float128(+Inf)\n");
      exit (1);
    }

  /* check -Inf */
  f = -1.0 / 0.0;
  mpfr_set_float128 (x, f, MPFR_RNDN);
  if (mpfr_inf_p (x) == 0 || mpfr_sgn (x) > 0)
    {
      printf ("Error in mpfr_set_float128(x, -Inf)\n");
      exit (1);
    }
  f = mpfr_get_float128 (x, MPFR_RNDN);
  if (f != (-1.0 / 0.0))
    {
      printf ("Error in mpfr_get_float128(-Inf)\n");
      exit (1);
    }

  /* check +0 */
  f = 0.0;
  mpfr_set_float128 (x, f, MPFR_RNDN);
  if (mpfr_zero_p (x) == 0 || mpfr_sgn (x) < 0)
    {
      printf ("Error in mpfr_set_float128(x, +0)\n");
      exit (1);
    }
  f = mpfr_get_float128 (x, MPFR_RNDN);
  if (f != 0.0 || 1 / f != 1 / 0.0)
    {
      printf ("Error in mpfr_get_float128(+0.0)\n");
      exit (1);
    }

  /* check -0 */
  f = -0.0;
  mpfr_set_float128 (x, f, MPFR_RNDN);
  if (mpfr_zero_p (x) == 0 || mpfr_sgn (x) > 0)
    {
      printf ("Error in mpfr_set_float128(x, -0)\n");
      exit (1);
    }
  f = mpfr_get_float128 (x, MPFR_RNDN);
  if (f != -0.0 || 1 / f != 1 / -0.0)
    {
      printf ("Error in mpfr_get_float128(-0.0)\n");
      exit (1);
    }

  mpfr_clear (x);
}

static void
check_large (void)
{
  __float128 f, e;
  int i;
  mpfr_t x, y;
  mpfr_rnd_t r;

  mpfr_init2 (x, 113);
  mpfr_init2 (y, 113);

  /* check with the largest float128 number 2^16384*(1-2^(-113)) */
  for (f = 1.0, i = 0; i < 113; i++)
    f = f + f;
  f = f - (__float128) 1.0;
  mpfr_set_ui (y, 1, MPFR_RNDN);
  mpfr_mul_2exp (y, y, 113, MPFR_RNDN);
  mpfr_sub_ui (y, y, 1, MPFR_RNDN);
  for (i = 113; i < 16384; i++)
    {
      for (r = 0; r < MPFR_RND_MAX; r++)
        {
          mpfr_set_float128 (x, f, r);
          if (! mpfr_equal_p (x, y))
            {
              printf ("mpfr_set_float128 failed for 2^%d*(1-2^(-113)) rnd=%s\n",
                      i, mpfr_print_rnd_mode (r));
              printf ("got ");
              mpfr_dump (x);
              exit (1);
            }
          e =  mpfr_get_float128 (x, r);
          if (e != f)
            {
              printf ("mpfr_get_float128 failed for 2^%d*(1-2^(-113)) rnd=%s\n",
                      i, mpfr_print_rnd_mode (r));
              exit (1);
            }
        }

      /* check with opposite number */
      f = -f;
      mpfr_neg (y, y, MPFR_RNDN);
      for (r = 0; r < MPFR_RND_MAX; r++)
        {
          mpfr_set_float128 (x, f, r);
          if (! mpfr_equal_p (x, y))
            {
              printf ("mpfr_set_float128 failed for -2^%d*(1-2^(-113)) rnd=%s\n",
                      i, mpfr_print_rnd_mode (r));
              printf ("got ");
              mpfr_dump (x);
              exit (1);
            }
          e =  mpfr_get_float128 (x, r);
          if (e != f)
            {
              printf ("mpfr_get_float128 failed for -2^%d*(1-2^(-113)) rnd=%s\n",
                      i, mpfr_print_rnd_mode (r));
              exit (1);
            }
        }

      f = -f;
      mpfr_neg (y, y, MPFR_RNDN);
      f = f + f;
      mpfr_add (y, y, y, MPFR_RNDN);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
check_small (void)
{
  __float128 f, e;
  int i;
  mpfr_t x, y;
  mpfr_rnd_t r;

  mpfr_init2 (x, 113);
  mpfr_init2 (y, 113);

  /* check with the smallest positive normal float128 number */
  f = 1.0;
  mpfr_set_ui (y, 1, MPFR_RNDN);
  for (i = 0; f != 0.0; i--)
    {
      for (r = 0; r < MPFR_RND_MAX; r++)
        {
          mpfr_set_float128 (x, f, r);
          if (! mpfr_equal_p (x, y))
            {
              printf ("mpfr_set_float128 failed for 2^%d rnd=%s\n", i,
                      mpfr_print_rnd_mode (r));
              printf ("got ");
              mpfr_dump (x);
              exit (1);
            }
          e =  mpfr_get_float128 (x, r);
          if (e != f)
            {
              printf ("mpfr_get_float128 failed for 2^%d rnd=%s\n",
                      i, mpfr_print_rnd_mode (r));
              exit (1);
            }

          /* check with opposite number */
          mpfr_set_float128 (x, -f, r);
          mpfr_neg  (y, y, MPFR_RNDN);
          if (! mpfr_equal_p (x, y))
            {
              printf ("mpfr_set_float128 failed for -2^%d rnd=%s\n", i,
                      mpfr_print_rnd_mode (r));
              printf ("got ");
              mpfr_dump (x);
              exit (1);
            }
          if (e != f)
            {
              printf ("mpfr_get_float128 failed for -2^%d rnd=%s\n",
                      i, mpfr_print_rnd_mode (r));
              exit (1);
            }

          mpfr_neg  (y, y, MPFR_RNDN);
        }
      f =  0.5 * f;
      mpfr_div_2exp (y, y, 1, MPFR_RNDN);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
#ifdef WITH_FPU_CONTROL
  fpu_control_t cw;

  /* cw=895 (0x037f): round to double extended precision
     cw=639 (0x027f): round to double precision
     cw=127 (0x007f): round to single precision */
  if (argc > 1)
    {
      cw = strtol(argv[1], NULL, 0);
      printf ("FPU control word: 0x%x\n", (unsigned int) cw);
      _FPU_SETCW (cw);
    }
#endif

  tests_start_mpfr ();

  check_special ();

  check_large ();

  check_small ();

  tests_end_mpfr ();

  return 0;
}

#else /* MPFR_WANT_FLOAT128 */

/* dummy main to say this test is ignored */
int
main (void)
{
  return 77;
}

#endif /* MPFR_WANT_FLOAT128 */
