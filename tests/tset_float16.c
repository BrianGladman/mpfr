/* Test file for mpfr_set_float16 and mpfr_get_float16.

Copyright 2012-2025 Free Software Foundation, Inc.
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

/* Needed due to the test on MPFR_WANT_FLOAT16 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef MPFR_WANT_FLOAT16

#include "mpfr-test.h"
#include <stdint.h>

typedef union { _Float16 x; uint16_t n; } b16u16;

static void
check_special (mpfr_rnd_t rnd)
{
  b16u16 f;
  mpfr_t x;
  int i;

  mpfr_init2 (x, 11);

  /* check all encodings of NaN */

  for (i = 0x7c01; i < 0x8000; i++)
    {
      f.n = i;
      mpfr_set_float16 (x, f.x, rnd);
      if (! mpfr_nan_p (x))
        {
          printf ("Error in mpfr_set_float16(x, NaN = 0x%x)\n", f.n);
          printf ("got ");
          mpfr_dump (x);
          exit (1);
        }
      f.x = mpfr_get_float16 (x, rnd);
      if (! DOUBLE_ISNAN (f.x))
        {
          printf ("Error in mpfr_get_float16(NaN)\n");
          printf ("got %f\n", (double) f.x);
          exit (1);
        }
      /* also check with sign bit set */
      f.n = 0x8000 | i;
      mpfr_set_float16 (x, f.x, rnd);
      if (! mpfr_nan_p (x))
        {
          printf ("Error in mpfr_set_float16(x, NaN = 0x%x)\n", f.n);
          printf ("got ");
          mpfr_dump (x);
          exit (1);
        }
      f.x = mpfr_get_float16 (x, rnd);
      if (! DOUBLE_ISNAN (f.x))
        {
          printf ("Error in mpfr_get_float16(NaN)\n");
          printf ("got %f\n", (double) f.x);
          exit (1);
        }
    }

  /* check +Inf */
  f.n = 0x7c00;
  mpfr_set_float16 (x, f.x, rnd);
  if (! mpfr_inf_p (x) || MPFR_IS_NEG (x))
    {
      printf ("Error in mpfr_set_float16(x, +Inf)\n");
      printf ("got ");
      mpfr_dump (x);
      exit (1);
    }
  f.x = mpfr_get_float16 (x, rnd);
  if (f.n != 0x7c00)
    {
      printf ("Error in mpfr_get_float16(+Inf)\n");
      printf ("got %f\n", (double) f.x);
      exit (1);
    }

  /* check -Inf */
  f.n = 0xfc00;
  mpfr_set_float16 (x, f.x, rnd);
  if (! mpfr_inf_p (x) || MPFR_IS_POS (x))
    {
      printf ("Error in mpfr_set_float16(x, -Inf)\n");
      printf ("got ");
      mpfr_dump (x);
      exit (1);
    }
  f.x = mpfr_get_float16 (x, rnd);
  if (f.n != 0xfc00)
    {
      printf ("Error in mpfr_get_float16(-Inf)\n");
      printf ("got %f\n", (double) f.x);
      exit (1);
    }

  /* check +0 */
  f.n = 0;
  mpfr_set_float16 (x, f.x, rnd);
  if (! mpfr_zero_p (x) || MPFR_IS_NEG (x))
    {
      printf ("Error in mpfr_set_float16(x, +0)\n");
      printf ("got ");
      mpfr_dump (x);
      exit (1);
    }
  f.x = mpfr_get_float16 (x, rnd);
  if (f.n != 0)
    {
      printf ("Error in mpfr_get_float16(+0.0)\n");
      printf ("got %f\n", (double) f.x);
      exit (1);
    }

  /* check -0 */
  f.n = 0x8000;
  mpfr_set_float16 (x, f.x, rnd);
  if (! mpfr_zero_p (x))
    {
      printf ("Error in mpfr_set_float16(x, -0)\n");
      printf ("got ");
      mpfr_dump (x);
      exit (1);
    }
#if defined(HAVE_SIGNEDZ)
  if (MPFR_IS_POS (x))
    {
      printf ("Error in mpfr_set_float16(x, -0)\n");
      printf ("got ");
      mpfr_dump (x);
      exit (1);
    }
#endif
  f.x = mpfr_get_float16 (x, rnd);
  if (f.n != 0x8000)
    {
      printf ("Error in mpfr_get_float16(-0.0)\n");
      printf ("got %f\n", (double) f.x);
      exit (1);
    }

  mpfr_clear (x);
}

/* check all subnormal and normal numbers */
static void
check_normal (mpfr_rnd_t rnd)
{
  b16u16 f;
  _Float16 g;
  mpfr_t x;
  int i, m, e;

  mpfr_init2 (x, 11);

  m = 1;
  e = -24;
  for (i = 1; i < 0x7c00; i++)
    {
      f.n = i;
      /* Invariant: f.x = m*2^e */
      mpfr_set_float16 (x, f.x, rnd);
      if (mpfr_cmp_ui_2exp (x, m, e) != 0)
        {
          printf ("Error in mpfr_set_float16(x, %a)\n", (float) f.x);
          printf ("got ");
          mpfr_dump (x);
          exit (1);
        }
      g = mpfr_get_float16 (x, rnd);
      if (g != f.x)
        {
          printf ("Error in mpfr_get_float16(%a)\n", (float) f.x);
          printf ("got %a\n", (float) g);
          exit (1);
        }
      /* also check with sign bit set */
      f.n = 0x8000 | i;
      mpfr_set_float16 (x, f.x, rnd);
      if (mpfr_cmp_si_2exp (x, -m, e) != 0)
        {
          printf ("Error in mpfr_set_float16(x, %a)\n", (float) f.x);
          printf ("got ");
          mpfr_dump (x);
          exit (1);
        }
      g = mpfr_get_float16 (x, rnd);
      if (g != f.x)
        {
          printf ("Error in mpfr_get_float16(%a)\n", (float) f.x);
          printf ("got %a\n", (float) g);
          exit (1);
        }

      /* update invariant */
      if (++m == 0x800)
        {
          e++;
          m = m >> 1;
        }
    }

  mpfr_clear (x);
}

int
main (int argc, char *argv[])
{
  int rnd;

  tests_start_mpfr ();

  RND_LOOP (rnd)
    {
      check_special ((mpfr_rnd_t) rnd);
      check_normal ((mpfr_rnd_t) rnd);
    }

  tests_end_mpfr ();

  return 0;
}

#else /* MPFR_WANT_FLOAT16 */

/* dummy main to say this test is ignored */
int
main (void)
{
  return 77;
}

#endif /* MPFR_WANT_FLOAT16 */
