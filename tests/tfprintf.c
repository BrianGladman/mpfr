/* tfprintf.c -- test file for mpfr_fprintf and mpfr_vfprintf

Copyright 2008 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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

#ifdef HAVE_STDARG

#include <stdarg.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_LOCALE_H
#include <locale.h>
#endif

#include <gmp.h>

#include "mpfr-test.h"

#if MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0)

const int prec_max_printf = 5000; /* limit for random precision in
                                     random_double() */

static void
check (FILE *fout, char *fmt, mpfr_t x)
{
  if (mpfr_fprintf (fout, fmt, x) == -1)
    {
      mpfr_printf ("Error in mpfr_fprintf(fout, \"%s\", %Re)\n",
                   fmt, x);
      exit (1);
    }
}

static int
check_special (FILE *fout)
{
  mpfr_t x;

  mpfr_init (x);

  mpfr_set_inf (x, 1);
  check (fout, "%Ra", x);
  check (fout, "%Rb", x);
  check (fout, "%Re", x);
  check (fout, "%Rf", x);
  check (fout, "%Rg", x);

  mpfr_set_inf (x, -1);
  check (fout, "%Ra", x);
  check (fout, "%Rb", x);
  check (fout, "%Re", x);
  check (fout, "%Rf", x);
  check (fout, "%Rg", x);

  mpfr_set_nan (x);
  check (fout, "%Ra", x);
  check (fout, "%Rb", x);
  check (fout, "%Re", x);
  check (fout, "%Rf", x);
  check (fout, "%Rg", x);

  mpfr_clear (x);
  return 0;
}

static int
check_random (FILE *fout, int nb_tests)
{
  int i;
  mpfr_t x;
  mp_rnd_t rnd;
  char flag[] =
    {
      '-',
      '+',
      ' ',
      '#',
      '0', /* no ambiguity: first zeros are flag zero*/
      '\''
    };
  char specifier[] =
    {
      'a',
      'b',
      'e',
      'f',
      'g'
    };
  mp_exp_t old_emin, old_emax;

  old_emin = mpfr_get_emin ();
  old_emax = mpfr_get_emax ();

  mpfr_init (x);

  for (i = 0; i < nb_tests; ++i)
    {
      int ret;
      int j, jmax;
      int spec, prec;
      const int fmt_size = 13;
      char fmt[fmt_size]; /* at most something like "%-+ #0'.*R*f" */
      char *ptr = fmt;

      tests_default_random (x, 256, MPFR_EMIN_MIN, MPFR_EMAX_MAX);
      rnd = RND_RAND ();

      *ptr++ = '%';
      spec = (int) (randlimb () % 5);
      jmax = (spec == 3 || spec == 4) ? 6 : 5; /* ' flag only with %f or %g */
      for (j = 0; j < jmax; j++)
        {
          if (randlimb () % 3 == 0)
            *ptr++ = flag[j];
        }
      *ptr++ = '.';
      *ptr++ = '*';
      *ptr++ = 'R';
      *ptr++ = '*';
      *ptr++ = specifier[spec];
      *ptr = '\0';
      MPFR_ASSERTD (ptr - fmt < fmt_size);

      /* advantage small precision */
      prec = (int) (randlimb () % ((randlimb () % 2) ? 10 : prec_max_printf));

      mpfr_fprintf (fout, "mpfr_fprintf(fout, \"%s\", %d, %s, %Re)\n",
                    fmt, prec, mpfr_print_rnd_mode (rnd), x);
      ret = mpfr_fprintf (fout, fmt, prec, rnd, x);
      if (ret == -1)
        {
          if (spec == 3
              && (MPFR_GET_EXP (x) > INT_MAX || MPFR_GET_EXP (x) < -INT_MAX))
            /* normal failure: x is too large to be output with full precision */
            {
              mpfr_fprintf (fout, "too large !");
            }
          else
            {
              mpfr_printf ("Error in mpfr_fprintf(fout, \"%s\", %d, %s, %Re)\n",
                           fmt, prec, mpfr_print_rnd_mode (rnd), x);
              exit (1);
            }
        }
      mpfr_fprintf (fout, "\n");
    }

  mpfr_set_emin (old_emin);
  mpfr_set_emax (old_emax);

  mpfr_clear (x);
  return 0;
}

int
main (int argc, char *argv[])
{
  FILE *fout;
  int N;

  tests_start_mpfr ();

  /* with no argument: prints to /dev/null,
     tfprintf N: prints N tests to stdout */
  if (argc == 1)
    {
      N = 1000;
      fout = fopen ("/dev/null", "w");
      /* If we failed to open this device, try with a dummy file */
      if (fout == NULL)
        {
          fout = fopen ("mpfrtest.txt", "w");

          if (fout == NULL)
            {
              printf ("Can't open /dev/null or a temporary file\n");
              exit (1);
            }
        }

      check_special (fout);
    }
  else
    {
      fout = stdout;
      N = atoi (argv[1]);
      if (fout == NULL)
        {
          printf ("Can't open stdout\n");
          exit (1);
        }
    }

  check_random (fout, N);

  fclose (fout);
  tests_end_mpfr ();
  return 0;
}

#else  /* MPFR_VERSION */

int
main (void)
{
  printf ("Warning! Test disabled for this MPFR version.\n");
  return 0;
}

#endif  /* MPFR_VERSION */

#else  /* HAVE_STDARG */

int
main (void)
{
  /* We have nothing to test. */
  return 0;
}

#endif  /* HAVE_STDARG */
