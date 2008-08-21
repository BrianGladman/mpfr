/* tprintf.c -- test file for mpfr_printf and mpfr_vprintf

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

#if defined HAVE_STDARG
#include <stdarg.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#ifdef _MPFR_H_HAVE_INTMAX_T
#include <stdint.h>
#endif

#include "mpfr-test.h"
#define STDOUT_FILENO 1

#if MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0)

/* limit for random precision in random() */
const int prec_max_printf = 5000;
int stdout_redirect;

static void
check (char *fmt, mpfr_t x)
{
  if (mpfr_printf (fmt, x) == -1)
    {
      printf ("Error in mpfr_printf(\"%s\", ...)\n", fmt);

      if (stdout_redirect)
        {
          if ((fflush (stdout) == EOF) || (fclose (stdout) == -1))
            {
              perror ("check");
              exit (1);
            }
        }
      exit (1);
    }
  putchar ('\n');
}

static void
check_vprintf (char *fmt, ...)
{
  va_list ap;

  va_start (ap, fmt);
  if (mpfr_vprintf (fmt, ap) == -1)
    {
      printf ("Error in mpfr_vprintf(\"%s\", ...)\n", fmt);

      if (stdout_redirect)
        {
          if ((fflush (stdout) == EOF) || (fclose (stdout) == -1))
            {
              perror ("check_vprintf");
              exit (1);
            }
        }
      va_end (ap);
      exit (1);
    }
  putchar ('\n');
  va_end (ap);
}

static void
check_special ()
{
  mpfr_t x;

  mpfr_init (x);

  mpfr_set_inf (x, 1);
  check ("%Ra", x);
  check ("%Rb", x);
  check ("%Re", x);
  check ("%Rf", x);
  check ("%Rg", x);
  check_vprintf ("%Ra", x);
  check_vprintf ("%Rb", x);
  check_vprintf ("%Re", x);
  check_vprintf ("%Rf", x);
  check_vprintf ("%Rg", x);

  mpfr_set_inf (x, -1);
  check ("%Ra", x);
  check ("%Rb", x);
  check ("%Re", x);
  check ("%Rf", x);
  check ("%Rg", x);
  check_vprintf ("%Ra", x);
  check_vprintf ("%Rb", x);
  check_vprintf ("%Re", x);
  check_vprintf ("%Rf", x);
  check_vprintf ("%Rg", x);

  mpfr_set_nan (x);
  check ("%Ra", x);
  check ("%Rb", x);
  check ("%Re", x);
  check ("%Rf", x);
  check ("%Rg", x);
  check_vprintf ("%Ra", x);
  check_vprintf ("%Rb", x);
  check_vprintf ("%Re", x);
  check_vprintf ("%Rf", x);
  check_vprintf ("%Rg", x);

  mpfr_clear (x);
}

static void
check_mixed ()
{
  char ch = 'a';
  signed char sch = -1;
  unsigned char uch = 1;
  short sh = -1;
  unsigned short ush = 1;
  int i = -1;
  unsigned int ui = 1;
  long lo = -1;
  unsigned long ulo = 1;
  float f = -1.25;
  double d = -1.25;
  long double ld = -1.25;

  ptrdiff_t p = 1;
  size_t si = 1;

#ifdef HAVE_LONG_LONG
  long long llo = -1;
  unsigned long long ullo = -1;
#endif

#ifdef _MPFR_H_HAVE_INTMAX_T
  intmax_t im = -1;
  uintmax_t uim = 1;
#endif

  mpz_t mpz;
  mpq_t mpq;
  mpf_t mpf;
  mp_rnd_t rnd = GMP_RNDN;

  mpfr_t mpfr;
  mpfr_prec_t prec;

  mpz_init (mpz);
  mpz_set_ui (mpz, ulo);
  mpq_init (mpq);
  mpq_set_si (mpq, lo, ulo);
  mpf_init (mpf);
  mpf_set_q (mpf, mpq);
  mpfr_init (mpfr);
  mpfr_set_f (mpfr, mpf, GMP_RNDN);
  prec = mpfr_get_prec (mpfr);

  check_vprintf ("a. %Ra, b. %hhu, c. %u, d. %lx%hhn", mpfr, uch, ui,
                  ulo, &uch);
  MPFR_ASSERTN (uch == 28);
  check_vprintf ("a. %hhi, b. %Rb, c. %u, d. %li%ln", sch, mpfr, i,
                  lo, &ulo);
  MPFR_ASSERTN (ulo == 37);
  check_vprintf ("a. %hi, b. %*f, c. %Re%hn", ush, 3, f, mpfr, &ush);
  MPFR_ASSERTN (ush == 29);
  check_vprintf ("a. %hi, b. %e, c. %#.2Rf%n", sh, d, mpfr, &i);
  MPFR_ASSERTN (i == 33);
  check_vprintf ("a. %R*A, b. %Fe, c. %i%zn", rnd, mpfr, mpf, si,
                  &si);
  MPFR_ASSERTN (si == 34);
  check_vprintf ("a. %Pu, b. %c, c. %Lf, d. %Zi%Zn", prec, ch, ld,
                  mpz, &mpz);
  MPFR_ASSERTN (mpz_cmp_ui (mpz, 31) == 0);
  check_vprintf ("%% a. %#.0RNg, b. %Qx%Rn, c. %td, d. %p", mpfr, mpq,
                  &mpfr, p, &i);
  MPFR_ASSERTN (mpfr_cmp_ui (mpfr, 16) == 0);

#ifdef HAVE_LONG_LONG
  check_vprintf ("a. %Re, b. %llx%Qn", mpfr, ullo, &mpq);
  MPFR_ASSERTN (mpq_cmp_ui (mpq, 31, 1) == 0);
  check_vprintf ("a. %lli, b. %Rf%Fn", llo, mpfr, &mpf);
  MPFR_ASSERTN (mpf_cmp_ui (mpf, 12) == 0);
  check_vprintf ("a. %qi, b. %Rf%qn", llo, mpfr, &ullo);
  MPFR_ASSERTN (ullo == 12);
#endif

#ifdef _MPFR_H_HAVE_INTMAX_T
  check_vprintf ("a. %*RA, b. %ji%Qn", 10, mpfr, im, &mpq);
  MPFR_ASSERTN (mpq_cmp_ui (mpq, 20, 1) == 0);
  check_vprintf ("a. %.*Re, b. %jx%Fn", 10, mpfr, uim, &mpf);
  MPFR_ASSERTN (mpf_cmp_ui (mpf, 25) == 0);
#endif

  mpfr_clear (mpfr);
  mpf_clear (mpf);
  mpq_clear (mpq);
  mpz_clear (mpz);
}

static void
check_random (int nb_tests)
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
#define FMT_SIZE 13
      char fmt[FMT_SIZE]; /* at most something like "%-+ #0'.*R*f" */
      char *ptr = fmt;

      tests_default_random (x, 256, MPFR_EMIN_MIN, MPFR_EMAX_MAX);
      rnd = (mp_rnd_t) RND_RAND ();

      spec = (int) (randlimb () % 5);
      jmax = (spec == 3 || spec == 4) ? 6 : 5; /* ' flag only with %f or %g */
      /* advantage small precision */
      prec = (randlimb () % 2) ? 10 : prec_max_printf;
      prec = (int) (randlimb () % prec);
      if (spec == 3 && mpfr_get_exp (x) > prec_max_printf)
        /*  change style 'f' to style 'e' when number x is large */
        --spec;

      *ptr++ = '%';
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
      MPFR_ASSERTD (ptr - fmt < FMT_SIZE);

      mpfr_printf ("mpfr_printf(\"%s\", %d, %s, %Re)\n", fmt, prec,
                   mpfr_print_rnd_mode (rnd), x);
      ret = mpfr_printf (fmt, prec, rnd, x);
      if (ret == -1)
        {
          if (spec == 3
              && (MPFR_GET_EXP (x) > INT_MAX || MPFR_GET_EXP (x) < -INT_MAX))
            /* normal failure: x is too large to be output with full precision */
            {
              mpfr_printf ("too large !");
            }
          else
            {
              printf ("Error in mpfr_printf(\"%s\", %d, %s, ...)",
                      fmt, prec, mpfr_print_rnd_mode (rnd));

              if (stdout_redirect)
                {
                  if ((fflush (stdout) == EOF) || (fclose (stdout) == -1))
                    {
                      perror ("check_random");
                      exit (1);
                    }
                }
              exit (1);
            }
        }
      putchar ('\n');
    }

  mpfr_set_emin (old_emin);
  mpfr_set_emax (old_emax);

  mpfr_clear (x);
}

int
main (int argc, char *argv[])
{
  int N;

  tests_start_mpfr ();

  /* with no argument: prints to /dev/null,
     tprintf N: prints N tests to stdout */
  if (argc == 1)
    {
      N = 1000;
      stdout_redirect = -1;
      if (freopen ("/dev/null", "w", stdout) == NULL)
        {
          /* We failed to open this device, try with a dummy file */
          if (freopen ("mpfrtest.txt", "w", stdout) == NULL)
            {
              printf ("Can't open /dev/null or a temporary file\n");
              exit (1);
            }
        }
    }
  else
    {
      stdout_redirect = 0;
      N = atoi (argv[1]);
    }

  check_special ();
  check_mixed ();
  check_random (N);

  if (stdout_redirect)
    {
      if ((fflush (stdout) == EOF) || (fclose (stdout) == -1))
        perror ("main");
    }
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
