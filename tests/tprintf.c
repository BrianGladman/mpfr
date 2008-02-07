/* mpfr_tprintf -- test file for mpfr_sprintf mpfr_vsprintf

Copyright 2007, 2008 Free Software Foundation, Inc.
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
const int buf_size = 1024;

const char pinf_str[] = "inf";
const char pinf_uc_str[] = "INF";
const char minf_str[] = "-inf";
const char minf_uc_str[] = "-INF";
const char nan_str[] = "nan";
const char nan_uc_str[] = "NAN";

/* compare expected string with the return of mpfr_sprintf(fmt, x)*/
static int
check_sprintf (const char *expected, const char *fmt, mpfr_srcptr x)
{
  int n;
  char buffer[buf_size];

  n = mpfr_sprintf (buffer, fmt, x);
  if (strcmp (buffer, expected) != 0)
    {
      printf ("Error in mpfr_sprintf (s, \"%s\", x);\n", fmt);
      printf ("expected: \"%s\"\ngot:      \"%s\"\n", expected, buffer);

      exit (1);
    }
  return n;
}

/* compare expected string with the return of mpfr_vsprintf(fmt, ...)*/
static int
check_vsprintf (const char *expected, const char *fmt, ...)
{
  int n;
  char buffer[buf_size];
  va_list ap;
  va_start (ap, fmt);

  n = mpfr_vsprintf (buffer, fmt, ap);
  if (strcmp (buffer, expected) != 0)
    {
      printf ("Error in mpfr_vsprintf (s, \"%s\", ...);\n", fmt);
      printf ("expected: \"%s\"\ngot:      \"%s\"\n", expected, buffer);

      va_end (ap);
      exit (1);
    }

  va_end (ap);
  return n;
}

static int
decimal (void)
{
  mpfr_t x;
  mpfr_t z;
  mpfr_init (z);
  mpfr_init2 (x, 128);

  /* special numbers */
  mpfr_set_inf (x, 1);
  check_sprintf (pinf_str, "%Re", x);
  check_sprintf (pinf_uc_str, "%RE", x);
  check_sprintf (pinf_str, "%Rf", x);
  check_sprintf (pinf_uc_str, "%RF", x);
  check_sprintf (pinf_str, "%Rg", x);
  check_sprintf (pinf_uc_str, "%RG", x);

  mpfr_set_inf (x, -1);
  check_sprintf (minf_str, "%Re", x);
  check_sprintf (minf_uc_str, "%RE", x);
  check_sprintf (minf_str, "%Rf", x);
  check_sprintf (minf_uc_str, "%RF", x);
  check_sprintf (minf_str, "%Rg", x);
  check_sprintf (minf_uc_str, "%RG", x);

  mpfr_set_nan (x);
  check_sprintf (nan_str, "%Re", x);
  check_sprintf (nan_uc_str, "%RE", x);
  check_sprintf (nan_str, "%Rf", x);
  check_sprintf (nan_uc_str, "%RF", x);
  check_sprintf (nan_str, "%Rg", x);
  check_sprintf (nan_uc_str, "%RG", x);

  /* positive numbers */
  mpfr_set_ui (x, 279296875, GMP_RNDN);
  mpfr_div_ui (x, x, 1000000000, GMP_RNDN);
  mpfr_add_ui (x, x, 1899347461, GMP_RNDN);
  mpfr_div_ui (x, x, 100, GMP_RNDN); /* x = 18993474.61279296875 */
  mpfr_set_ui (z, 0, GMP_RNDD);

  /* simplest case right justified */
  check_sprintf ("      1.899347461279296875e+07", "%30Re", x);
  check_sprintf ("                         2e+07", "%30.0Re", x);
  check_sprintf ("          18993474.61279296875", "%30Rf", x);
  check_sprintf ("              18993474.6127930", "%30.7Rf", x);
  check_sprintf ("                   1.89935e+07", "%30Rg", x);
  check_sprintf ("                         2e+07", "%30.0Rg", x);
  check_sprintf ("          18993474.61279296875", "%30.19Rg", x);
  check_sprintf ("                         0e+00", "%30.0Re", z);
  check_sprintf ("                             0", "%30.0Rf", z);
  check_sprintf ("                        0.0000", "%30.4Rf", z);
  check_sprintf ("                             0", "%30.0Rg", z);
  check_sprintf ("                        0.0000", "%30.4Rg", z);
  /* sign or space, pad with leading zeros */
  check_sprintf (" 000001.899347461279296875E+07", "% 030RE", x);
  check_sprintf (" 0000000000000000001.89935E+07", "% 030RG", x);
  check_sprintf (" 0000000000000000000000002E+07", "% 030.0RE", x);
  check_sprintf (" 0000000000000000000000000E+00", "% 030.0RE", z);
  check_sprintf (" 00000000000000000000000000000", "% 030.0RF", z);
  /* sign + or -, left justified */
  check_sprintf ("+1.899347461279296875e+07     ", "%+-30Re", x);
  check_sprintf ("+2e+07                        ", "%+-30.0Re", x);
  check_sprintf ("+0e+00                        ", "%+-30.0Re", z);
  check_sprintf ("+0                            ", "%+-30.0Rf", z);
  /* decimal point, left justified, precision and rounding parameter */
  check_vsprintf ("1.9E+07   ", "%#-10.*R*E", 1, GMP_RNDN, x);
  check_vsprintf ("2.E+07    ", "%#-10.*R*E", 0, GMP_RNDN, x);
  check_vsprintf ("2.E+07    ", "%#-10.*R*G", 0, GMP_RNDN, x);
  check_vsprintf ("0.E+00    ", "%#-10.*R*E", 0, GMP_RNDN, z);
  check_vsprintf ("0.        ", "%#-10.*R*F", 0, GMP_RNDN, z);
  check_vsprintf ("0.        ", "%#-10.*R*G", 0, GMP_RNDN, z);
  /* sign or space */
  check_sprintf (" 1.899e+07", "% .3RNe", x);
  check_sprintf (" 2e+07",     "% .0RNe", x);
  /* sign + or -, decimal point, pad with leading zeros */
  check_sprintf ("+0001.8E+07", "%0+#11.1RZE", x);
  check_sprintf ("+00001.E+07", "%0+#11.0RZE", x);
  check_sprintf ("+0000.0E+00", "%0+#11.1RZE", z);
  check_sprintf ("+00000000.0", "%0+#11.1RZF", z);
  /* pad with leading zero */
  check_sprintf ("0000001.899347461279296875e+07", "%030RDe", x);
  check_sprintf ("00000000000000000000000001e+07", "%030.0RDe", x);
  /* sign or space, decimal point, left justified */
  check_sprintf (" 1.8E+07   ", "%- #11.1RDE", x);
  check_sprintf (" 1.E+07    ", "%- #11.0RDE", x);

  /* negative numbers */
  mpfr_mul_si (x, x, -1, GMP_RNDD);
  mpfr_mul_si (z, z, -1, GMP_RNDD);

  /* sign + or - */
  check_sprintf ("  -1.8e+07", "%+10.1RUe", x);
  check_sprintf ("    -1e+07", "%+10.0RUe", x);
  check_sprintf ("    -0e+00", "%+10.0RUe", z);
  check_sprintf ("        -0", "%+10.0RUf", z);


  /* neighborhood of 1 */
  mpfr_set_ui (x, 9999, GMP_RNDN);
  mpfr_div_ui (x, x, 10000, GMP_RNDN); /* x=0.9999 */
  check_sprintf ("1E+00     ", "%-10.0RE", x);
  check_sprintf ("1.0E+00   ", "%-10.1RE", x);
  check_sprintf ("9.9990E-01", "%-10.4RE", x);
  check_sprintf ("1.0       ", "%-10.1RF", x);
  check_sprintf ("0.9999    ", "%-10.4RF", x);
  check_sprintf ("1         ", "%-10.0RG", x);
  check_sprintf ("1         ", "%-10.1RG", x);
  check_sprintf ("0.9999    ", "%-10.4RG", x);
  check_sprintf ("1.        ", "%-#10.0RG", x);
  check_sprintf ("1.        ", "%-#10.1RG", x);
  check_sprintf ("1.0       ", "%-#10.2RG", x);
  check_sprintf ("0.9999    ", "%-#10.4RG", x);

  /* multiple of 10 */
  mpfr_set_ui (x, 17, GMP_RNDN);
  mpfr_exp10 (x, x, GMP_RNDN); /* x=1e17 */
  check_sprintf ("1e+17", "%Re", x);
  check_sprintf ("1.000e+17", "%.3Re", x);
  check_sprintf ("100000000000000000", "%Rf", x);
  check_sprintf ("100000000000000000.0", "%.1Rf", x);

  mpfr_ui_div (x, 1, x, GMP_RNDN); /* x=1e-17 */
  check_sprintf ("1e-17", "%Re", x);
  check_sprintf ("0.00000000000000001", "%Rf", x);
  check_sprintf ("1e-17", "%Rg", x);
  check_sprintf ("0.0", "%.1RDf", x);
  check_sprintf ("0.1", "%.1RUf", x);
  check_sprintf ("0", "%.0RDf", x);
  check_sprintf ("1", "%.0RUf", x);

  /* check rounding mode */
  mpfr_set_ui (x, 76, GMP_RNDN);
  mpfr_div_ui (x, x, 10000, GMP_RNDN); /* x= 0.0076 */
  check_sprintf ("0.007", "%.3RDF", x);
  check_sprintf ("0.007", "%.3RZF", x);
  check_sprintf ("0.008", "%.3RF", x);
  check_sprintf ("0.008", "%.3RUF", x);

  /* limit test for the choice beetwen %f-style and %g-style */
  mpfr_set_ui (x, 999, GMP_RNDN);
  mpfr_div_ui (x, x, 10000000, GMP_RNDN); /* x=0.0000999 */
  check_sprintf ("0.0001", "%.0Rg", x);
  check_sprintf ("9e-05", "%.0RDg", x);
  check_sprintf ("0.0001", "%.2Rg", x);

  mpfr_set_si_2exp (x, -1, -15, GMP_RNDN); /* x=-2^-15 */
  check_sprintf ("-3.0517578125e-05", "%.300Rg", x);

  mpfr_clears (x, z, (void *)0);
  return 0;
}

static int
hexadecimal (void)
{
  mpfr_t x, z;
  mpfr_inits2 (64, x, z, (void *)0);

  /* special */
  mpfr_set_inf (x, 1);
  check_sprintf (pinf_str, "%Ra", x);
  check_sprintf (pinf_uc_str, "%RA", x);

  mpfr_set_inf (x, -1);
  check_sprintf (minf_str, "%Ra", x);
  check_sprintf (minf_uc_str, "%RA", x);

  mpfr_set_nan (x);
  check_sprintf (nan_str, "%Ra", x);
  check_sprintf (nan_uc_str, "%RA", x);

  /* regular numbers */
  mpfr_set_str (x, "FEDCBA9.87654321", 16, GMP_RNDN);
  mpfr_set_ui (z, 0, GMP_RNDZ);

  /* simplest case right justified */
  check_sprintf ("   0xf.edcba987654321p+24", "%25Ra", x);
  check_sprintf ("                  0x8p+25", "%25.0Ra", x);
  check_sprintf ("                   0x0p+0", "%25.0Ra", z);
  /* sign or space, pad with leading zeros */
  check_sprintf (" 0X00F.EDCBA987654321P+24", "% 025RA", x);
  check_sprintf (" 0X000000000000000008P+25", "% 025.0RA", x);
  check_sprintf (" 0X0000000000000000000P+0", "% 025.0RA", z);
  /* sign + or -, left justified */
  check_sprintf ("+0xf.edcba987654321p+24  ", "%+-25Ra", x);
  check_sprintf ("+0x8p+25                 ", "%+-25.0Ra", x);
  check_sprintf ("+0x0p+0                  ", "%+-25.0Ra", z);
  /* decimal point, left justified, precision and rounding parameter */
  check_vsprintf ("0XF.FP+24 ", "%#-10.*R*A", 1, GMP_RNDN, x);
  check_vsprintf ("0X8.P+25  ", "%#-10.*R*A", 0, GMP_RNDN, x);
  check_vsprintf ("0X0.P+0   ", "%#-10.*R*A", 0, GMP_RNDN, z);
  /* sign or space */
  check_sprintf (" 0xf.eddp+24", "% .3RNa", x);
  check_sprintf (" 0x8p+25",     "% .0RNa", x);
  /* sign + or -, decimal point, pad with leading zeros */
  check_sprintf ("+0X0F.EP+24", "%0+#11.1RZA", x);
  check_sprintf ("+0X00F.P+24", "%0+#11.0RZA", x);
  check_sprintf ("+0X000.0P+0", "%0+#11.1RZA", z);
  /* pad with leading zero */
  check_sprintf ("0x0000f.edcba987654321p+24", "%026RDa", x);
  check_sprintf ("0x0000000000000000000fp+24", "%026.0RDa", x);
  /* sign or space, decimal point, left justified */
  check_sprintf (" 0XF.EP+24 " , "%- #11.1RDA", x);
  check_sprintf (" 0XF.P+24  " , "%- #11.0RDA", x);

  mpfr_mul_si (x, x, -1, GMP_RNDD);
  mpfr_mul_si (z, z, -1, GMP_RNDD);

  /* sign + or - */
  check_sprintf ("-0xf.ep+24", "%+10.1RUa", x);
  check_sprintf ("  -0xfp+24", "%+10.0RUa", x);
  check_sprintf ("   -0x0p+0", "%+10.0RUa", z);

  mpfr_clears (x, z, (void *)0);
  return 0;
}

static int
binary (void)
{
  mpfr_t x;
  mpfr_t z;
  mpfr_inits2 (64, x, z, (void *)0);

  /* special */
  mpfr_set_inf (x, 1);
  check_sprintf (pinf_str, "%Rb", x);

  mpfr_set_inf (x, -1);
  check_sprintf (minf_str, "%Rb", x);

  mpfr_set_nan (x);
  check_sprintf (nan_str, "%Rb", x);

  /* regular numbers */
  mpfr_set_str (x, "1110010101.1001101", 2, GMP_RNDN);
  mpfr_set_ui (z, 0, GMP_RNDN);

  /* simplest case: right justified */
  check_sprintf ("    1.1100101011001101p+9", "%25Rb", x);
  check_sprintf ("                     0p+0", "%25Rb", z);
  /* sign or space, pad with leading zeros */
  check_sprintf (" 0001.1100101011001101p+9", "% 025Rb", x);
  check_sprintf (" 000000000000000000000p+0", "% 025Rb", z);
  /* sign + or -, left justified */
  check_sprintf ("+1.1100101011001101p+9   ", "%+-25Rb", x);
  check_sprintf ("+0p+0                    ", "%+-25Rb", z);
  /* sign or space */
  check_sprintf (" 1.110p+9",  "% .3RNb", x);
  check_sprintf (" 1.1101p+9", "% .4RNb", x);
  check_sprintf (" 0.0000p+0", "% .4RNb", z);
  /* sign + or -, decimal point, pad with leading zeros */
  check_sprintf ("+00001.1p+9", "%0+#11.1RZb", x);
  check_sprintf ("+0001.0p+10", "%0+#11.1RNb", x);
  check_sprintf ("+000000.p+0", "%0+#11.0RNb", z);
  /* pad with leading zero */
  check_sprintf ("00001.1100101011001101p+9", "%025RDb", x);
  /* sign or space, decimal point (unused), left justified */
  check_sprintf (" 1.1p+9    " , "%- #11.1RDb", x);
  check_sprintf (" 1.1p+9    " , "%- #11.0RDb", x);

  mpfr_mul_si (x, x, -1, GMP_RNDD);
  mpfr_mul_si (z, z, -1, GMP_RNDD);

  /* sign + or - */
  check_sprintf ("   -1.1p+9", "%+10.1RUb", x);
  check_sprintf ("   -0.0p+0", "%+10.1RUb", z);

  mpfr_clears (x, z, (void *)0);
  return 0;
}

static int
mixed (void)
{
  int n1;
  int n2;
  int i = 121;
  long double d = 1. / 31.;
  mpf_t mpf;
  mpq_t mpq;
  mpz_t mpz;
  mpfr_t x;
  mp_rnd_t rnd;

  mpf_init (mpf);
  mpf_set_ui (mpf, 40);
  mpf_div_ui (mpf, mpf, 31); /* mpf = 40.0 / 31.0 */
  mpq_init (mpq);
  mpq_set_ui (mpq, 123456, 4567890);
  mpz_init (mpz);
  mpz_fib_ui (mpz, 64);
  mpfr_init (x);
  mpfr_set_si (x, -875, GMP_RNDN);
  mpfr_div_ui (x, x, 1000, GMP_RNDN);
  mpfr_add_si (x, x, -12345678, GMP_RNDN); /* x = -1.2345678875e7 */
  rnd = GMP_RNDD;

  check_vsprintf ("121", "%i", i);
  check_vsprintf ("121, -12345679", "%i, %.0Rf", i, x);
  check_vsprintf ("10610209857723, -1.2345678875e+07", "%Zi, %R*e", mpz, rnd,
                  x);
  check_vsprintf ("-12345678.9, 121", "%.1Rf, %i", x, i);
  check_vsprintf ("-12345678, 1e240/45b352", "%.0R*f, %Qx", GMP_RNDZ, x, mpq);
  check_vsprintf ("121, -12345678.875000000000, 1.290323", "%i, %.*Rf, %Ff",
                  i, 12, x, mpf);
  n1 = check_vsprintf ("00000010610209857723, -1.2345678875e+07, 0.032258",
                       "%.*Zi, %R*e, %Lf%n", 20, mpz, rnd, x, d, &n2);

  if (n1 != n2)
    {
      printf ("error in number of characters written by mpfr_vsprintf\n");
      printf ("expected: %d\n", n2);
      printf ("     got: %d\n", n1);
      exit (1);
    }
  mpf_clear (mpf);
  mpq_clear (mpq);
  mpz_clear (mpz);
  mpfr_clear (x);
  return 0;
}

static int
random_double (void)
{
  int i;
  mpfr_t x;
  double y;
  char flag[] =
    {
      '-',
      '+',
      ' ',
      '#',
      '0' /* no ambiguity: first zeros are flag zero*/
    };
  /* no 'a': mpfr and glibc do not have the same semantic */
  char specifier[] =
    {
      'e',
      'f',
      'g',
      'E',
      'f', /* SUSv2 doesn't accept %F, but %F and %f are the same for
              regular numbers */
      'G',
    };

  mpfr_init2 (x, 53);

  for (i=0; i<1000; )
    {
      y = DBL_RAND ();
      if (!Isnan(y))
        {
          int j, spec, prec;
          char fmt_mpfr[11];
          char *ptr_mpfr = fmt_mpfr;
          char fmt[10];
          char *ptr = fmt;
          int xi;
          char *xs;
          int yi;
          char *ys;

          i++;

          if ((rand() / (RAND_MAX + 1.0)) > .5)
            y = -y;
          mpfr_set_d (x, y, GMP_RNDN);

          *ptr_mpfr++ = *ptr++ = '%';
          for (j = 0; j < 5; j++)
            {
              if ((rand() / (RAND_MAX + 1.0)) < .3)
                *ptr_mpfr++ = *ptr++ = flag[j];
            }
          *ptr_mpfr++ = *ptr++ = '.';
          *ptr_mpfr++ = *ptr++ = '*';
          *ptr_mpfr++ = 'R';
          spec = (int) (6.0 * (rand() / (RAND_MAX + 1.0)));
          *ptr_mpfr++ = *ptr++ = specifier[spec];
          *ptr_mpfr = *ptr = '\0';

          /* advantage small precision */
          if ((rand() / (RAND_MAX + 1.0)) > .5)
            prec = (int) (10 * (rand() / (RAND_MAX + 1.0)));
          else
            prec = (int) (prec_max_printf * (rand() / (RAND_MAX + 1.0)));

          if (y != mpfr_get_d (x, GMP_RNDN))
            /* conversion error: skip this one */
            continue;

          xi = mpfr_asprintf (&xs, fmt_mpfr, prec, x);
          yi = mpfr_asprintf (&ys, fmt, prec, y);

          if (xi != yi || strcmp (xs, ys))
            {
              mpfr_printf ("Error in mpfr_asprintf(\"%%%s\", %d, %Re)\n" \
                      "expected: ", fmt_mpfr, prec, x);
              printf ("%s", ys);
              printf ("\n     got: ");
              printf ("%s", xs);
              putchar ('\n');

              exit (1);
            }

          mpfr_free_str (xs);
          mpfr_free_str (ys);
        }
    }

  mpfr_clear (x);
  return 0;
}

int
main (int argc, char **argv)
{
  char *locale;

  tests_start_mpfr ();

#if defined(HAVE_LOCALE_H) && defined(HAVE_SETLOCALE)
  /* currently, we just check with 'C' locale */
  locale = setlocale (LC_NUMERIC, "C");
#endif

  hexadecimal ();
  binary ();
  decimal ();
  mixed ();
  random_double ();

#if defined(HAVE_LOCALE_H) && defined(HAVE_SETLOCALE)
  setlocale (LC_NUMERIC, locale);
#endif

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
