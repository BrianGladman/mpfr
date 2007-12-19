/* mpfr_tprintf -- test file for mpfr_sprintf mpfr_vsprintf

Copyright 2007 Free Software Foundation, Inc.
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

const int buf_size = 1024;
const char pinf_str[] = "inf";
const char pinf_uc_str[] = "INF";
const char minf_str[] = "-inf";
const char minf_uc_str[] = "-INF";
const char nan_str[] = "nan";
const char nan_uc_str[] = "NAN";

/* compare expected string with the return of mpfr_sprintf(fmt, x)*/
static void
check_sprintf (const char *expected, const char *fmt, mpfr_srcptr x)
{
  char buffer[buf_size];

  mpfr_sprintf (buffer, fmt, x);
  if (strcmp (buffer, expected) != 0)
    {
      printf ("Error in mpfr_sprintf (s, \"%s\", x);\n", fmt);
      printf ("expected: \"%s\"\ngot:      \"%s\"\n", expected, buffer);

      exit (1);
    }
}

/* compare expected string with the return of mpfr_vsprintf(fmt, ...)*/
static void
check_vsprintf (const char *expected, const char *fmt, ...)
{
  char buffer[buf_size];
  va_list ap;
  va_start (ap, fmt);

  mpfr_vsprintf (buffer, fmt, ap);
  if (strcmp (buffer, expected) != 0)
    {
      printf ("Error in mpfr_vsprintf (s, \"%s\", ...);\n", fmt);
      printf ("expected: \"%s\"\ngot:      \"%s\"\n", expected, buffer);

      va_end (ap);
      exit (1);
    }

  va_end (ap);
}

static int
integer (void)
{
  mpfr_t x;
  mpfr_init (x);

  /* special values */
  mpfr_set_inf (x, 1);
  check_sprintf (pinf_str, "%Rd", x);
  check_sprintf (pinf_uc_str, "%RX", x);

  mpfr_set_inf (x, -1);
  check_sprintf (minf_str, "%Ro", x);
  check_sprintf (minf_uc_str, "%RX", x);

  mpfr_set_nan (x);
  check_sprintf (nan_str, "%Rd", x);
  check_sprintf (nan_uc_str, "%RX", x);

  /* regular numbers */
  mpfr_set_d (x, 1895485593474.61279296875, GMP_RNDD);

  /* base ten */
  check_sprintf ("1895485593474", "%RDd", x);
  check_sprintf ("1895485593475", "%RNi", x);
  check_sprintf ("1895485593475", "%RUu", x);

  /* base sixteen */
  check_sprintf ("1b953bed782", "%RZx", x);
  check_sprintf ("0X1B953BED783", "%#RNX", x);

  /* base eight */
  check_sprintf ("33452357553603", "%RNo", x);

  /* flags, width, and precision checking */
  mpfr_set_si (x, -1641, GMP_RNDD);
  check_sprintf ("-00000001641", "%+012RDd", x);

  mpfr_neg (x, x, GMP_RNDD);
  check_vsprintf ("001641      ", "%-*.*RDd", 12, 6, x);

  mpfr_clear (x);
  return 0;
}

static int
floating_point (void)
{
  mpfr_t x;
  mpfr_init (x);

  /* special numbers */
  mpfr_set_inf (x, 1);
  check_sprintf (pinf_str, "%Rf", x);
  check_sprintf (pinf_uc_str, "%RF", x);

  mpfr_set_inf (x, -1);
  check_sprintf (minf_str, "%Rf", x);
  check_sprintf (minf_uc_str, "%RF", x);

  mpfr_set_nan (x);
  check_sprintf (nan_str, "%Rf", x);
  check_sprintf (nan_uc_str, "%RF", x);

  /* regular numbers */
  mpfr_set_d (x, 1895485593474.61279296875, GMP_RNDD);

  /* default: 6 decimal digits */
  check_sprintf ("1895485593474.612792", "%RDf", x);
  /* test rounding */
  check_sprintf ("1895485593474.6128", "%.4RNf", x);
  /* decimal point, no decimal digit */
  check_sprintf ("1895485593475.", "%#.0RUf", x);
  /* request the significant digits */
  check_sprintf ("1.8954855934746127e+12", "%.RDe", x);
  /* test width field */
  check_sprintf (" 1.895E+12", "%10.3RNE", x);
  /* show sign */
  check_sprintf ("+1.8954856e+12", "%+10.7RUe", x);
  /* 'g' conversion specifier */
  check_sprintf ("=1.90E+12", "=%.3RNG", x);

  /* this case need to be fix in mpfr_vasprintf */
/*   mpfr_set_d (x, 11.875, GMP_RNDD); */
/*   check_sprintf ("11.875", "%RUg", x); */

  mpfr_clear (x);
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
  int i = 121;
  long double d = 1. / 31.;
  mpf_t mpf;
  mpq_t mpq;
  mpz_t mpz;
  mpfr_t x;
  mp_rnd_t rnd;

  mpf_init (mpf);
  mpf_set_d (mpf, 40.0 * d);
  mpq_init (mpq);
  mpq_set_ui (mpq, 123456, 4567890);
  mpz_init (mpz);
  mpz_fib_ui (mpz, 64);
  mpfr_init (x);
  mpfr_set_d (x, -1.2345678875e7, GMP_RNDN);
  rnd = GMP_RNDD;

  check_vsprintf ("121", "%i", i);
  check_vsprintf ("121, -12345679", "%i, %.0Rf", i, x);
  check_vsprintf ("10610209857723, -1.2345678875e+07", "%Zi, %R*e", mpz, rnd,
                  x);
  check_vsprintf ("-12345678.9, 121", "%.1Rf, %i", x, i);
  check_vsprintf ("-12345678, 1e240/45b352", "%.0R*f, %Qx", GMP_RNDZ, x, mpq);
  check_vsprintf ("121, -12345678.875000000000, 1.290323", "%i, %.*Rf, %Ff",
                  i, 12, x, mpf);
  check_vsprintf ("00000010610209857723, -1.2345678875e+07, 0.032258",
                  "%.*Zi, %R*e, %Lf", 20, mpz, rnd, x, d);

  mpf_clear (mpf);
  mpq_clear (mpq);
  mpz_clear (mpz);
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

  integer ();
  hexadecimal ();
  binary ();
  floating_point ();  /* [TODO] */
  mixed ();           /* [TODO] */

#if defined(HAVE_LOCALE_H) && defined(HAVE_SETLOCALE)
  setlocale (LC_NUMERIC, locale);
#endif

  tests_end_mpfr ();
  return 0;
}

#else  /* HAVE_STDARG */

int
main (int argc, char **argv)
{
  /* We have nothing to test. */
  return 0;
}

#endif  /* HAVE_STDARG */
