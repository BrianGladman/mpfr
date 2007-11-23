/* mpfr_tprintf -- test file for mpfr_sprintf

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_LOCALE_H
#include <locale.h>
#endif

#include <gmp.h>

#include "mpfr-test.h"

const int buf_size = 1024;
__gmp_const char pinf_str[] = "inf";
__gmp_const char pinf_uc_str[] = "INF";
__gmp_const char minf_str[] = "-inf";
__gmp_const char minf_uc_str[] = "-INF";
__gmp_const char nan_str[] = "nan";
__gmp_const char nan_uc_str[] = "NAN";

static void
special ()
{
  char buffer[buf_size];
  mpfr_t x;
  mpfr_init (x);

  mpfr_set_inf (x, 1);
  mpfr_sprintf (buffer, "%Rf", x);
  if (strcmp (buffer, pinf_str) != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf\nexpected: %s\ngot:      %s\n",
               pinf_str, buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%RF", x);
  if (strcmp (buffer, pinf_uc_str) != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf\nexpected: %s\ngot:      %s\n",
               pinf_uc_str, buffer);
      exit (1);
    }

  mpfr_set_inf (x, -1);
  mpfr_sprintf (buffer, "%Rf", x);
  if (strcmp (buffer, minf_str) != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf\nexpected: %s\ngot:      %s\n",
               minf_str, buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%RF", x);
  if (strcmp (buffer, minf_uc_str) != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf\nexpected: %s\ngot:      %s\n",
               minf_uc_str, buffer);
      exit (1);
    }

  mpfr_set_nan (x);
  mpfr_sprintf (buffer, "%Rf", x);
  if (strcmp (buffer, nan_str) != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf\nexpected: %s\ngot:      %s\n",
               nan_str, buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%RF", x);
  if (strcmp (buffer, nan_uc_str) != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf\nexpected: %s\ngot:      %s\n",
               nan_uc_str, buffer);
      exit (1);
    }

  mpfr_clear (x);
}

static int
integer ()
{
  char buffer[buf_size];
  mpfr_t x;
  mpfr_init (x);

  mpfr_set_d (x, 1895485593474.61279296875, GMP_RNDD);
  /* base ten */
  mpfr_sprintf (buffer, "%RDd", x);
  if (strcmp (buffer, "1895485593474") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%RDd\", x)\n");
      fprintf (stderr, "expected: 1895485593474\ngot:      %s\n", buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%RNi", x);
  if (strcmp (buffer, "1895485593475") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%RNi\", x)\n");
      fprintf (stderr, "expected: 1895485593475\ngot:      %s\n", buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%RUu", x);
  if (strcmp (buffer, "1895485593475") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%RUu\", x)\n");
      fprintf (stderr, "expected: 1895485593475\ngot:      %s\n", buffer);
      exit (1);
    }
  /* base sixteen */
  mpfr_sprintf (buffer, "%RZx", x);
  if (strcmp (buffer, "1b953bed782") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%RZx\", x)\n");
      fprintf (stderr, "expected: 1b953bed782\ngot:      %s\n", buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%#RNX", x);
  if (strcmp (buffer, "0X1B953BED783") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%#RNX\", x)\n");
      fprintf (stderr, "expected: 0X1B953BED783\ngot:      %s\n", buffer);
      exit (1);
    }
  /* base eight */
  mpfr_sprintf (buffer, "%RNo", x);
  if (strcmp (buffer, "33452357553603") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%RNo\", x)\n");
      fprintf (stderr, "expected: 33452357553603\ngot:      %s\n", buffer);
      exit (1);
    }

  /* flags, width, and precision checking */
  mpfr_set_si (x, -1641, GMP_RNDD);
  mpfr_sprintf (buffer, "%+012RDd", x);
  if (strcmp (buffer, "-00000001641") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%+012RDd\", x)\n");
      fprintf (stderr, "expected: \"-00000001641\"\ngot:      %s\n", buffer);
      exit (1);
    }
  mpfr_neg (x, x, GMP_RNDD);
  mpfr_sprintf (buffer, "%-*.*RDd", 12, 6, x);
  if (strcmp (buffer, "001641      ") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%-*.*RDd\", 12, 6, x)\n");
      fprintf (stderr, "expected: \"001641      \"\ngot:      \"%s\"\n",
	       buffer);
      exit (1);
    }

  mpfr_clear (x);
  return 0;
}

static int
floating_point ()
{
  char buffer[buf_size];
  mpfr_t x;
  mpfr_init (x);
  mpfr_set_d (x, 1895485593474.61279296875, GMP_RNDD);

  /* default: 6 decimal digits */
  mpfr_sprintf (buffer, "%RDf", x);
  if (strcmp (buffer, "1895485593474.612792") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%RDf\", x)\n");
      fprintf (stderr, "expected: 1895485593474.612792\ngot:      %s\n",
	       buffer);
      exit (1);
    }
  /* test rounding */
  mpfr_sprintf (buffer, "%.4RNf", x);
  if (strcmp (buffer, "1895485593474.6128") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%.4RNf\", x)\n");
      fprintf (stderr, "expected: 1895485593475.6128\ngot:      %s\n", buffer);
      exit (1);
    }
  /* decimal point, no decimal digit */
  mpfr_sprintf (buffer, "%#.0RUf", x);
  if (strcmp (buffer, "1895485593475.") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%#.0RUf\", x)\n");
      fprintf (stderr, "expected: 1895485593475.\ngot:      %s\n", buffer);
      exit (1);
    }
  /* request the significant digits */
  mpfr_sprintf (buffer, "%.RDe", x);
  if (strcmp (buffer, "1.8954855934746127e+12") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%.RDe\", x)\n");
      fprintf (stderr,
	       "expected: 1.8954855934746127e+12\ngot:      %s\n",
	       buffer);
      exit (1);
    }
  /* test width field */
  mpfr_sprintf (buffer, "%10.3RNE", x);
  if (strcmp (buffer, " 1.895E+12") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%10.3RNE\", x)\n");
      fprintf (stderr, "expected: \" 1.895E+12\"\ngot:      \"%s\"\n", buffer);
      exit (1);
    }
  /* show sign */
  mpfr_sprintf (buffer, "%+10.7RUe", x);
  if (strcmp (buffer, "+1.8954856e+12") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%+10.7RUe\", x)\n");
      fprintf (stderr, "expected: +1.8954856e+12\ngot:      %s\n", buffer);
      exit (1);
    }
  /* 'g' conversion specifier */
  mpfr_sprintf (buffer, "=%.3RNG", x);
  if (strcmp (buffer, "=1.90E+12") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"=%%.3RNG\", x)\n");
      fprintf (stderr, "expected: \"=1.90E+12\"\ngot:      \"%s\"\n", buffer);
      exit (1);
    }
  mpfr_set_d (x, 11.875, GMP_RNDD);
  mpfr_sprintf (buffer, "%RUg", x);
  if (strcmp (buffer, "11.88") != 0) /* [FIXME] a verifier */
    {
      fprintf (stderr, "Error in mpfr_sprintf(s, \"%%RUg\", x)\n");
      fprintf (stderr, "expected: 11.88\ngot:      %s\n", buffer);
      exit (1);
    }


  mpfr_clear (x);
  return 0;
}

static int
mixed ()
{
  int i = 121;
  long double d = 1. / 31.;
  mpf_t mpf;
  mpq_t mpq;
  mpz_t mpz;
  char buffer[buf_size];
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

  mpfr_sprintf (buffer, "%i", i, x);
  if (strcmp (buffer, "121") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf (s, \"%%i\", i, x);\n");
      fprintf (stderr, "expected: 121\ngot:      %s\n", buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%i, %.0Rf", i, x);
  if (strcmp (buffer, "121, -12345679") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf (s, \"%%i, %%.0Rf\", i, x);\n");
      fprintf (stderr, "expected: 121, -12345679\ngot:      %s\n",
	       buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%Zi, %R*e", mpz, rnd, x);
  if (strcmp (buffer, "10610209857723, -1.2345678875e+07") != 0)
    {
      fprintf (stderr,
	       "Error in mpfr_sprintf (s, \"%%Zi, %%R*e\", mpz, rnd, x);\n");
      fprintf (stderr,
	       "expected: 10610209857723, -1.2345678875e+07\ngot:      %s\n",
	       buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%.1Rf, %i", x, i);
  if (strcmp (buffer, "-12345678.9, 121") != 0)
    {
      fprintf (stderr, "Error in mpfr_sprintf (s, \"%%.1Rf, %%i\", x, i);\n");
      fprintf (stderr, "expected: -12345678.9, 121\ngot:      %s\n", buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%.0R*f, %Qx", GMP_RNDZ, x, mpq);
  if (strcmp (buffer, "-12345678, 1e240/45b352") != 0)
    {
      fprintf (stderr,
	       "Error in mpfr_sprintf (s, \"%%R*e, %%Qx\", GMP_RNDZ, x, mpq)\n");
      fprintf (stderr, "expected: -12345678, 1e240/45b352\ngot:      %s\n",
	       buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%i, %.*Rf, %Ff", i, 12, x, mpf);
  if (strcmp (buffer, "121, -12345678.875000000000, 1.290323") != 0)
    {
      fprintf (stderr,
	       "Error in mpfr_sprintf (s, \"%%i, %%.*Rf, %%Ff\", i, 12, x, mpf)\n");
      fprintf (stderr, \
	       "expected: 121, -12345678.875000000000, 1.290323\ngot:      %s\n",
	       buffer);
      exit (1);
    }
  mpfr_sprintf (buffer, "%.*Zi, %R*e, %Lf", 20, mpz, rnd, x, d);
  if (strcmp (buffer, "00000010610209857723, -1.2345678875e+07, 0.032258") != 0)
    {
      fprintf (stderr,
	       "Error in mpfr_sprintf (s,\"%%.*Zi, %%R*e, %%Lf\", 20, mpz, GMP_RNDD, x, d)\n");
      fprintf (stderr, \
	       "expected: 00000010610209857723, -1.2345678875e+07, 0.032258\ngot:      %s\n",
	       buffer);
      exit (1);
    }

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

#ifdef HAVE_LOCALE_H && HAVE_SETLOCALE
  /* currently, we just check with 'C' locale */
  locale = setlocale (LC_NUMERIC, "C");
#endif

  special ();
  integer ();
  floating_point ();
  mixed ();

#ifdef HAVE_LOCALE_H && HAVE_SETLOCALE
  setlocale (LC_NUMERIC, locale);
#endif

  tests_end_mpfr ();
  return 0;
}
