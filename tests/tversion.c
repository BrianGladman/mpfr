/* Test file for mpfr_version.

Copyright 2004-2015 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

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

#include "mpfr-test.h"

int
main (void)
{
  int err = 0;

  /* Test the GMP and MPFR versions. */
  if (test_version ())
    exit (1);

  printf ("[tversion] MPFR %s\n", MPFR_VERSION_STRING);

  /* TODO: We may want to output info for non-GNUC-compat compilers too. See:
   * http://sourceforge.net/p/predef/wiki/Compilers/
   * http://nadeausoftware.com/articles/2012/10/c_c_tip_how_detect_compiler_name_and_version_using_compiler_predefined_macros
   */
#if defined(__GNUC__) && defined(__VERSION__)
  printf ("[tversion] GCC: %s\n", __VERSION__);
#endif

#ifdef __MPIR_VERSION
  printf ("[tversion] MPIR: header %d.%d.%d, library %s\n",
          __MPIR_VERSION, __MPIR_VERSION_MINOR, __MPIR_VERSION_PATCHLEVEL,
          mpir_version);
#else
#ifdef MPFR_USE_MINI_GMP
  printf ("[tversion] mini-gmp\n");
#else
  printf ("[tversion] GMP: header %d.%d.%d, library %s\n",
          __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL,
          gmp_version);
#endif
#endif

  if (
#ifdef MPFR_USE_THREAD_SAFE
      !
#endif
      mpfr_buildopt_tls_p ())
    {
      printf ("ERROR! mpfr_buildopt_tls_p() and macros"
              " do not match!\n");
      err = 1;
    }

  if (
#ifdef MPFR_WANT_DECIMAL_FLOATS
      !
#endif
      mpfr_buildopt_decimal_p ())
    {
      printf ("ERROR! mpfr_buildopt_decimal_p() and macros"
              " do not match!\n");
      err = 1;
    }

  if (
#if defined(MPFR_HAVE_GMP_IMPL) || defined(WANT_GMP_INTERNALS)
      !
#endif
      mpfr_buildopt_gmpinternals_p ())
    {
      printf ("ERROR! mpfr_buildopt_gmpinternals_p() and macros"
              " do not match!\n");
      err = 1;
    }

  printf ("[tversion] TLS = %s, decimal = %s, GMP internals = %s\n",
          mpfr_buildopt_tls_p () ? "yes" : "no",
          mpfr_buildopt_decimal_p () ? "yes" : "no",
          mpfr_buildopt_gmpinternals_p () ? "yes" : "no");

  if (strcmp (mpfr_buildopt_tune_case (), MPFR_TUNE_CASE) != 0)
    {
      printf ("ERROR! mpfr_buildopt_tune_case() and MPFR_TUNE_CASE"
              " do not match!\n  %s\n  %s\n",
              mpfr_buildopt_tune_case (), MPFR_TUNE_CASE);
      err = 1;
    }
  else
    printf ("[tversion] MPFR tuning parameters from %s\n", MPFR_TUNE_CASE);

  if (strcmp (mpfr_get_patches (), "") != 0)
    printf ("[tversion] MPFR patches: %s\n", mpfr_get_patches ());

  tests_start_mpfr ();
  if (locale != NULL)
    printf ("[tversion] Locale: %s\n", locale);
  tests_end_mpfr ();

  return err;
}
