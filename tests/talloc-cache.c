/* talloc-cache  -- test file concerning memory allocation and cache

Copyright 2016-2017 Free Software Foundation, Inc.
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

/* WARNING! This is not the final code and this code will be incorrect
   once some issues with memory allocation are resolved. It is currently
   for temporary testing only.
   In the 3.1 branch, it crashes before r9467, and r9467 fixes the crash
   (among other issues).
   In the trunk, it currently fails due to the mpz_t pool (this will be
   fixed by the change of memory allocation method).
*/

#include <stdlib.h>

#include "mpfr-test.h"

#define A 4096
#define I 10000

static void *
my_alloc1 (size_t s)
{
  void *p = malloc (s + A);
  *(int *) p = 1;
  return (void *) ((char *) p + A);
}

static void *
my_realloc1 (void *p, size_t t, size_t s)
{
  p = (void *) ((char *) p - A);
  if (*(int *) p != 1)
    abort ();
  return (void *) ((char *) realloc (p, s + A) + A);
}

static void
my_free1 (void *p, size_t t)
{
  p = (void *) ((char *) p - A);
  if (*(int *) p != 1)
    abort ();
  free (p);
}

static void *
my_alloc2 (size_t s)
{
  void *p = malloc (s + A);
  *(int *) p = 2;
  return (void *) ((char *) p + A);
}

static void *
my_realloc2 (void *p, size_t t, size_t s)
{
  p = (void *) ((char *) p - A);
  if (*(int *) p != 2)
    abort ();
  return (void *) ((char *) realloc (p, s + A) + A);
}

static void
my_free2 (void *p, size_t t)
{
  p = (void *) ((char *) p - A);
  if (*(int *) p != 2)
    abort ();
  free (p);
}

int
main (void)
{
  mpfr_t x;

  tests_memory_disabled = 2;
  tests_start_mpfr ();

  mp_set_memory_functions (my_alloc1, my_realloc1, my_free1);

  mpfr_init2 (x, 53);
  mpfr_set_ui (x, I, MPFR_RNDN);
  mpfr_sin (x, x, MPFR_RNDN);
  mpfr_clear (x);

  mp_set_memory_functions (my_alloc2, my_realloc2, my_free2);

  mpfr_init2 (x, 1000);
  mpfr_set_ui (x, I, MPFR_RNDN);
  mpfr_sin (x, x, MPFR_RNDN);
  mpfr_clear (x);

  tests_end_mpfr ();
  return 0;
}
