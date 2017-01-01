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

/* This test program is for temporary testing only. It is not in the
   usual test suite (not in check_PROGRAMS from Makefile.am).
   In the 3.1 branch, it crashes before r9467, and r9467 fixes the crash
   (among other issues).
   In the trunk, it crashes, which is not expected with the current code.
*/

#include <stdlib.h>

#include "mpfr-test.h"

#define A 4096
#define I 10000

static void *
my_alloc1 (size_t s)
{
return malloc (s);
}

static void *
my_realloc1 (void *p, size_t t, size_t s)
{
return realloc (p, s);
}

static void
my_free1 (void *p, size_t t)
{
  return free (p);
}

static void *
my_alloc2 (size_t s)
{
  return (void *) ((char *) malloc (s + A) + A);
}

static void *
my_realloc2 (void *p, size_t t, size_t s)
{
  return (void *) ((char *) realloc ((char *) p - A, s + A) + A);
}

static void
my_free2 (void *p, size_t t)
{
  return free ((char *) p - A);
}

int
main (void)
{
  mpfr_t x;

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

  return 0;
}
