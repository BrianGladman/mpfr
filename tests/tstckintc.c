/* Test file for mpfr_stack_*

Copyright 2005 Free Software Foundation.

This file is part of the MPFR Library.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "mpfr-test.h"

#define BUFFER_SIZE 1000
#define PREC_TESTED 200

char Buffer[BUFFER_SIZE];
char *stack = Buffer;
mp_prec_t p = PREC_TESTED;

#define ALIGNED(s) ( ((s)+sizeof (long)-1) / sizeof (long) * sizeof (long))

static void *
new (size_t s)
{
  void *p = (void*) stack;
  stack += ALIGNED (s);
  if (MPFR_UNLIKELY (stack > &Buffer[BUFFER_SIZE]))
    {
      printf ("Stack overflow.\n");
      exit (1);
    }
  return p;
}

 /* Alloc a new mpfr_t on the main stack */
static mpfr_ptr
new_mpfr (mp_prec_t p)
{
  mpfr_ptr x = (mpfr_ptr) new (sizeof (mpfr_t));
  void *mantissa = new (mpfr_stack_get_size (p));
  mpfr_stack_init (mantissa, p);
  mpfr_stack_init_set (x, 0, 0, p, mantissa);
  return x;
}

 /* Garbage the stack by keeping only x  */
static mpfr_ptr
return_mpfr (mpfr_ptr x, void *old_stack)
{
  void *mantissa       = mpfr_stack_get_mantissa (x);
  size_t size_mantissa = mpfr_stack_get_size (mpfr_get_prec (x));
  mpfr_ptr newx;

  memmove (old_stack, x, sizeof (mpfr_t));
  memmove (old_stack + ALIGNED (sizeof (mpfr_t)), mantissa, size_mantissa);
  newx = (mpfr_ptr) old_stack;
  mpfr_stack_move (newx, old_stack + ALIGNED (sizeof (mpfr_t)));
  stack = old_stack + ALIGNED (sizeof (mpfr_t)) + ALIGNED (size_mantissa);
  return newx;
}

static void
test1 (void)
{
  mpfr_ptr x, y;
  void *org;

  org = stack;
  x = new_mpfr (p);
  y = new_mpfr (p);
  mpfr_set_ui (x, 42, GMP_RNDN);
  mpfr_set_ui (y, 17, GMP_RNDN);
  mpfr_add (y, x, y, GMP_RNDN);
  y = return_mpfr (y, org);
  if (y != x || mpfr_cmp_ui (y, 59) != 0)
    {
      printf ("Compact (1) failed!\n");
      exit (1);
    }
  stack = org;
}

/* We build the MPFR variable each time it is needed */
/* a[0] is the kind, a[1] is the exponent, &a[2] is the mantissa */
static long *
dummy_new (void)
{
  long *r;

  r = new (ALIGNED (2*sizeof (long)) + ALIGNED (mpfr_stack_get_size (p)));
  MPFR_ASSERTN (r != NULL);
  (mpfr_stack_init) (&r[2], p);
  r[0] = (int) MPFR_NAN_KIND;
  r[1] = 0;
  return r;
}

static long *
dummy_set_si (long si)
{
  mpfr_t x;
  long * r = dummy_new ();
  (mpfr_stack_init_set) (x, 0, 0, p, &r[2]);
  mpfr_set_si (x, si, GMP_RNDN);
  r[0] = mpfr_stack_get_kind (x);
  r[1] = mpfr_stack_get_exp (x);
  return r;
}

static long *
dummy_add (long *a, long *b)
{
  mpfr_t x, y, z;
  long *r = dummy_new ();
  mpfr_stack_init_set (x, 0, 0, p, &r[2]);
  (mpfr_stack_init_set) (y, a[0], a[1], p, &a[2]);
  mpfr_stack_init_set (z, b[0], b[1], p, &b[2]);
  mpfr_add (x, y, z, GMP_RNDN);
  r[0] = (mpfr_stack_get_kind) (x);
  r[1] = (mpfr_stack_get_exp) (x);
  return r;
}

static long *
dummy_compact (long *r, void *org_stack)
{
  memmove (org_stack, r,
           ALIGNED (2*sizeof (long)) + ALIGNED ((mpfr_stack_get_size) (p)));
  return org_stack;
}

static void
test2 (void)
{
  mpfr_t x;
  void *org = stack;
  long *a, *b, *c;

  a = dummy_set_si (42);
  b = dummy_set_si (17);
  c = dummy_add (a, b);
  c = dummy_compact (c, org);
  (mpfr_stack_init_set) (x, c[0], c[1], p, &c[2]);
  if (c != a || mpfr_cmp_ui (x, 59) != 0)
    {
      printf ("Compact (2) failed! c=%p a=%p\n", c, a);
      mpfr_dump (x);
      exit (1);
    }
  stack = org;
}

int
main (void)
{
  tests_start_mpfr ();
  /* We test iff long = mp_limb_t */
  if (sizeof (long) == sizeof (mp_limb_t))
    {
      test1 ();
      test2 ();
    }
  tests_end_mpfr ();
  return 0;
}
