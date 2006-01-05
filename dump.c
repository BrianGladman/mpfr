/* mpfr_dump -- Dump a float to stdout.

Copyright 1999, 2001, 2004, 2006 Free Software Foundation, Inc.
(Copied from the GNU MP Library.)

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
#include <string.h>

#include "mpfr-impl.h"

void
mpfr_dump (mpfr_srcptr u)
{
  mpfr_print_binary(u);
  putchar('\n');
}
