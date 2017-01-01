/* mpfr_dump -- Dump a float to stdout.

Copyright 1999, 2001, 2004, 2006-2017 Free Software Foundation, Inc.
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

#include "mpfr-impl.h"

/* mpfr_dump is mainly for debugging purpose. It outputs a MPFR number
 * in some unspecified format, then a newline character. This function
 * is part of the API, but the output format may change without breaking
 * the ABI.
 *
 * TODO: Since this is for debugging, it should support invalid data,
 * such as:
 *   - out-of-range exponents (with some form of warning);
 *   - numbers that are not normalized (MSB = 0);
 *   - numbers with trailing bits that are not all zero's
 *     (e.g. output in square brackets in such a case).
 * The format may depend on MPFR_* environment variables. For instance,
 * some environment variables could specify prefix and suffix strings
 * to colorize parts of the output that correspond to invalid data.
 *
 * Moreover, get rid of mpfr_print_binary? It should be replaced
 * either by mpfr_dump or by mpfr_out_str in the MPFR code.
 */

void
mpfr_dump (mpfr_srcptr u)
{
  mpfr_print_binary(u);
  putchar('\n');
}
