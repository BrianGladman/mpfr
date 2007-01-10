/* mpfr_set_default_rounding_mode -- set the default rounding mode
   mpfr_get_default_rounding_mode -- get the default rounding mode

Copyright 1999, 2001, 2004, 2005, 2006, 2007 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include "mpfr-impl.h"

mp_rnd_t MPFR_THREAD_ATTR __gmpfr_default_rounding_mode = GMP_RNDN;

void
mpfr_set_default_rounding_mode (mp_rnd_t rnd_mode)
{
  if (rnd_mode >= GMP_RNDN && rnd_mode <= GMP_RNDD)
    __gmpfr_default_rounding_mode = rnd_mode;
}

#undef mpfr_get_default_rounding_mode
mp_rnd_t
mpfr_get_default_rounding_mode (void)
{
  return __gmpfr_default_rounding_mode;
}
