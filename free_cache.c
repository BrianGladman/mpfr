/* mpfr_free_cache - Free the cache used by MPFR for internal consts.

Copyright 2004 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include "mpfr-impl.h"

void
mpfr_free_cache (void)
{
  if (__gmpfr_const_pi_prec != 0)
    {
      mpfr_clear (__mpfr_const_pi);
      __gmpfr_const_pi_prec = 0;
    }

  if (__gmpfr_const_log2_prec != 0)
    {
      mpfr_clear (__mpfr_const_log2);
      __gmpfr_const_log2_prec = 0;
    }
}
