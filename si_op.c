/* mpfr_add_si -- add a floating-point number with a machine integer
   mpfr_sub_si -- sub  a floating-point number with a machine integer
   mpfr_si_sub -- sub  a machine number with a floating-point number

Copyright 2004, 2005, 2006 Free Software Foundation, Inc.

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

#include "mpfr-impl.h"

int
mpfr_add_si (mpfr_ptr y, mpfr_srcptr x, long int u, mp_rnd_t rnd_mode)
{
  if (u >= 0)
    return mpfr_add_ui (y, x, u, rnd_mode);
  else
    return mpfr_sub_ui (y, x, -u, rnd_mode);
}

int
mpfr_sub_si (mpfr_ptr y, mpfr_srcptr x, long int u, mp_rnd_t rnd_mode)
{
  if (u >= 0)
    return mpfr_sub_ui (y, x, u, rnd_mode);
  else
    return mpfr_add_ui (y, x, -u, rnd_mode);
}

int
mpfr_si_sub (mpfr_ptr y, long int u, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  if (u >= 0)
    return mpfr_ui_sub (y, u, x, rnd_mode);
  else
    {
    int res = -mpfr_add_ui (y, x, -u, MPFR_INVERT_RND (rnd_mode));
    MPFR_CHANGE_SIGN (y);
    return res;
    }
}

