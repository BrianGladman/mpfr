/* Test file for mpfr_dump.

Copyright (C) 2000 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdlib.h>
#include <unistd.h>
#include "gmp.h"
#include "mpfr.h"

int main()
{
  mpfr_t z;

  mpfr_init2(z, 100);
  mpfr_set_ui(z, 0, GMP_RNDN);
  mpfr_dump(z, GMP_RNDD);
  printf("   ^--- 0.e1 printed above is ok\n");
  mpfr_clear(z);
  return 0;
}

