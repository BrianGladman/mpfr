/* Test file for mpfr_version.

Copyright 2004, 2005 Free Software Foundation.

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
#include <string.h>

#include "mpfr-test.h"

int
main (void)
{
  char buffer[16];
  const char *version;

  sprintf (buffer, "%d.%d.%d", MPFR_VERSION_MAJOR, MPFR_VERSION_MINOR,
           MPFR_VERSION_PATCHLEVEL);
  version = mpfr_get_version ();
  if (strcmp (buffer, version) != 0)
    {
      printf ("Incorrect version (%s vs %s)\n", buffer, version);
      exit (1);
    }
  if (strcmp (MPFR_VERSION_STRING, version) != 0)
    {
      printf ("Incorrect version (%s vs %s)\n", buffer, version);
      exit (1);
    }
  return 0;
}
