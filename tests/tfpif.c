/* Test file for mpfr_fpif.

Copyright 2012 Free Software Foundation, Inc.
Contributed by Olivier Demengeon.

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

#include "mpfr-test.h"

#define FILE_NAME_RW "mpfrtest.txt" /* temporary name (written then read) */
#define FILE_NAME_R  "mpfrtest.dat" /* fixed file name (read only) */

int
main (int argc, char *argv[])
{
  char *filenameCompressed = FILE_NAME_RW;
  char *data = FILE_NAME_R;
  int status;
  FILE *fh;
  mpfr_t x[9];
  mpfr_t y;
  int i;

  if (argc != 1)
    {
      printf ("Usage: %s\n", argv[0]);
      exit (1);
    }

  tests_start_mpfr ();

  mpfr_init2 (x[0], 130);
  mpfr_init2 (x[8], 130);
  mpfr_inits2 (2048, x[1], x[2], x[3], x[4], x[5], x[6], x[7], (mpfr_ptr) 0);
  mpfr_set_d (x[0], 45.2564215, MPFR_RNDN);
  mpfr_set_d (x[1], 45.2564215, MPFR_RNDN);
  mpfr_set_d (x[2], 45.2564215, MPFR_RNDN);
  mpfr_set_exp (x[2], -48000);
  mpfr_set_inf (x[3], -1);
  mpfr_set_zero (x[4], 0);
  mpfr_set_nan (x[5]);
  mpfr_set_d (x[6], 104348, MPFR_RNDN);
  mpfr_set_d (x[7], 33215, MPFR_RNDN);
  mpfr_div (x[8], x[6], x[7], MPFR_RNDN);
  mpfr_div (x[6], x[6], x[7], MPFR_RNDN);

  /* we first write to file FILE_NAME_RW the numbers x[i] */
  filenameCompressed = FILE_NAME_RW;
  fh = fopen (filenameCompressed, "w");
  if (fh == NULL)
    {
      printf ("Failed to open for writing %s, exiting...\n",
              filenameCompressed);
      exit (1);
    }

  for (i = 0; i < 9; i++)
    {
      status = mpfr_fpif_export (fh, x[i]);
      if (status != 0)
        {
          fclose (fh);
          printf ("Failed to export number %d, exiting...\n", i);
          exit (1);
        }
    }

  fclose (fh);

  /* we then read back FILE_NAME_RW and check we get the same numbers x[i] */
  fh = fopen (filenameCompressed, "r");
  if (fh == NULL)
    {
      printf ("Failed to open for reading %s, exiting...\n",
              filenameCompressed);
      exit (1);
    }

  for (i = 0; i < 9; i++)
    {
      mpfr_init2 (y, 2);
      mpfr_fpif_import (y, fh);
      if (mpfr_cmp(x[i], y) != 0)
        {
          printf ("mpfr_cmp failed on written number %d, exiting...\n", i);
          printf ("expected "); mpfr_dump (x[i]);
          printf ("got      "); mpfr_dump (y);
          exit (1);
        }
      mpfr_clear (y);
    }
  fclose (fh);

  /* we do the same for the fixed file FILE_NAME_R, this ensures
     we get same results with different word size or endianness */
  fh = fopen (data, "r");
  if (fh == NULL)
    {
      printf ("Failed to open for reading %s, exiting...\n", data);
      exit (1);
    }

  for (i = 0; i < 9; i++)
    {
      mpfr_init2 (y, 2);
      mpfr_fpif_import (y, fh);
      if (mpfr_cmp (x[i], y) != 0)
        {
          printf ("mpfr_cmp failed on data number %d, exiting...\n", i);
          printf ("expected "); mpfr_dump (x[i]);
          printf ("got      "); mpfr_dump (y);
          exit (1);
        }
      mpfr_clear (y);
    }
  fclose (fh);

  for (i = 0; i < 9; i++)
    mpfr_clear (x[i]);

  remove (FILE_NAME_RW);
  tests_end_mpfr ();
  return 0;
}
