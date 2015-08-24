/* Test file for:
 mpfr_fits_sint_p, mpfr_fits_slong_p, mpfr_fits_sshort_p,
 mpfr_fits_uint_p, mpfr_fits_ulong_p, mpfr_fits_ushort_p

Copyright 2004-2015 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

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

#ifdef HAVE_CONFIG_H
# include "config.h"       /* for a build within gmp */
#endif

#include "mpfr-intmax.h"
#include "mpfr-test.h"

#define ERROR1(N)                                               \
  do                                                            \
    {                                                           \
      printf("Error %d for rnd = %s and x = ", N,               \
             mpfr_print_rnd_mode ((mpfr_rnd_t) r));             \
      mpfr_dump(x);                                             \
      exit(1);                                                  \
    }                                                           \
  while (0)

#define FTEST(N,NOT,FCT)                                        \
  do                                                            \
    {                                                           \
      __gmpfr_flags = ex_flags;                                 \
      if (NOT FCT (x, (mpfr_rnd_t) r))                          \
        {                                                       \
          printf ("Error %d for %s, rnd = %s and x = ",         \
                  N, #FCT,                                      \
                  mpfr_print_rnd_mode ((mpfr_rnd_t) r));        \
          mpfr_dump (x);                                        \
          exit (1);                                             \
        }                                                       \
      if (__gmpfr_flags != ex_flags)                            \
        {                                                       \
          mpfr_flags_t flags = __gmpfr_flags;                   \
          printf ("Flags error %d for %s, rnd = %s and x = ",   \
                  N, #FCT,                                      \
                  mpfr_print_rnd_mode ((mpfr_rnd_t) r));        \
          mpfr_dump(x);                                         \
          printf ("Expected flags:");                           \
          flags_out (ex_flags);                                 \
          printf ("Got flags:     ");                           \
          flags_out (flags);                                    \
          exit (1);                                             \
        }                                                       \
    }                                                           \
  while (0)

#define CHECK_ALL(N,NOT)                                        \
  do                                                            \
    {                                                           \
      FTEST (N, NOT, mpfr_fits_ulong_p);                        \
      FTEST (N, NOT, mpfr_fits_slong_p);                        \
      FTEST (N, NOT, mpfr_fits_uint_p);                         \
      FTEST (N, NOT, mpfr_fits_sint_p);                         \
      FTEST (N, NOT, mpfr_fits_ushort_p);                       \
      FTEST (N, NOT, mpfr_fits_sshort_p);                       \
    }                                                           \
  while (0)

static void check_intmax (void);

int
main (void)
{
  mpfr_t x, y;
  mpfr_flags_t flags[2] = { 0, MPFR_FLAGS_ALL }, ex_flags;
  int i, r, fi;

  tests_start_mpfr ();

  mpfr_init2 (x, 256);
  mpfr_init2 (y, 8);

  RND_LOOP (r)
    for (fi = 0; fi < numberof (flags); fi++)
      {
        ex_flags = flags[fi];

        /* Check NaN */
        mpfr_set_nan (x);
        CHECK_ALL (1, !!);

        /* Check +Inf */
        mpfr_set_inf (x, 1);
        CHECK_ALL (2, !!);

        /* Check -Inf */
        mpfr_set_inf (x, -1);
        CHECK_ALL (3, !!);

        /* Check +0 */
        mpfr_set_zero (x, 1);
        CHECK_ALL (4, !);

        /* Check -0 */
        mpfr_set_zero (x, -1);
        CHECK_ALL (5, !);

        /* Check small positive op */
        mpfr_set_str1 (x, "1@-1");
        CHECK_ALL (6, !);

        /* Check 17 */
        mpfr_set_ui (x, 17, MPFR_RNDN);
        CHECK_ALL (7, !);

        /* Check large values (no fit) */
        mpfr_set_ui (x, ULONG_MAX, MPFR_RNDN);
        mpfr_mul_2exp (x, x, 1, MPFR_RNDN);
        CHECK_ALL (8, !!);
        mpfr_mul_2exp (x, x, 40, MPFR_RNDN);
        CHECK_ALL (9, !!);

        /* Check min/max of the types (except 0 for unsigned types) */
        mpfr_set_ui (x, ULONG_MAX, MPFR_RNDN);
        FTEST (10, !, mpfr_fits_ulong_p);
        mpfr_set_si (x, LONG_MAX, MPFR_RNDN);
        FTEST (11, !, mpfr_fits_slong_p);
        mpfr_set_si (x, LONG_MIN, MPFR_RNDN);
        FTEST (12, !, mpfr_fits_slong_p);
        mpfr_set_ui (x, UINT_MAX, MPFR_RNDN);
        FTEST (13, !, mpfr_fits_uint_p);
        mpfr_set_si (x, INT_MAX, MPFR_RNDN);
        FTEST (14, !, mpfr_fits_sint_p);
        mpfr_set_si (x, INT_MIN, MPFR_RNDN);
        FTEST (15, !, mpfr_fits_sint_p);
        mpfr_set_ui (x, USHRT_MAX, MPFR_RNDN);
        FTEST (16, !, mpfr_fits_ushort_p);
        mpfr_set_si (x, SHRT_MAX, MPFR_RNDN);
        FTEST (17, !, mpfr_fits_sshort_p);
        mpfr_set_si (x, SHRT_MIN, MPFR_RNDN);
        FTEST (18, !, mpfr_fits_sshort_p);

        /* Check negative op */
        for (i = 1; i <= 4; i++)
          {
            int inv;

            mpfr_set_si_2exp (x, -i, -2, MPFR_RNDN);
            mpfr_rint (y, x, (mpfr_rnd_t) r);
            inv = MPFR_NOTZERO (y);
            FTEST (19, inv ^ !, mpfr_fits_ulong_p);
            FTEST (19,       !, mpfr_fits_slong_p);
            FTEST (19, inv ^ !, mpfr_fits_uint_p);
            FTEST (19,       !, mpfr_fits_sint_p);
            FTEST (19, inv ^ !, mpfr_fits_ushort_p);
            FTEST (19,       !, mpfr_fits_sshort_p);
          }
      }

  mpfr_clear (x);
  mpfr_clear (y);

  check_intmax ();

  tests_end_mpfr ();
  return 0;
}

static void
check_intmax (void)
{
#ifdef _MPFR_H_HAVE_INTMAX_T
  mpfr_t x, y;
  int i, r;

  mpfr_init2 (x, sizeof (uintmax_t) * CHAR_BIT);
  mpfr_init2 (y, 8);

  RND_LOOP (r)
    {
      /* Check NAN */
      mpfr_set_nan (x);
      if (mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (52);
      if (mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (53);

      /* Check INF */
      mpfr_set_inf (x, 1);
      if (mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (54);
      if (mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (55);

      /* Check Zero */
      MPFR_SET_ZERO (x);
      if (!mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (56);
      if (!mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (57);

      /* Check positive small op */
      mpfr_set_str1 (x, "1@-1");
      if (!mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (58);
      if (!mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (59);

      /* Check 17 */
      mpfr_set_ui (x, 17, MPFR_RNDN);
      if (!mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (60);
      if (!mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (61);

      /* Check hugest */
      mpfr_set_ui_2exp (x, 42, sizeof (uintmax_t) * 32, MPFR_RNDN);
      if (mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (62);
      if (mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (63);

      /* Check all other values */
      mpfr_set_uj (x, MPFR_UINTMAX_MAX, MPFR_RNDN);
      mpfr_add_ui (x, x, 1, MPFR_RNDN);
      if (mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (64);
      mpfr_set_uj (x, MPFR_UINTMAX_MAX, MPFR_RNDN);
      if (!mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (65);
      mpfr_set_sj (x, MPFR_INTMAX_MAX, MPFR_RNDN);
      mpfr_add_ui (x, x, 1, MPFR_RNDN);
      if (mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (66);
      mpfr_set_sj (x, MPFR_INTMAX_MAX, MPFR_RNDN);
      if (!mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (67);
      mpfr_set_sj (x, MPFR_INTMAX_MIN, MPFR_RNDN);
      if (!mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (68);
      mpfr_sub_ui (x, x, 1, MPFR_RNDN);
      if (mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
        ERROR1 (69);

      /* Check negative op */
      for (i = 1; i <= 4; i++)
        {
          int inv;

          mpfr_set_si_2exp (x, -i, -2, MPFR_RNDN);
          mpfr_rint (y, x, (mpfr_rnd_t) r);
          inv = MPFR_NOTZERO (y);
          if (!mpfr_fits_uintmax_p (x, (mpfr_rnd_t) r) ^ inv)
            ERROR1 (70);
          if (!mpfr_fits_intmax_p (x, (mpfr_rnd_t) r))
            ERROR1 (71);
        }
    }

  mpfr_clear (x);
  mpfr_clear (y);
#endif
}
