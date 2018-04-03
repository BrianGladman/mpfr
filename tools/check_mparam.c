/* to check some mparam.h table:
   1) make a symbolic link to the corresponding mparam.h or
      provide -DMPARAM='"..."' with a path to the mparam.h
      file when compiling this program
   2) compile and run this program */

#include <stdio.h>

#ifndef MPARAM
# define MPARAM "mparam.h"
#endif

#include MPARAM

#define numberof_const(x)  (sizeof (x) / sizeof ((x)[0]))

static short mulhigh_ktab[] = {MPFR_MULHIGH_TAB};
#define MPFR_MULHIGH_TAB_SIZE (numberof_const (mulhigh_ktab))

static short sqrhigh_ktab[] = {MPFR_SQRHIGH_TAB};
#define MPFR_SQRHIGH_TAB_SIZE (numberof_const (sqrhigh_ktab))

static short divhigh_ktab[] = {MPFR_DIVHIGH_TAB};
#define MPFR_DIVHIGH_TAB_SIZE (numberof_const (divhigh_ktab))

int main (void)
{
  int err = 0, n;

  for (n = 0; n < MPFR_MULHIGH_TAB_SIZE; n++)
    if (mulhigh_ktab[n] >= n)
      {
        printf ("Error, mulhigh_ktab[%d] = %d\n", n, mulhigh_ktab[n]);
        err = 1;
      }

  for (n = 0; n < MPFR_SQRHIGH_TAB_SIZE; n++)
    if (sqrhigh_ktab[n] >= n)
      {
        printf ("Error, sqrhigh_ktab[%d] = %d\n", n, sqrhigh_ktab[n]);
        err = 1;
      }

  for (n = 2; n < MPFR_DIVHIGH_TAB_SIZE; n++)
    if (divhigh_ktab[n] >= n-1)
      {
        printf ("Error, divhigh_ktab[%d] = %d\n", n, divhigh_ktab[n]);
        err = 1;
      }

  return err;
}
