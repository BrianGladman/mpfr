/* functions to set/get the machine rounding mode,
   this should use the interface from the ISO C9X draft
   (file <fenv.h>) */

#include <stdio.h>
#include "gmp.h"
#include "mpfr.h"

#ifdef IRIX64
#include <sys/fpu.h>
extern int swapRM();
#define TOZERO swapRM(ROUND_TO_ZERO)
#define TOINFP swapRM(ROUND_TO_PLUS_INFINITY)
#define TONEAREST swapRM(ROUND_TO_NEAREST)
#define TOINFM swapRM(ROUND_TO_MINUS_INFINITY)
#elif (defined (solaris) || defined (sun4))
#include <ieeefp.h>
#define TOZERO fpsetround(FP_RZ)
#define TOINFP fpsetround(FP_RP)
#define TONEAREST fpsetround(FP_RN)
#define TOINFM fpsetround(FP_RM)
#elif alpha
#ifdef __GNUC__
/* GCC patched include files forget to define those... */
#define FP_RND_RZ       0
#define FP_RND_RN       1
#define FP_RND_RP       2
#define FP_RND_RM       3
#endif
#include <float.h>
#define TOZERO write_rnd(FP_RND_RZ)
#define TOINFP write_rnd(FP_RND_RP)
#define TONEAREST write_rnd(FP_RND_RN)
#define TOINFM write_rnd(FP_RND_RM)
#elif AIX
#include <float.h>
#define TOZERO fp_swap_rnd(FP_RND_RZ)
#define TOINFP fp_swap_rnd(FP_RND_RP)
#define TONEAREST fp_swap_rnd(FP_RND_RN)
#define TOINFM fp_swap_rnd(FP_RND_RM)
#elif sunos
#include <floatingpoint.h>
char *out;
#define TOZERO ieee_flags("set","direction","tozero",&out)
#define TOINFP ieee_flags("set","direction","positive",&out)
#define TONEAREST ieee_flags("set","direction","nearest",&out)
#define TOINFM ieee_flags("set","direction","negative",&out)
#elif hp700
#define TOZERO fpsetround(FP_RZ)
#define TOINFP fpsetround(FP_RP)
#define TONEAREST fpsetround(FP_RN)
#define TOINFM fpsetround(FP_RM)
#elif (defined (__i386__) || defined (__i486__) || defined (linux))
#include <fpu_control.h>
#ifdef LIBC211
#define __setfpucw(cw) __asm__ ("fldcw %0" : : "m" (cw))
#endif
/* be careful to put bits 9-8 (Precision control) 
   to 10 (round to double precision) instead of the
   default 11 (round to extended precision) */
#define TOZERO __setfpucw(0x1e7f)
#define TOINFP __setfpucw(0x1a7f)
#define TOINFM __setfpucw(0x167f)
#define TONEAREST __setfpucw(0x127f)
#endif

/* sets the machine rounding mode to the value rnd_mode */
void 
#if __STDC__
mpfr_set_machine_rnd_mode(unsigned char rnd_mode)
#else
mpfr_set_machine_rnd_mode(rnd_mode)
     unsigned char rnd_mode; 
#endif
{
  switch (rnd_mode) {
  case GMP_RNDN: TONEAREST; break;
  case GMP_RNDZ: TOZERO; break;
  case GMP_RNDU: TOINFP; break;
  case GMP_RNDD: TOINFM; break;
  default: fprintf(stderr,"invalid rounding mode\n"); exit(1);
  }
}

