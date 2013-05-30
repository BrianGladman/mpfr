/*
Copyright 2005-2009 Free Software Foundation, Inc.
Contributed by Patrick Pelissier, INRIA.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#ifndef __TIMP__H__
#define __TIMP__H__

/* Usage:
 *  Before doing the measure, call TIMP_OVERHEAD ();
 *  Then unsigned long long t = TIMP_MEASURE (f(x));
 *  to measure the # of cycles taken by the call to f(x).
 */

#define TIMP_VERSION 1*100+0*10+0

#ifndef __GNUC__
# error  CC != GCC 
#endif

/* High accuracy timing */
#if defined (USE_CLOCK_MONOTONIC)

/* Needs to include -lrt in the library section */
#include <time.h> 

#define timp_rdtsc()                                           \
({unsigned long long int x;				       \
  struct timespec ts;                                          \
  clock_gettime(CLOCK_MONOTONIC, &ts);                         \
  x = ts.tv_sec * 1000000000UL + ts.tv_nsec;                   \
 x; })
#define timp_rdtsc_before(time) (time = timp_rdtsc())
#define timp_rdtsc_after(time)  (time = timp_rdtsc())

#elif defined (__i386__) || defined(__amd64__)

#define timp_rdtsc_before(time)           \
        __asm__ __volatile__(             \
                ".align 64\n\t"           \
                "xorl %%eax,%%eax\n\t"    \
                "cpuid\n\t"               \
                "rdtsc\n\t"               \
                "movl %%eax,(%0)\n\t"     \
                "movl %%edx,4(%0)\n\t"    \
                "xorl %%eax,%%eax\n\t"    \
                "cpuid\n\t"               \
                : /* no output */         \
                : "S"(&time)              \
                : "eax", "ebx", "ecx", "edx", "memory")

#define timp_rdtsc_after(time)            \
        __asm__ __volatile__(             \
                "xorl %%eax,%%eax\n\t"    \
                "cpuid\n\t"               \
                "rdtsc\n\t"               \
                "movl %%eax,(%0)\n\t"     \
                "movl %%edx,4(%0)\n\t"    \
                "xorl %%eax,%%eax\n\t"    \
                "cpuid\n\t"               \
                : /* no output */         \
                : "S"(&time)              \
                : "eax", "ebx", "ecx", "edx", "memory")

#elif defined (__ia64)

#define timp_rdtsc()                                           \
({ unsigned long long int x;                                   \
  __asm__ __volatile__("mov %0=ar.itc" : "=r"(x) :: "memory"); \
  x; })
#define timp_rdtsc_before(time) (time = timp_rdtsc())
#define timp_rdtsc_after(time)  (time = timp_rdtsc())

#elif defined (__alpha)

#define timp_rdtsc()                              \
({ unsigned long long int x;                      \
   __asm__ volatile ("rpcc %0\n\t" : "=r" (x));   \
   x; })
#define timp_rdtsc_before(time) (time = tpp_rdtsc())
#define timp_rdtsc_after(time)  (time = tpp_rdtsc())

#else
# error Unsupported CPU
#endif

/* We do several measures and keep the minimum to avoid counting
 * hardware interruption cycles.
 * The filling of the CPU cache is done because we do several loops,
 * and get the minimum.
 * Declaring num_cycle as "volatile" is to avoid optimisation when it is
 * possible (To properly calcul overhead).
 * overhead is calculated outside by a call to:
 *   overhead = MEASURE("overhead", ;)
 * Use a lot the preprocessor.
 * It is a macro to be very flexible.
 */
static unsigned long long int timp_overhead = 0;

#define TIMP_NUM_TRY  4327
#define TIMP_MAX_WAIT_FOR_MEASURE 10000000ULL

#define TIMP_MEASURE(CODE)                                            \
  ({                                                                  \
  volatile unsigned long long int num_cycle, num_cycle2;              \
  unsigned long long min_num_cycle, start_num_cycle;                  \
  int _i;                                                             \
  timp_rdtsc_before (start_num_cycle);                                \
  min_num_cycle = 0xFFFFFFFFFFFFFFFFLL;                               \
  for(_i = 0 ; _i < TIMP_NUM_TRY ; _i++) {                            \
    timp_rdtsc_before(num_cycle);                                     \
    CODE;                                                             \
    timp_rdtsc_after(num_cycle2);                                     \
    num_cycle =  num_cycle2 - num_cycle;                              \
    if (num_cycle < min_num_cycle)                                    \
      min_num_cycle = num_cycle;                                      \
    if (num_cycle2 - start_num_cycle > TIMP_MAX_WAIT_FOR_MEASURE)     \
      break;                                                          \
  }                                                                   \
  min_num_cycle - timp_overhead; })

#define TIMP_OVERHEAD()                                               \
  (timp_overhead = 0, timp_overhead = TIMP_MEASURE((void) 0) )

#endif /* __TIMP__H__ */
