/*
Copyright 2005, 2006, 2007, 2008, 2009 Free Software Foundation, Inc.
Contributed by Patrick Pelissier, INRIA.

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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "timp.h"

/* To avoid __gmpf_cmp to be declared as pure */
#define __GMP_NO_ATTRIBUTE_CONST_PURE
#include "gmp.h"
#include "mpfr.h"

#ifdef SCS_SUPPORT
# define SCS(x) x
# include "scs.h"
# define EXTRA_TEST_LIST \
 BENCH("SCSLIB( Compiled var ):::::::::::", ; ); \
 BENCH("scs_add", scs_add(sc1,sc2,sc3)); \
 BENCH("scs_sub", scs_sub(sc1,sc2,sc3)); \
 BENCH("scs_mul", scs_mul(sc1,sc2,sc3)); \
 BENCH("scs_div", scs_div(sc1,sc2,sc3)); \
 BENCH("scs_set", scs_set(sc1,sc2)); \
 BENCH("scs_set0", scs_set_si(sc1,0)); \
 BENCH("scs_set1", scs_set_si(sc1,1));
#else
# define SCS(x) ((void) 0)
# define EXTRA_TEST_LIST ((void)0);
#endif

#undef EXTRA_TEST_LIST
# define EXTRA_TEST_LIST \
  BENCH("mpfr_exp", mpfr_exp(a,b,GMP_RNDN)); \
  BENCH("mpfr_log", mpfr_log(a,b,GMP_RNDN)); \
  BENCH("mpfr_sin", mpfr_sin(a,b,GMP_RNDN)); \
  BENCH("mpfr_cos", mpfr_cos(a,b,GMP_RNDN)); \
  BENCH("mpfr_tan", mpfr_tan(a,b,GMP_RNDN)); \
  BENCH("mpfr_asin", mpfr_asin(a,b,GMP_RNDN)); \
  BENCH("mpfr_acos", mpfr_acos(a,b,GMP_RNDN)); \
  BENCH("mpfr_atan", mpfr_atan(a,b,GMP_RNDN)); \
  BENCH("mpfr_agm", mpfr_agm(a,b,c,GMP_RNDN)); \
  BENCH("mpfr_const_log2", (mpfr_const_log2) (a, GMP_RNDN)); \
  BENCH("mpfr_const_pi", (mpfr_const_pi)(a, GMP_RNDN)); \
  BENCH("mpfr_sinh", mpfr_sinh(a,b,GMP_RNDN)); \
  BENCH("mpfr_cosh", mpfr_cosh(a,b,GMP_RNDN)); \
  BENCH("mpfr_tanh", mpfr_tanh(a,b,GMP_RNDN)); \
  BENCH("mpfr_asinh", mpfr_asinh(a,b,GMP_RNDN)); \
  BENCH("mpfr_acosh", mpfr_acosh(a,b,GMP_RNDN)); \
  BENCH("mpfr_atanh", mpfr_atanh(a,b,GMP_RNDN)); 



/* Theses macros help the compiler to determine if a test is likely*/
/* or unlikely. */
#if __GNUC__ >= 3
# define LIKELY(x)   (__builtin_expect(!!(x),1))
# define UNLIKELY(x) (__builtin_expect((x),0))
#else
# define LIKELY(x)   (x)
# define UNLIKELY(x) (x)
#endif

/*
 * List of all the tests to do.
 * Macro "Bench" is defined futhermore
 */
#define TEST_LIST \
  BENCH("MPFR::::::::::", ; ); \
  BENCH("mpfr_add", mpfr_add(a,b,c,GMP_RNDN)); \
  BENCH("mpfr_sub", mpfr_sub(a,b,c,GMP_RNDN)); \
  BENCH("mpfr_mul", mpfr_mul(a,b,c,GMP_RNDN)); \
  BENCH("mpfr_div", mpfr_div(a,b,c,GMP_RNDN)); \
  BENCH("mpfr_sqrt", mpfr_sqrt(a,b,GMP_RNDN)); \
  BENCH("mpfr_cmp", mpfr_cmp(b,c)); \
  BENCH("mpfr_sgn", mpfr_sgn(b)); \
  BENCH("mpfr_set", mpfr_set(a,b, GMP_RNDN)); \
  BENCH("mpfr_set0", mpfr_set_si(a,0,GMP_RNDN)); \
  BENCH("mpfr_set1", mpfr_set_si(a,1,GMP_RNDN)); \
  BENCH("mpfr_swap", mpfr_swap(b,c)); \
  BENCH("MPF:::::::::::", ; ); \
  BENCH("mpf_add", mpf_add(x,y,z)); \
  BENCH("mpf_sub", mpf_sub(x,y,z)); \
  BENCH("mpf_mul", mpf_mul(x,y,z)); \
  BENCH("mpf_div", mpf_div(x,y,z)); \
  BENCH("mpf_sqrt", mpf_sqrt(x,y)); \
  BENCH("mpf_cmp", mpf_cmp(y,z)); \
  BENCH("mpf_set", mpf_set(x,y)); \
  BENCH("mpf_set0", mpf_set_si(x,0)); \
  BENCH("mpf_set1", mpf_set_si(x,1)); \
  BENCH("mpf_swap", mpf_swap(y,z)); \
  EXTRA_TEST_LIST

#define USAGE                                                                \
 "Get the graph of the low-level functions of Mpfr (gnuplot).\n"             \
  __FILE__" " __DATE__" " __TIME__" GCC "__VERSION__ "\n"                    \
 "Usage: mpfr-gfx [-bPREC_BEGIN] [-ePREC_END] [-sPREC_STEP] [-mSTAT_SIZE] \n"\
 "       [-oFILENAME] [-xFUNCTION_NUM] [-yFUNCTION_NUM] [-c] [-fSMOOTH]\n"

unsigned long num;
mpf_t *xt, *yt, *zt;
int smooth = 3;

void lets_start(unsigned long n, mp_prec_t p)
{
  unsigned long i;
  gmp_randstate_t state;

  num = n;
  xt = malloc(sizeof(mpf_t) * num);
  yt = malloc(sizeof(mpf_t) * num);
  zt = malloc(sizeof(mpf_t) * num);
  if (xt==NULL || yt==NULL || zt==NULL)
    {
      fprintf(stderr, "Can't allocate tables!\n");
      abort();
    }

  gmp_randinit_lc_2exp_size (state, 128);
  gmp_randseed_ui (state, 1452369);
  for(i = 0 ; i < num ; i++)
    {
      mpf_init2(xt[i], p);
      mpf_init2(yt[i], p);
      mpf_init2(zt[i], p);
      mpf_urandomb(yt[i], state, p);
      mpf_urandomb(zt[i], state, p);
    }  
  gmp_randclear(state);
}

void lets_end(void)
{
  unsigned long i;

  for(i = 0 ; i < num ; i++)
    {
      mpf_clear(xt[i]);
      mpf_clear(yt[i]);
      mpf_clear(zt[i]);
    }
  free (xt);
  free (yt);
  free (zt);
}

double get_speed(mp_prec_t p, int select)
{
  unsigned long long mc[num], m;
  mpfr_t a,b,c;
  mpf_t x,y,z;
  unsigned long long moy;
  int i,j=0, op, cont, print_done = 0;
  const char *str = "void";
  SCS(( scs_t sc1, sc2, sc3 ));

  mpf_init2(x, p);  mpf_init2(y, p);  mpf_init2(z, p);
  mpfr_init2(a, p); mpfr_init2(b, p); mpfr_init2(c, p);

  for(i = 0 ; i < num ; i++)
    {
      //      yt[i][0]._mp_exp = (rand() % p) / GMP_NUMB_BITS;
      //zt[i][0]._mp_exp = (rand() % p) / GMP_NUMB_BITS;      
      mc[i] = 0xFFFFFFFFFFFFFFFLL;
    }

  TIMP_OVERHEAD ();

  for(j = 0, cont = smooth ; cont ; j++, cont--)
    {
      for(i = 0 ; i < num ; i++)
	{
	  /* Set var for tests */
	  mpf_set(y, yt[i]);
	  mpf_set(z, zt[i]);
	  mpfr_set_f(b, yt[i], GMP_RNDN);
	  mpfr_set_f(c, zt[i], GMP_RNDN);
	  SCS(( scs_set_mpfr(sc2, b), scs_set_mpfr(sc3, c) ));
#undef BENCH
#define BENCH(TEST_STR, TEST)                           \
	  if (op++ == select) {                         \
	      m = TIMP_MEASURE(TEST);                   \
              str = TEST_STR;                           \
	      if (m < mc[i]) {mc[i] = m; cont = smooth;}\
	    }
	  op = 0;
	  TEST_LIST;
	  if (print_done == 0 && strcmp (str, "void") != 0 )
	    {
	      printf("Prec=%4.4lu Func=%20.20s", p, str);
	      fflush (stdout);
	      print_done = 1;
	    }
	}
    }
  mpfr_clear(a);  mpfr_clear(b);  mpfr_clear(c);
  mpf_clear(x);   mpf_clear(y);   mpf_clear(z);
  /* End */
  /* Save result */
  moy = mc[0];
  for(i = 1 ; i < num ; i++) moy += mc[i];
  printf(" Pass=%4.4d..................%Lu.%Lu\n",
	 j+1, moy/num, (moy*100LL/num)%100LL);
  return (double) (moy) / (double) num;
}

int write_data(const char *filename, 
	       unsigned long num,
	       mp_prec_t p1, mp_prec_t p2, mp_prec_t ps,
	       int select1, int select2)
{
  FILE *f;
  mp_prec_t p;

  lets_start (num, p2);
  f = fopen(filename, "w");
  if (f==NULL)
    {
      fprintf(stderr, "Can't open %s!\n", filename);
      lets_end ();
      abort();
    }
  for (p = p1 ; p < p2 ; p+=ps)
    fprintf(f, "%lu\t%1.20e\t%1.20e\n", p, 
	    get_speed(p, select1),
	    get_speed(p, select2));
  fclose(f);
  lets_end ();
  return 0;
}

int write_data2(const char *filename, 
		unsigned long num,
		mp_prec_t p_begin, mp_prec_t p_end, mp_prec_t p_step,
		int s_begin, int s_end)
{
  FILE *f;
  mp_prec_t p;
  int s;

  lets_start (num, p_end);
  f = fopen (filename, "w");
  if (f == NULL)
    {
      fprintf (stderr, "Can't open %s!\n", filename);
      lets_end ();
      exit (1);
    }
  for (p = p_begin ; p < p_end ; p += p_step)
    {
      fprintf (f, "%lu", p);
      for (s = s_begin ; s <= s_end ; s++)
	fprintf (f, "\t%1.20e", get_speed (p, s));
      fprintf (f, "\n");
    }
  fclose (f);
  lets_end ();
  return 0;
}

int op_num (void)
{
  int op;
#undef BENCH
#define BENCH(TEST_STR, TEST) op++;
  op = 0;
  TEST_LIST;
  return op;
}

int main(int argc, const char *argv[])
{
  mp_prec_t p1, p2, ps;
  int i;
  unsigned long stat;
  const char *filename = "plot.data";
  int select1, select2, max_op, conti;

  printf (USAGE);

  max_op = op_num ();
  select1 = 1; select2 = 13;
  p1 = 2; p2 = 500; ps = 4;
  stat = 500;
  conti = 0;

  for(i = 1 ; i < argc ; i++)
    {
      if (argv[i][0] == '-')
	{
	  switch (argv[i][1])
	    {
	    case 'b':
	      p1 = atol(argv[i]+2);
	      break;
	    case 'e':
	      p2 = atol(argv[i]+2);
	      break;
            case 's':
              ps = atol(argv[i]+2);
              break;
 	    case 'm':
	      stat = atol(argv[i]+2);
	      break;
	    case 'x':
	      select1 = atoi (argv[i]+2);
	      select2 = select1 + 12;
	      break;
            case 'y':
              select2 = atoi (argv[i]+2);
              break;
	    case 'o':
	      filename = argv[i]+2;
	      break;
	    case 'c':
	      conti = 1;
	      break;
	    case 'f':
	      smooth = atoi  (argv[i]+2);
	      break;
	    default:
	      fprintf(stderr, "Unkwown option: %s\n", argv[i]);
	      abort ();
	    }
	}
    }
  /* Set low priority */
  setpriority(PRIO_PROCESS,0,14);
  printf("GMP:%s MPFR:%s From p=%lu to %lu by %lu Output: %s N=%ld\n",
	 gmp_version, mpfr_get_version(), p1,p2,ps, filename,stat);
  if (select2 >= max_op)
    select2 = max_op-1;
  if (select1 >= max_op)
    select1 = max_op-1;

  if (conti == 0)
    write_data (filename, stat, p1, p2, ps, select1, select2);
  else
    write_data2 (filename, stat, p1, p2, ps, select1, select2);

  return 0;
}
