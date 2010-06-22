#include <stdlib.h>
#include <time.h>

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#undef _PROTO
#define _PROTO __GMP_PROTO
#include "speed.h"

#define SMALL_PRECISION 32

/* Let f be a function for which we have several implementations f1, f2... */
/* We wish to have a quick overview of which implementation is the best    */
/* in function of the point x where f(x) is computed and of the prectision */
/* prec requested by the user.                                             */
/* This is performed by drawing a 2D graphic with color indicating which   */
/* method is the best.                                                     */
/* For building this graphic, the following structur is used (see the core */
/* of generate_2D_sample for an explanation of each field.                 */
struct speed_params2D
{
  /* x-window: [min_x, max_x] or [2^min_x, 2^max_x]                        */
  /*           or [-2^(max_x), -2^(min_x)] U [2^min_x, 2^max_x]            */
  /* depending on the value of logarithmic_scale_x                         */
  double min_x;
  double max_x;

  /* prec-window: [min_prec, max_prec] */
  mpfr_prec_t min_prec;
  mpfr_prec_t max_prec;

  /* number of sample points for the x-axis and the prec-axis */
  int nb_points_x;
  int nb_points_prec;

  /* should the sample points be logarithmically scaled or not */
  int logarithmic_scale_x;
  int logarithmic_scale_prec;

  /* list of functions g1, g2... measuring the execution time of f1, f2...  */
  /* when given a point x and a precision prec stored in s.                 */
  /* We use s->xp to store the significant of x, s->r to store its exponent */
  /* s->align_xp to store its sign, and s->size to store prec.              */
  double (**speed_funcs) (struct speed_params *s);
};

/* Given an array t of nb_functions double indicating the timings of several */
/* implementations, return i, such that t[i] is the best timing.             */
int
find_best (double *t, int nb_functions)
{
  int i, ibest;
  double best;

  if (nb_functions<=0)
    {
      fprintf (stderr, "There is no function\n");
      abort ();
    }

  ibest = 0;
  best = t[0];
  for (i=1; i<nb_functions; i++)
    {
      if (t[i]<best)
        {
          best = t[i];
          ibest = i;
        }
    }

  return ibest;
}

void generate_2D_sample (FILE *output, struct speed_params2D param)
{
  mpfr_t temp;
  double incr_prec;
  mpfr_t incr_x;
  mpfr_t x, x2;
  double prec;
  struct speed_params s;
  int i;
  int test;
  int nb_functions;
  double *t; /* store the timing of each implementation */

  /* We first determine how many implementations we have */
  nb_functions = 0;
  while (param.speed_funcs[nb_functions] != NULL)
    nb_functions++;

  t = malloc (nb_functions * sizeof (double));
  if (t == NULL)
    {
      fprintf (stderr, "Can't allocate memory.\n");
      abort ();
    }


  mpfr_init2 (temp, SMALL_PRECISION);

  /* The precision is sampled from min_prec to max_prec with        */
  /* approximately nb_points_prec points. If logarithmic_scale_prec */
  /* is true, the precision is multiplied by incr_prec at each      */
  /* step. Otherwise, incr_prec is added at each step.              */
  if (param.logarithmic_scale_prec)
    {
      mpfr_set_ui (temp, (unsigned long int)param.max_prec, MPFR_RNDU);
      mpfr_div_ui (temp, temp, (unsigned long int)param.min_prec, MPFR_RNDU);
      mpfr_root (temp, temp,
                 (unsigned long int)param.nb_points_prec, MPFR_RNDU);
      incr_prec = mpfr_get_d (temp, MPFR_RNDU);
    }
  else
    {
      incr_prec = (double)param.max_prec - (double)param.min_prec;
      incr_prec = incr_prec/((double)param.nb_points_prec);
    }

  /* The points x are sampled according to the following rule:             */
  /* If logarithmic_scale_x = 0:                                           */
  /*    nb_points_x points are equally distributed between min_x and max_x */
  /* If logarithmic_scale_x = 1:                                           */
  /*    nb_points_x points are sampled from 2^(min_x) to 2^(max_x). At     */
  /*    each step, the current point is multiplied by incr_x.              */
  /* If logarithmic_scale_x = -1:                                          */
  /*    nb_points_x/2 points are sampled from -2^(max_x) to -2^(min_x)     */
  /*    (at each step, the current point is divided by incr_x);  and       */
  /*    nb_points_x/2 points are sampled from 2^(min_x) to 2^(max_x)       */
  /*    (at each step, the current point is multiplied by incr_x).         */
  mpfr_init2 (incr_x, param.max_prec);
  if (param.logarithmic_scale_x == 0)
    {
      mpfr_set_d (incr_x,
                  (param.max_x - param.min_x)/(double)param.nb_points_x,
                  MPFR_RNDU);
    }
  else if (param.logarithmic_scale_x == -1)
    {
      mpfr_set_d (incr_x,
                  2.*(param.max_x - param.min_x)/(double)param.nb_points_x,
                  MPFR_RNDU);
      mpfr_exp2 (incr_x, incr_x, MPFR_RNDU);
    }
  else
    { /* other values of param.logarithmic_scale_x are considered as 1 */
      mpfr_set_d (incr_x,
                  (param.max_x - param.min_x)/(double)param.nb_points_x,
                  MPFR_RNDU);
      mpfr_exp2 (incr_x, incr_x, MPFR_RNDU);
    }

  /* Main loop */
  mpfr_init2 (x, param.max_prec);
  mpfr_init2 (x2, param.max_prec);
  prec = (double)param.min_prec;
  while (prec <= param.max_prec)
    {
      if (param.logarithmic_scale_x == 0)
        mpfr_set_d (temp, param.min_x, MPFR_RNDU);
      else if (param.logarithmic_scale_x == -1)
        {
          mpfr_set_d (temp, param.max_x, MPFR_RNDD);
          mpfr_exp2 (temp, temp, MPFR_RNDD);
          mpfr_neg (temp, temp, MPFR_RNDU);
        }
      else
        {
          mpfr_set_d (temp, param.min_x, MPFR_RNDD);
          mpfr_exp2 (temp, temp, MPFR_RNDD);
        }

      /* We perturb x a little bit, in order to avoid trailing zeros that */
      /* might change the behavior of algorithms.                         */
      mpfr_const_pi (x, MPFR_RNDN);
      mpfr_div_2ui (x, x, 7, MPFR_RNDN);
      mpfr_add_ui (x, x, 1, MPFR_RNDN);
      mpfr_mul (x, x, temp, MPFR_RNDN);

      test = 1;
      while (test)
        {
          mpfr_fprintf (output, "%e\t", mpfr_get_d (x, MPFR_RNDN));
          mpfr_fprintf (output, "%Pu\t", (mpfr_prec_t)prec);

          s.r = (mp_limb_t)mpfr_get_exp (x);
          s.size = (mpfr_prec_t)prec;
          s.align_xp = (mpfr_sgn (x) > 0)?1:2;
          mpfr_set_prec (x2, (mpfr_prec_t)prec);
          mpfr_set (x2, x, MPFR_RNDU);
          s.xp = x2->_mpfr_d;

          for (i=0; i<nb_functions; i++)
            {
              t[i] = speed_measure (param.speed_funcs[i], &s);
              mpfr_fprintf (output, "%e\t", t[i]);
            }
          fprintf (output, "%d\n", 1 + find_best (t, nb_functions));

          if (param.logarithmic_scale_x == 0)
            {
              mpfr_add (x, x, incr_x, MPFR_RNDU);
              if (mpfr_cmp_d (x, param.max_x) > 0)
                test=0;
            }
          else
            {
              if (mpfr_sgn (x) < 0 )
                { /* if x<0, it means that logarithmic_scale_x=-1 */
                  mpfr_div (x, x, incr_x, MPFR_RNDU);
                  mpfr_abs (temp, x, MPFR_RNDD);
                  mpfr_log2 (temp, temp, MPFR_RNDD);
                  if (mpfr_cmp_d (temp, param.min_x) < 0)
                    mpfr_neg (x, x, MPFR_RNDN);
                }
              else
                {
                  mpfr_mul (x, x, incr_x, MPFR_RNDU);
                  mpfr_set (temp, x, MPFR_RNDD);
                  mpfr_log2 (temp, temp, MPFR_RNDD);
                  if (mpfr_cmp_d (temp, param.max_x) > 0)
                    test=0;
                }
            }
        }

      prec = ( (param.logarithmic_scale_prec) ? (prec * incr_prec)
               : (prec + incr_prec) );
      fprintf (output, "\n");
    }

  free (t);
  mpfr_clear (incr_x);
  mpfr_clear (x);
  mpfr_clear (x2);
  mpfr_clear (temp);

  return;
}

#define SPEED_MPFR_FUNC_2D(mean_func)                   \
  do                                                    \
    {                                                   \
      double t;                                         \
      unsigned i;                                       \
      mpfr_t w, x;                                      \
      mp_size_t size;                                   \
                                                        \
      SPEED_RESTRICT_COND (s->size >= MPFR_PREC_MIN);   \
      SPEED_RESTRICT_COND (s->size <= MPFR_PREC_MAX);   \
                                                        \
      size = (s->size-1)/GMP_NUMB_BITS+1;               \
      s->xp[size-1] |= MPFR_LIMB_HIGHBIT;               \
      MPFR_TMP_INIT1 (s->xp, x, s->size);               \
      MPFR_SET_EXP (x, (mp_exp_t)s->r);                 \
      if (s->align_xp == 2) MPFR_SET_NEG (x);           \
                                                        \
      mpfr_init2 (w, s->size);                          \
      speed_starttime ();                               \
      i = s->reps;                                      \
                                                        \
      do                                                \
        mean_func (w, x, MPFR_RNDN);                    \
      while (--i != 0);                                 \
      t = speed_endtime ();                             \
                                                        \
      mpfr_clear (w);                                   \
      return t;                                         \
    }                                                   \
  while (0)

mpfr_prec_t mpfr_exp_2_threshold;
mpfr_prec_t old_threshold = MPFR_EXP_2_THRESHOLD;
#undef  MPFR_EXP_2_THRESHOLD
#define MPFR_EXP_2_THRESHOLD mpfr_exp_2_threshold
#include "exp_2.c"

double
timing_exp1 (struct speed_params *s)
{
  mpfr_exp_2_threshold = s->size+1;
  SPEED_MPFR_FUNC_2D (mpfr_exp_2);
}

double
timing_exp2 (struct speed_params *s)
{
  mpfr_exp_2_threshold = s->size-1;
  SPEED_MPFR_FUNC_2D (mpfr_exp_2);
}

double
timing_exp3 (struct speed_params *s)
{
  SPEED_MPFR_FUNC_2D (mpfr_exp_3);
}

int
main (void)
{
  char filename[256] = "tune.dat";
  FILE *output;
  struct speed_params2D param;
  double (*speed_funcs[4]) (struct speed_params *s);
  speed_funcs[0] = timing_exp1;
  speed_funcs[1] = timing_exp2;
  speed_funcs[2] = timing_exp3;
  speed_funcs[3] = NULL;
  output = fopen (filename, "w");
  if (output == NULL)
    {
      fprintf (stderr, "Can't open file '%s' for writing.\n", filename);
      abort ();
    }

  param.min_x = -13.;
  param.max_x = 27.;
  param.min_prec = 50;
  param.max_prec = 100000;
  param.nb_points_x = 400;
  param.nb_points_prec = 400;
  param.logarithmic_scale_x  = -1;
  param.logarithmic_scale_prec = 1;
  param.speed_funcs = speed_funcs;

  generate_2D_sample (output, param);

  fclose (output);
  mpfr_free_cache ();
  return 0;
}
