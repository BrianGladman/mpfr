#include <math.h>

/* generate a random double using the whole range of possible values,
   including denormalized numbers, NaN, infinities, ... */
double drand()
{
  double d; int *i;

  i = (int*) &d;
  i[0] = lrand48();
  i[1] = lrand48();
  if (lrand48()%2) d=-d; /* generates negative numbers */
  return d;
}

/* returns the number of ulp's between a and b */
int ulp(a,b) double a,b;
{
  double eps=1.1102230246251565404e-16; /* 2^(-53) */
  if (a==0.0) {
    if (b==0.0) return 0;
    else if (b<0.0) return 2147483647;
    else return -2147483647;
  }
  b = (a-b)/a;
  if (b>0)
    return (int) floor(b/eps);
  else
    return (int) ceil(b/eps);
}

/* return double m*2^e */
double dbl(m,e) double m; int e;
{
  if (e>=0) while (e-->0) m *= 2.0;
  else while (e++<0) m /= 2.0;
  return m;
}
