/* Example file to test mpfr addition against NTL.
   Usage:
   0) compile this file with NTL
   1) compile tadd.c with -DCHECK_EXTERNAL
   2) ./tadd | egrep -v 'Seed|Inf|NaN' > /tmp/add.log
      (Warning, this produces a large file.)
   3) ./RRTest < /tmp/add.log
*/

#include <NTL/mat_RR.h>

NTL_CLIENT

void
ReadBinary (ZZ &x)
{
  int s = 1;
  ZZ b, y;

  cin >> b;

  if (b < 0)
    {
      s = -1;
      b = -b;
    }

  x = 0;
  y = 1;
  while (b != 0)
    {
      x += (b % 10) * y;
      y *= 2;
      b /= 10;
    }
  if (s < 0)
    x = -x;
}

long
ReadRR (RR &a)
{
  long p;
  ZZ x;
  long e;

  cin >> p;
  ReadBinary (x);
  cin >> e;
  MakeRRPrec (a, x, e, p);
  return p;
}

void
Output (RR a)
{
  cout << a.mantissa() << "*2^(" << a.exponent() << ")" << endl;
}

int main()
{
  RR a, b, c, d;
  long line = 0;
  long p;

  while (!feof(stdin))
    {
      if (++line % 1000 == 0)
        cout << "line " << line << endl;
      ReadRR (b);
      //      ReadRR (c);
      p = ReadRR (a);
      SqrRootPrec (d, b, p);
      if (d != a)
        {
          cerr << "error at line " << line << " for b="; Output(b);
          //          cerr << " c="; Output(c);
          cerr << "expected "; Output(a);
          cerr << "got      "; Output(d);
          cerr << "prec(d)=" << p << endl;
        }
    }
}
