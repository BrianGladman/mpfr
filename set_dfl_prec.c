#include "gmp.h"
#include "gmp-impl.h"

/* default is 53 bits */
mp_size_t __gmp_default_fp_bit_precision = 53;

void
#if __STDC__
mpfr_set_default_prec (unsigned long int prec_in_bits)
#else
mpfr_set_default_prec (prec_in_bits)
     unsigned long int prec_in_bits;
#endif
{
  __gmp_default_fp_bit_precision = prec_in_bits;
}
