#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

void
#if __STDC__
mpfr_clear (mpfr_ptr m)
#else
mpfr_init (m)
     mpfr_ptr m;
#endif
{
  (*_mp_free_func) (m->_mp_d, ((m->_mp_prec>>3) + 1));
}
