#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include <sys/time.h>
#include <math.h>

const double BETA = 0.7;
const mp_size_t SHORT_MUL_THRESHOLD = 1; // 2*KARATSUBA_MUL_THRESHOLD;

#include <sys/resource.h>

int
cputime (void)
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

/*
Low short multiplication using the classical multiplication scheme.
Adds the h least significant limbs of the product of (s1p,size1) and
(s2p,size2) to (rp,h). An eventual carry is ignored.
size1>=size2
h<=size1+size2
*/
mpn_addlowshortmul_clas(rp, s1p, size1, s2p, size2, h)
     mp_ptr rp, s1p, s2p;
     mp_size_t size1, size2, h;
{
  mp_ptr p2=s2p, p3=rp, s2end = s2p+size2, msl, mslp;
  mp_size_t s;
  if (size2>=h) {
    for ( ; h; p2++, p3++, h--)
      mpn_addmul_1(p3, s1p, h, *p2);
    return;
  }
  if (size1>=h) {
    for ( ; p2<s2end; p2++, p3++, h--)
      mpn_addmul_1(p3, s1p, h, *p2);
    return;
  }
  s = h-size1;
  TMP_DECL (marker);
  TMP_MARK (marker);
  msl = (mp_ptr) TMP_ALLOC(s*BYTES_PER_MP_LIMB);
  for (mslp=msl; s; p2++, p3++, s--, mslp++)
    *mslp = mpn_addmul_1(p3, s1p, size1, *p2);
  mpn_add_n(rp+size1, rp+size1, msl, h-size1);
  TMP_FREE (marker);
  for ( ; p2<s2end; p2++, p3++, size1--) 
    mpn_addmul_1(p3, s1p, size1, *p2);
}

/*
Low short multiplication using fast multiplication.
Adds the size least significant limbs of the product of (s1p,size) and
(s2p,size) in (rp,size). An eventual carry is ignored.
*/
mpn_addlowshortmul_fast_1(rp, s1p, s2p, size)
     mp_ptr rp, s1p, s2p;
     mp_size_t size;
{
  mp_ptr fp;
  if (size<SHORT_MUL_THRESHOLD)
    mpn_addlowshortmul_clas(rp, s1p, size, s2p, size, size);
  else {
    mp_size_t newsize = ceil(BETA*size);
    TMP_DECL (marker);
    TMP_MARK (marker);
    fp = (mp_ptr) TMP_ALLOC(2*newsize*BYTES_PER_MP_LIMB);
    mpn_mul(fp, s1p, newsize, s2p, newsize);
    mpn_add(rp, rp, size, fp, size);
    TMP_FREE (marker);

    mpn_addlowshortmul_fast_1(rp+newsize, s1p+newsize, s2p, size-newsize);
    mpn_addlowshortmul_fast_1(rp+newsize, s1p, s2p+newsize, size-newsize);
  }
}

mpn_addlowshortmul_fast_2(rp, s1p, s2p, size, newsize)
     mp_ptr rp, s1p, s2p;
     mp_size_t size;
     mp_size_t newsize;
{
  mp_ptr fp;
  int dif = size-newsize;
  
  TMP_DECL (marker);
  TMP_MARK (marker);
  fp = (mp_ptr) TMP_ALLOC(2*newsize*BYTES_PER_MP_LIMB);
  mpn_mul(fp, s1p, newsize, s2p, newsize);
  mpn_add(rp, rp, size, fp, size);
  TMP_FREE (marker);

  mpn_addlowshortmul_clas(rp+newsize, s1p+newsize, dif, s2p, dif, dif);
  // mpn_addlowshortmul_fast_1(rp+newsize, s1p+newsize, s2p, size-newsize);
  mpn_addlowshortmul_clas(rp+newsize, s1p, dif, s2p+newsize, dif, dif);
  // mpn_addlowshortmul_fast_1(rp+newsize, s1p, s2p+newsize, size-newsize);
}

/*
Low short multiplication using fast multiplication.
Adds the h least significant limbs of the product of (s1p,size1) and
(s2p,size2) to (rp,h). An eventual carry is ignored.
size1>=size2
h<=size1+size2
*/
mpn_addlowshortmul_fast(rp, s1p, size1, s2p, size2, h)
     mp_ptr rp, s1p, s2p;
     mp_size_t size1, size2, h;
{
  mp_size_t ds;
  mp_ptr prod;
  if (h<=size2) {
    mpn_addlowshortmul_fast_1(rp, s1p, s2p, h);
    return;
  }
  if (h<=size1) {
    ds = h-size2;
    TMP_DECL (marker);
    TMP_MARK (marker);
    prod = (mp_ptr) TMP_ALLOC(h*BYTES_PER_MP_LIMB);
    if (ds>=size2)
      mpn_mul(prod, s1p, ds, s2p, size2);
    else
      mpn_mul(prod, s2p, size2, s1p, ds);
    mpn_addlowshortmul_fast_1(prod+ds, s1p+ds, s2p, size2);
    mpn_add_n(rp, rp, prod, h);
    TMP_FREE (marker);
    return;
  }
  ds = h-size2;
  TMP_DECL (marker);
  TMP_MARK (marker);
  prod = (mp_ptr) TMP_ALLOC(h*BYTES_PER_MP_LIMB);
  if (ds>=size2)
    mpn_mul(prod, s1p, ds, s2p, size2);
  else
    mpn_mul(prod, s2p, size2, s1p, ds);
  mpn_addlowshortmul_fast(prod+ds, s2p, size2, s1p+ds, size1-ds, size2);
  mpn_add_n(rp, rp, prod, h);
  TMP_FREE (marker);
}

/*
Low short multiplication.
Returns the h least significant limbs of the product of (s1p,size1) and
(s2p,size2). sp has size at least h.
*/
void mpn_lowshortmul(sp, s1p, size1, s2p, size2, h)
     mp_ptr sp;
     mp_srcptr s1p, s2p;
     mp_size_t size1, size2, h;
{
  if (h>=size1+size2)
    if (size1>=size2)
      mpn_mul(sp, s1p, size1, s2p, size2);
    else
      mpn_mul(sp, s2p, size2, s1p, size1);
  else {
    MPN_ZERO(sp, h);
    if (size1>=size2)
      mpn_addlowshortmul_fast(sp, s1p, size1, s2p, size2, h);
    else
      mpn_addlowshortmul_fast(sp, s2p, size2, s1p, size1, h);
  }
}

/*
High short multiplication using the classical multiplication scheme.
Let hsp be the high short product of (s1p,size1) and (s2p,size2) with
bound h, i.e. hsp is the sum of all s1p[i]*s2p[j]*beta^(i+j) with i+j>=h.
Adds hsp[h]+hsp[h+1]*beta+hsp[h+2]*beta^2+... to (rp,size1+size2-h).
Returns carry, i.e. 0 or 1.
size1>=size2
h<size1+size2
*/
mp_limb_t mpn_addhighshortmul_clas(rp, s1p, size1, s2p, size2, h)
     mp_ptr rp, s1p, s2p;
     mp_size_t size1, size2, h;
{
  mp_ptr p1, p2, pr, s2end = s2p+size2, msl, mslp;
  mp_limb_t retval;
  mp_size_t s;
  if (h>=size1-1) {
    TMP_DECL (marker);
    TMP_MARK (marker);
    msl = (mp_ptr) TMP_ALLOC((size2-(h-(size1-1)))*BYTES_PER_MP_LIMB);
    mslp = msl;
    
    for (p1=s1p+size1-1, p2=s2p+(h-(size1-1)), s=1; p2<s2end; p1--, p2++, s++, mslp++) {
      *mslp = mpn_addmul_1(rp, p1, s, *p2);
    }
    retval =  mpn_add_n(rp+1, rp+1, msl, size2-(h-(size1-1)));
    TMP_FREE (marker);
    return retval;
  }
  if (h>=size2-1) {
    TMP_DECL (marker);
    TMP_MARK (marker);
    msl = (mp_ptr) TMP_ALLOC(size2*BYTES_PER_MP_LIMB);

    for (mslp=msl, p1=s1p+h, p2=s2p, s=size1-h; p2<s2end; p1--, p2++, s++, mslp++)
      *mslp = mpn_addmul_1(rp, p1, s, *p2);
    
    retval = mpn_add_n(rp+size1-h, rp+size1-h, msl, size2);
    TMP_FREE (marker);
    return retval;
  }

  TMP_DECL (marker);
  TMP_MARK (marker);
  msl = (mp_ptr) TMP_ALLOC(size2*BYTES_PER_MP_LIMB);
  for (mslp=msl, p1=s1p+h, p2=s2p, s=size1-h; p1>=s1p; p1--, p2++, mslp++, s++)
      *mslp = mpn_addmul_1(rp, p1, s, *p2);
  for (pr=rp+1; p2<s2end; p2++, mslp++, pr++) 
    *mslp = mpn_addmul_1(pr, s1p, size1, *p2);
  retval = mpn_add_n(rp+size1-h, rp+size1-h, msl, size2);
  TMP_FREE (marker);
  return retval;
}

/*
High short multiplication using fast multiplication.
Let p be the product of (s1p,size) and (s2p,size).
Let hsp be the high short product of (s1p,size) and (s2p,size), i.e.
hsp is the sum of all s1p[i]*s2p[j]*beta^(i+j) with i+j>=size.
Adds N[size]+N[size+1]*beta+N[size+2]*beta^2+... to (rp,size), where
hsp <= N <= p.
Returns carry, i.e. 0 or 1.
Here, the formulation with N is needed because you have no control
over the carries in the recursive calls.
*/
mp_limb_t mpn_addhighshortmul_fast_1(rp, s1p, s2p, size)
     mp_ptr rp, s1p, s2p;
     mp_size_t size;
{
  mp_ptr fp;
  if (size<SHORT_MUL_THRESHOLD)
    return mpn_addhighshortmul_clas(rp, s1p, size, s2p, size, size);
  else {
    mp_limb_t msl1, msl2;
    
    mp_size_t newsize = ceil(BETA*size);
    TMP_DECL (marker);
    TMP_MARK (marker);
    fp = (mp_ptr) TMP_ALLOC(2*newsize*BYTES_PER_MP_LIMB);
    mpn_mul(fp, s1p+size-newsize, newsize, s2p+size-newsize, newsize);
    msl1 = mpn_add_n(rp, rp, fp+2*newsize-size, size);
    TMP_FREE (marker);

    msl2 = mpn_addhighshortmul_fast_1(rp, s1p+newsize, s2p, size-newsize);
    msl1 += mpn_add_1(rp+size-newsize, rp+size-newsize, newsize, msl2); 
    msl2 = mpn_addhighshortmul_fast_1(rp, s1p, s2p+newsize, size-newsize);
    msl1 += mpn_add_1(rp+size-newsize, rp+size-newsize, newsize, msl2);
    return msl1;
  }
}

/*
High short multiplication using fast multiplication.
Let p be the product of (s1p,size) and (s2p,size).
Let hsp be the high short product of (s1p,size1) and (s2p,size2) with
bound h, i.e. hsp is the sum of all s1p[i]*s2p[j]*beta^(i+j) with i+j>=h.
Adds N[h]+N[h+1]*beta+N[h+2]*beta^2+... to (rp,size1+size2-h), where
hsp <= N <= p.
Returns carry, i.e. 0 or 1.
size1>=size2
h<size1+size2
Here, the formulation with N is needed because you have no control
over the carries in the recursive calls.
*/
mp_limb_t mpn_addhighshortmul_fast(rp, s1p, size1, s2p, size2, h)
     mp_ptr rp, s1p, s2p;
     mp_size_t size1, size2, h;
{
  mp_size_t ds, s;
  mp_ptr prod;
  mp_limb_t msl, retval;
  if (h>=size1) {
    s = size1+size2-h;
    return mpn_addhighshortmul_fast_1(rp, s1p+size1-s, s2p+size2-s, s);
  }
  if (h>=size2) {
    ds = size1-h;
    TMP_DECL (marker);
    TMP_MARK (marker);
    prod = (mp_ptr) TMP_ALLOC((ds+size2)*BYTES_PER_MP_LIMB);
    if (ds>=size2)
      mpn_mul(prod, s1p+h, ds, s2p, size2);
    else
      mpn_mul(prod, s2p, size2, s1p+h, ds);
    msl = mpn_addhighshortmul_fast_1(prod, s1p+(h-size2), s2p, size2);
    mpn_add_1(prod+size2, prod+size2, ds, msl);
    retval = mpn_add_n(rp, rp, prod, size1+size2-h);
    TMP_FREE (marker);
    return retval;
  }
  ds = size1-h;
  TMP_DECL (marker);
  TMP_MARK (marker);
  prod = (mp_ptr) TMP_ALLOC((ds+size2)*BYTES_PER_MP_LIMB);
  if (ds>=size2)
    mpn_mul(prod, s1p+h, ds, s2p, size2);
  else
    mpn_mul(prod, s2p, size2, s1p+h, ds);
  msl = mpn_addhighshortmul_fast(prod, s2p, size2, s1p, h, h);
  mpn_add_1(prod+size2, prod+size2, ds, msl);
  retval = mpn_add_n(rp, rp, prod, size1+size2-h);
  TMP_FREE (marker);
  return retval;
}

/*
High short multiplication.
Let p be the product of (s1p,size1) and (s2p,size2).
Let hsp be the high short product of (s1p,size1) and (s2p,size2) with
bound h, i.e. hsp is the sum of all s1p[i]*s2p[j]*beta^(i+j) with i+j>=h.
Returns N[h]+N[h+1]*beta+N[h+2]*beta^2+..., where hsp <= N <= p.
The result contains size1+size2-h limbs.
If h>=size1+size2, the result contains one zero limb.
Here, the formulation with N is needed because you have no control
over the carries in the recursive calls.
*/
mp_ptr mpn_highshortmul_approx(s1p, size1, s2p, size2, h)
     mp_ptr s1p, s2p;
     mp_size_t size1, size2, h;
{
  mp_ptr sp;
  if (h>=size1+size2) {
    sp = alloca(BYTES_PER_MP_LIMB);
    MPN_ZERO(sp,1);
  }
  sp = alloca((size1+size2-h)*BYTES_PER_MP_LIMB);
  MPN_ZERO(sp,size1+size2-h);
  if (size1>=size2)
    mpn_addhighshortmul_fast(sp, s1p, size1, s2p, size2, h);
  else
    mpn_addhighshortmul_fast(sp, s2p, size2, s1p, size1, h);
  return sp;
}

/*
High short multiplication.
Let p be the product of (s1p,size1) and (s2p,size2).
Returns p[h]+p[h+1]*beta+p[h+2]*beta^2+....
The result contains size1+size2-h limbs.
If h>=size1+size2, the result contains one zero limb.
*/
mp_ptr mpn_highshortmul(s1p, size1, s2p, size2, h)
     mp_ptr s1p, s2p;
     mp_size_t size1, size2, h;
{
  mp_ptr sp, prod;
  mp_size_t s;
  if (h>=size1+size2) {
    sp = (mp_ptr) calloc(1, sizeof(mp_limb_t));
    return sp;
  }
  if (h<2) {
    sp = (mp_ptr) calloc(size1+size2, sizeof(mp_limb_t));
    if (size1>=size2)
      mpn_mul(sp, s1p, size1, s2p, size2);
    else
      mpn_mul(sp, s2p, size2, s1p, size1);
    return sp+h;
  }
  sp = (mp_ptr) calloc(size1+size2-h+2, sizeof(mp_limb_t));
  s = MAX(size1, size2);
  if (size1>=size2)
    mpn_addhighshortmul_fast(sp, s1p, size1, s2p, size2, h-2);
  else
    mpn_addhighshortmul_fast(sp, s2p, size2, s1p, size1, h-2);
  if (sp[1]+s<sp[1]) {
    // We cannot guarantee correctness, therefore perform
    // a full multiplication 
    sp = (mp_ptr) calloc(size1+size2, sizeof(mp_limb_t));
    if (size1>=size2)
      mpn_mul(sp, s1p, size1, s2p, size2);
    else
      mpn_mul(sp, s2p, size2, s1p, size1);
    return sp+h;
  }
  return sp+2;
}






// TESTING
test_mpn_addlowshortmul_clas(size1, size2, h)
     mp_size_t size1, size2, h;
{
  mp_ptr s1p, s2p, rp1, rp2, prod;
  mp_size_t i;

  printf("test_mpn_addlowshortmul_clas: "); fflush(stdout);

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  rp1 = alloca(h*BYTES_PER_MP_LIMB);
  mpn_random(rp1, h);
  rp2 = alloca(h*BYTES_PER_MP_LIMB);
  for (i=0; i<h; i++)
    *(rp2+i) = *(rp1+i);

  mpn_addlowshortmul_clas(rp1, s1p, size1, s2p, size2, h);

  prod = alloca((size1+size2)*BYTES_PER_MP_LIMB);
  mpn_mul(prod, s1p, size1, s2p, size2);
  mpn_add(rp2, rp2, h, prod, h);

  for (i=0; i<h; i++)
    if (*(rp1+i)-*(rp2+i)) {
      printf("ERROR\n");
      return;
      printf("%d  %u  %u\n", i, *(rp1+i), *(rp2+i));
    }

  printf("OK\n"); fflush(stdout);
}

test_mpn_addlowshortmul_fast(size1, size2, h)
     mp_size_t size1, size2, h;
{
  mp_ptr s1p, s2p, rp1, rp2, prod;
  mp_size_t i;

  printf("test_mpn_addlowshortmul_fast: ");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  rp1 = alloca(h*BYTES_PER_MP_LIMB);
  mpn_random(rp1, h);
  rp2 = alloca(h*BYTES_PER_MP_LIMB);
  for (i=0; i<h; i++)
    *(rp2+i) = *(rp1+i);

  mpn_addlowshortmul_fast(rp1, s1p, size1, s2p, size2, h);

  prod = alloca((size1+size2)*BYTES_PER_MP_LIMB);
  mpn_mul(prod, s1p, size1, s2p, size2);
  mpn_add(rp2, rp2, h, prod, h);

  for (i=0; i<h; i++)
    if (*(rp1+i)-*(rp2+i)) {
      //printf("ERROR\n");
      //return;
      printf("%d  %u  %u\n", i, *(rp1+i), *(rp2+i));
    }

  printf("OK\n");
}

test_mpn_addhighshortmul_clas(size1, size2, h)
     mp_size_t size1, size2, h;
{
  mp_ptr s1p, s2p, rp1, rp2, prod, p1, p2;
  mp_limb_t msl, carry; 
  mp_size_t i;

  printf("test_mpn_addhighshortmul_clas: ");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  rp1 = alloca((size1+size2-h)*BYTES_PER_MP_LIMB);
  MPN_ZERO(rp1, size1+size2-h);
  rp2 = alloca((size1+size2-h)*BYTES_PER_MP_LIMB);
  MPN_ZERO(rp2, size1+size2-h);

  mpn_addhighshortmul_clas(rp1, s1p, size1, s2p, size2, h);

  prod = alloca((size1+size2)*BYTES_PER_MP_LIMB);
  mpn_mul(prod, s1p, size1, s2p, size2);
  mpn_add_n(rp2, rp2, prod+h, size1+size2-h);

  for (i=2; i<size1+size2-h; i++)
    if (*(rp1+i)-*(rp2+i)) {
      printf("ERROR\n");
      return;
      printf("%d  %u  %u\n", i, *(rp1+i), *(rp2+i));
    }
  
  printf("OK\n");

  //  printf("%d  %u  %u\n", 0, rp1[0], rp2[0]);
  //  printf("%d  %u  %u\n", 1, rp1[1], rp2[1]);
}

test_mpn_addhighshortmul_fast(size1, size2, h)
     mp_size_t size1, size2, h;
{
  mp_ptr s1p, s2p, rp1, rp2, prod;
  mp_size_t i;

  printf("test_mpn_addhighshortmul_fast: ");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  rp1 = alloca((size1+size2-h)*BYTES_PER_MP_LIMB);
  MPN_ZERO(rp1, size1+size2-h);
  rp2 = alloca((size1+size2-h)*BYTES_PER_MP_LIMB);
  MPN_ZERO(rp2, size1+size2-h);

  mpn_addhighshortmul_fast(rp1, s1p, size1, s2p, size2, h);

  prod = alloca((size1+size2)*BYTES_PER_MP_LIMB);
  mpn_mul(prod, s1p, size1, s2p, size2);
  mpn_add_n(rp2, rp2, prod+h, size1+size2-h);

  for (i=2; i<size1+size2-h; i++)
    if (*(rp1+i)-*(rp2+i)) {
      printf("ERROR\n");
      return;
      printf("%d  %u  %u\n", i, *(rp1+i), *(rp2+i));
    }
  
  printf("OK\n");
  
  //  printf("%d  %u  %u\n", 0, rp1[0], rp2[0]);
  //  printf("%d  %u  %u\n", 1, rp1[1], rp2[1]);
}

test_mpn_lowshortmul(size1, size2, h)
     mp_size_t size1, size2, h;
{
  mp_ptr s1p, s2p, sp, prod;
  mp_size_t s, i;

  printf("test_mpn_lowshortmul: ");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  sp = alloca(h*BYTES_PER_MP_LIMB);

  mpn_lowshortmul(sp, s1p, size1, s2p, size2, h);

  s = MAX(h, size1+size2);
  prod = alloca(s*BYTES_PER_MP_LIMB);
  MPN_ZERO(prod,s);
  mpn_mul(prod, s1p, size1, s2p, size2);

  for (i=0; i<h; i++)
    if (sp[i]-prod[i]) {
      printf("ERROR\n");
      return;
      printf("%d  %u  %u\n", i, sp[i], prod[i]);
    }
  
  printf("OK\n");
}

test_mpn_highshortmul(size1, size2, h)
     mp_size_t size1, size2, h;
{
  mp_ptr s1p, s2p, sp, prod;
  mp_size_t i;

  printf("test_mpn_highshortmul: ");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);

  sp = mpn_highshortmul(s1p, size1, s2p, size2, h);

  prod = alloca((size1+size2)*BYTES_PER_MP_LIMB);
  mpn_mul(prod, s1p, size1, s2p, size2);

  for (i=0; i<size1+size2-h; i++)
    if (*(sp+i)-*(prod+h+i)) {
      printf("ERROR\n");
      return;
      printf("%d  %u  %u\n", i, *(sp+i), *(prod+i));
    }
  
  printf("OK\n");
}



// TIMING
void start_timer(struct timeval *timer)
{
  gettimeofday(timer, NULL);
}

struct timeval elapsed(struct timeval timer)
{
  struct timeval now, elapse;
  gettimeofday(&now, NULL);
  if (now.tv_usec-timer.tv_usec >=0) {
    elapse.tv_sec = now.tv_sec-timer.tv_sec;
    elapse.tv_usec = now.tv_usec-timer.tv_usec;
  }
  else {
    elapse.tv_sec = now.tv_sec-timer.tv_sec-1;
    elapse.tv_usec= 1000000+now.tv_usec-timer.tv_usec;
  }
  return elapse;
}

time_mpn_addlowshortmul_clas(size1, size2, h, c)
     mp_size_t size1, size2, h;
     int c;
{
  mp_ptr s1p, s2p, rp, prod;
  struct timeval timer, elapse;
  int i;

  printf("time_mpn_addlowshortmul_clas:\n");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  rp = alloca(h*BYTES_PER_MP_LIMB);
  mpn_random(rp, h);

  start_timer(&timer);
  for (i=c ; i; i--)
    mpn_addlowshortmul_clas(rp, s1p, size1, s2p, size2, h);
  elapse = elapsed(timer);
  printf("  addlowshortmul_clas: %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);

  start_timer(&timer);
  for (i=c ; i; i--) {
    TMP_DECL (marker);
    TMP_MARK (marker);
    prod = (mp_ptr) TMP_ALLOC((size1+size2)*BYTES_PER_MP_LIMB);
    mpn_mul(prod, s1p, size1, s2p, size2);
    TMP_FREE (marker);
  }
  elapse = elapsed(timer);
  printf("  addfullmul:          %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);
}

time_mpn_addlowshortmul_fast(size1, size2, h, c)
     mp_size_t size1, size2, h;
     int c;
{
  mp_ptr s1p, s2p, rp, prod;
  struct timeval timer, elapse;
  int i;

  printf("time_mpn_addlowshortmul_fast:\n");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  rp = alloca(h*BYTES_PER_MP_LIMB);
  mpn_random(rp, h);

  start_timer(&timer);
  for (i=c ; i; i--)
    mpn_addlowshortmul_fast(rp, s1p, size1, s2p, size2, h);
  elapse = elapsed(timer);
  printf("  addlowshortmul_fast: %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);

  start_timer(&timer);
  for (i=c ; i; i--) {
    TMP_DECL (marker);
    TMP_MARK (marker);
    prod = (mp_ptr) TMP_ALLOC((size1+size2)*BYTES_PER_MP_LIMB);
    mpn_mul(prod, s1p, size1, s2p, size2);
    TMP_FREE (marker);
  }
  elapse = elapsed(timer);
  printf("  addfullmul:          %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);
}

time_mpn_addhighshortmul_clas(size1, size2, h, c)
     mp_size_t h;
     int c;
{
  mp_ptr s1p, s2p, rp, prod;
  struct timeval timer, elapse;
  int i;

  printf("time_mpn_addhighshortmul_clas:\n");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  rp = alloca((size1+size2-h)*BYTES_PER_MP_LIMB);
  MPN_ZERO(rp,size1+size2-h);

  start_timer(&timer);
  for (i=c ; i; i--)
    mpn_addhighshortmul_clas(rp, s1p, size1, s2p, size2, h);
  elapse = elapsed(timer);
  printf("  addhighshortmul_clas: %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);

  start_timer(&timer);
  for (i=c ; i; i--) {
    TMP_DECL (marker);
    TMP_MARK (marker);
    prod = (mp_ptr) TMP_ALLOC((size1+size2)*BYTES_PER_MP_LIMB);
    mpn_mul(prod, s1p, size1, s2p, size2);
    TMP_FREE (marker);
  }
  elapse = elapsed(timer);
  printf("  addfullmul:           %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);
}

time_mpn_addhighshortmul_fast(size1, size2, h, c)
     mp_size_t h;
     int c;
{
  mp_ptr s1p, s2p, rp, prod;
  struct timeval timer, elapse;
  int i;

  printf("time_mpn_addhighshortmul_fast:\n");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);
  rp = alloca((size1+size2-h)*BYTES_PER_MP_LIMB);
  MPN_ZERO(rp,size1+size2-h);

  start_timer(&timer);
  for (i=c ; i; i--)
    mpn_addhighshortmul_fast(rp, s1p, size1, s2p, size2, h);
  elapse = elapsed(timer);
  printf("  addhighshortmul_fast: %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);

  start_timer(&timer);
  for (i=c ; i; i--) {
    TMP_DECL (marker);
    TMP_MARK (marker);
    prod = (mp_ptr) TMP_ALLOC((size1+size2)*BYTES_PER_MP_LIMB);
    mpn_mul(prod, s1p, size1, s2p, size2);
    TMP_FREE (marker);
  }
  elapse = elapsed(timer);
  printf("  addfullmul:           %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);
}

time_mpn_lowshortmul(size1, size2, h, c)
     mp_size_t size1, size2, h;
     int c;
{
  mp_ptr s1p, s2p, sp, prod;
  int i, st1, st2;

  printf("time_mpn_lowshortmul:\n");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);

  TMP_DECL (marker);
  TMP_MARK (marker);
  prod = (mp_ptr) TMP_ALLOC((size1+size2)*BYTES_PER_MP_LIMB);
  st1 = cputime();
  for (i=c ; i; i--) mpn_lowshortmul(prod, s1p, size1, s2p, size2, h);
  st1 = cputime()-st1;
  printf("  lowshortmul: %dms\n", st1);

  st2 = cputime();
  for (i=c ; i; i--) mpn_mul(prod, s1p, size1, s2p, size2);
  st2 = cputime()-st2;
  TMP_FREE (marker);
  printf("  fullmul:     %dms\n", st2);
  printf("lowshortmul/fullmul=%1.2f\n", (double) st1 / (double) st2);
}

time_mpn_lowshortmul_all(size1, size2, h, c)
     mp_size_t size1, size2, h;
     int c;
{
  mp_ptr s1p, s2p, sp, prod;
  int i, st1, st2, st3, newsize, opt_newsize;
  
  if (size1<h || size2<h) {
    fprintf(stderr, "Error: size1<h || size2<h\n"); exit(1);
  }

  printf("time_mpn_lowshortmul:\n");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);

  TMP_DECL (marker);
  TMP_MARK (marker);
  prod = (mp_ptr) TMP_ALLOC((size1+size2)*BYTES_PER_MP_LIMB);
  for (newsize=(h+1)/2;newsize<h;newsize++) {
    st1 = cputime();
    for (i=c ; i; i--) mpn_addlowshortmul_fast_2(prod, s1p, s2p, h, newsize);
    st1 = cputime()-st1;
    if (newsize==(h+1)/2) st2=st1+1;
    printf("  lowshortmul(%u): %dms\n", newsize, st1);
    if (st1<st2) { opt_newsize=newsize; st2=st1; }
  }
  printf("best newsize for h=%u is %u\n", h, opt_newsize);
  st1=st2;

  st3 = cputime();
  for (i=c ; i; i--) mpn_addlowshortmul_clas(prod, s1p, size1, s2p, size2, h);
  st3 = cputime()-st3;
  printf("  mpn_addlowshortmul_clas:     %dms\n", st3);

  st2 = cputime();
  for (i=c ; i; i--) mpn_mul(prod, s1p, size1, s2p, size2);
  st2 = cputime()-st2;
  TMP_FREE (marker);
  printf("fullmul:     %dms\n", st2);

  printf("shortclass/fullmul=%1.2f shortfast/fullmul=%1.2f\n",
	 (double) st3 / (double) st2, (double) st1 / (double) st2);
}

time_mpn_highshortmul(size1, size2, h, c)
     mp_size_t size1, size2, h;
     int c;
{
  mp_ptr s1p, s2p, sp, prod;
  struct timeval timer, elapse;
  int i;

  printf("time_mpn_highshortmul:\n");

  s1p = alloca(size1*BYTES_PER_MP_LIMB);
  mpn_random(s1p, size1);
  s2p = alloca(size2*BYTES_PER_MP_LIMB);
  mpn_random(s2p, size2);

  start_timer(&timer);
  for (i=c ; i; i--) {
    sp = mpn_highshortmul(s1p, size1, s2p, size2, h);
    free(sp);
  }
  elapse = elapsed(timer);
  printf("  highshortmul: %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);

  start_timer(&timer);
  for (i=c ; i; i--) {
    TMP_DECL (marker);
    TMP_MARK (marker);
    prod = (mp_ptr) TMP_ALLOC((size1+size2)*BYTES_PER_MP_LIMB);
    mpn_mul(prod, s1p, size1, s2p, size2);
    TMP_FREE (marker);
  }
  elapse = elapsed(timer);
  printf("  fullmul:      %d.%.6d\n",elapse.tv_sec,elapse.tv_usec);
}






main(int argc, char *argv[])
{
  mp_size_t size1, size2, h;
  int c;

  printf("SHORT_MUL_THRESHOLD=%u\n", SHORT_MUL_THRESHOLD);
#if 0
  scanf("%d", &size1);
  scanf("%d", &size2);
  scanf("%d", &h);
  scanf("%d", &c);
#else
  size1 = atoi(argv[1]);
  size2 = atoi(argv[2]);
  h = atoi(argv[3]);
  c = (argc<5) ? 1 : atoi(argv[4]);
#endif

  //  test_mpn_addlowshortmul_clas(size1, size2, h);
  //  time_mpn_addlowshortmul_clas(size1, size2, h, c);

  //  test_mpn_addhighshortmul_clas(size1, size2, h);
  //  time_mpn_addhighshortmul_clas(size1, size2, size1+size2-h, c);

  // test_mpn_addlowshortmul_fast(size1, size2, h);
  // time_mpn_addlowshortmul_fast(size1, size2, h, c);

  //  test_mpn_addhighshortmul_fast(size1, size2, h);
  //  time_mpn_addhighshortmul_fast(size1, size2, size1+size2-h, c);

  test_mpn_lowshortmul(size1, size2, h);
  time_mpn_lowshortmul_all(size1, size2, h, c);

  //  test_mpn_highshortmul(size1, size2, h);
  //  time_mpn_highshortmul(size1, size2, h, c);
}
