/* 10^6 additions on a PII-400:
precision       *    mpf_add mpfr_add(RNDZ/RNDN/RNDU)   maple   mupad
53 bits       0.003       0.45       0.64/0.59/0.63
100 bits                  0.52       0.76/0.75/0.77
225                       0.54       1.14
500 bits                  0.72       1.58/1.59/1.65
1000 bits                 1.10       2.34
2017                      1.87       3.31
5025                      4.17       4.87
10017                     7.69       7.52
20017                     15.0       13.4
50017                     57.8       37.8
100017                    124.       105.
*/  

/* #define DEBUG */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"

#define ABS(x) (((x)>0) ? (x) : (-x))

double drand()
{
  double d; long int *i;

  i = (long int*) &d;
  i[0] = lrand48();
  i[1] = lrand48();
  if (lrand48()%2) d=-d; /* generates negative numbers */
  return d;
}

/* checks that x+y gives the same results in double
   and with mpfr with 53 bits of precision */
void check(double x, double y, unsigned int rnd_mode, unsigned int px, 
unsigned int py, unsigned int pz, double res)
{
  double z1,z2; mpfr_t xx,yy,zz;

  /*  printf("x=%1.20e, y=%1.20e, rnd_mode=%u px=%u py=%u pz=%u\n",x,y,rnd_mode,
      px, py, pz); */
  mpfr_init2(xx, px);
  mpfr_init2(yy, py);
  mpfr_init2(zz, pz);
  mpfr_set_d(xx, x, rnd_mode);
  mpfr_set_d(yy, y, rnd_mode);
  mpfr_add(zz, xx, yy, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  z1 = (res==0.0) ? x+y : res;
  z2 = mpfr_get_d(zz);
  /* printf("x+y=%1.20e\n",z2); */
  if (px==53 && py==53 && pz==53) res=1.0;
  if (res!=0.0 && z1!=z2) {
    printf("expected sum is %1.20e, got %1.20e\n",z1,z2);
    printf("mpfr_add failed for x=%1.20e y=%1.20e with rnd_mode=%u\n",x,y,rnd_mode);
    exit(1);
  }
  mpfr_clear(xx); mpfr_clear(yy); mpfr_clear(zz);
}

/* returns the number of ulp's between a and b */
int ulp(a,b) double a,b;
{
  double eps=1.1102230246251565404e-16; /* 2^(-53) */
  b = (a-b)/a; if (b<0) b = -b;
  return (int) floor(b/eps);
}

void check2(x,px,y,py,pz,rnd_mode) double x,y; int px,py,pz,rnd_mode;
{
  mpfr_t xx, yy, zz; double z,z2; int u;

#ifdef DEBUG
printf("x=%1.20e,%d y=%1.20e,%d pz=%d,rnd=%d\n",x,px,y,py,pz,rnd_mode);
#endif
  mpfr_init2(xx,px); mpfr_init2(yy,py); mpfr_init2(zz,pz);
  mpfr_set_d(xx, x, rnd_mode);
  mpfr_set_d(yy, y, rnd_mode);
  mpfr_add(zz, xx, yy, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  z = x+y; z2=mpfr_get_d(zz); u=ulp(z,z2);
#ifdef DEBUG
  printf("x+y=%1.20e,%d (%d ulp) rnd_mode=%d\n",z2,pz,u,rnd_mode);
  mpfr_set_d(zz, z2, rnd_mode); 
  printf("i.e."); mpfr_print_raw(zz); putchar('\n');
#endif
    /* one ulp difference is possible due to composed rounding */
  if (px>=53 && py>=53 && pz>=53 && ABS(u)>1) { 
    printf("x=%1.20e,%d y=%1.20e,%d pz=%d,rnd=%d\n",x,px,y,py,pz,rnd_mode);
    printf("got %1.20e\n",z2);
    printf("result should be %1.20e (diff=%d ulp)\n",z,u);
    mpfr_set_d(zz, z, rnd_mode);
    printf("i.e."); mpfr_print_raw(zz); putchar('\n');
    exit(1); }
  mpfr_clear(xx); mpfr_clear(yy); mpfr_clear(zz);
}

check64()
{
  mpfr_t x, t, u;
  mpfr_init2(x, 65); mpfr_init2(t, 65); mpfr_init2(u, 65);
  mpfr_set_str_raw(x, "0.10011010101000110101010000000011001001001110001011101011111011101E623");
  mpfr_set_str_raw(t, "0.10011010101000110101010000000011001001001110001011101011111011100E623");
  mpfr_sub(u, x, t, GMP_RNDU);
  if (mpfr_get_d(u) != 9.4349060620538533806e167) { /* 2^558 */
    printf("Error (1) in mpfr_sub\n"); exit(1);
  }
  mpfr_init2(x, 64); mpfr_init2(t, 64); mpfr_init2(u, 64);
  mpfr_set_str_raw(x, "0.1000011110101111011110111111000011101011101111101101101100000100E-220");
  mpfr_set_str_raw(t, "0.1000011110101111011110111111000011101011101111101101010011111101E-220");
  mpfr_add(u, x, t, GMP_RNDU);
  if ((MANT(u)[0] & 1) != 1) { 
    printf("error in mpfr_add with rnd_mode=GMP_RNDU\n");
    printf("b=  "); mpfr_print_raw(x); putchar('\n');
    printf("c=  "); mpfr_print_raw(t); putchar('\n');
    printf("b+c="); mpfr_print_raw(u); putchar('\n');
    exit(1);
  }
  mpfr_clear(x); mpfr_clear(t); mpfr_clear(u);
}

main(argc,argv) int argc; char *argv[];
{
  double x,y; int i,prec,rnd_mode,px,py,pz;

  check64();
  check(1.16809465359248765399e+196, 7.92883212101990665259e+196,
	GMP_RNDU, 53, 53, 53, 0.0);
  check(3.14553393112021279444e-67, 3.14553401015952024126e-67,
	GMP_RNDU, 53, 53, 53, 0.0);
  printf("Checking random precisions\n");
  srand(getpid());
  check2(2.26531902208967707071e+168,99,-2.67795218510613988524e+168,67,94,2);
  check2(-1.21510626304662318398e+145,70,1.21367733647758957118e+145,65,61,3);
  check2(2.73028857032080744543e+155,83,-1.16446121423113355603e+163,59,125,1);
  check2(-4.38589520019641698848e+78,155,-1.09923643769309483415e+72,15,159,3);
  check2(-1.49963910666191123860e+265,76,-2.30915090591874527520e-191,8,75,1);
  check2(5.17945380930936917508e+112,119,1.11369077158813567738e+108,15,150,1);
  check2(-2.66910493504493276454e-52,117,1.61188644159592323415e-52,61,68,1);
  check2(-1.87427265794105342764e-57,175,1.76570844587489516446e+190,2,115,1);
  check2(6.85523243386777784171e+107,187,-2.78148588123699111146e+48,87,178,3);
  check2(-1.15706375390780417299e-135,94,-1.07455137477117851576e-129,66,111,2);
  check2(-1.15706375390780417299e-135,94,-1.07455137477117851576e-129,66,111,3);
  check2(-3.31624349995221499866e-22,107,-8.20150212714204839621e+156,79,99,3);
  check2(-1.15706375390780417299e-135,94,-1.07455137477117851576e-129,66,111,3);
  check2(-1.08007920352320089721e+150,63,1.77607317509426332389e+73,64,64,0);
  check2(4.49465557237618783128e+53,108,-2.45103927353799477871e+48,60,105,0);
  check2(3.25471707846623300604e-160,81,-7.93846654265839958715e-274,58,54,0);
  check2(-8.88471912490080158206e+253,79,-7.84488427404526918825e+124,95,53,3);
  check2(-2.18548638152863831959e-125,61,-1.22788940592412363564e+226,71,54,0);
  check2(-7.94156823309993162569e+77,74,-5.26820160805275124882e+80,99,101,3);
  check2(-3.85170653452493859064e+189,62,2.18827389706660314090e+158,94,106,3);
  check2(1.07966151149311101950e+46,88,1.13198076934147457400e+46,67,53,0);
  check2(3.36768223223409657622e+209,55,-9.61624007357265441884e+219,113,53,0);
  check2(-6.47376909368979326475e+159,111,5.11127211034490340501e+159,99,62,3);
  check2(-4.95229483271607845549e+220,110,-6.06992115033276044773e+213,109,55,0);
  check2(-6.47376909368979326475e+159,74,5.11127211034490340501e+159,111,75,2);
  check2(2.26531902208967707070e+168,99,-2.67795218510613988525e+168,67,94,2);
  check2(-2.28886326552077479586e-188,67,3.41419438647157839320e-177,60,110,2);
  check2(-2.66910493504493276454e-52,117,1.61188644159592323415e-52,61,68,1);
  check2(2.90983392714730768886e+50,101,2.31299792168440591870e+50,74,105,1);
  check2(2.72046257722708717791e+243,97,-1.62158447436486437113e+243,83,96,0);
  /* checks random precisions */
  for (i=0;i<1000000;i++) {
#ifdef DEBUG
printf("\nTest i=%d\n",i);
#endif
    px = 53 + (rand() % 64); 
    py = 53 + (rand() % 64); 
    pz = 53 + (rand() % 64); 
    rnd_mode = rand() % 4;
    do { x = drand(); } while (isnan(x));
    do { y = drand(); } while (isnan(y));
    check2(x,px,y,py,pz,rnd_mode);
  }
  printf("Checking double precision (53 bits)\n");
  prec = (argc<2) ? 53 : atoi(argv[1]);
  rnd_mode = (argc<3) ? -1 : atoi(argv[2]);
  check(-8.22183238641455905806e-19, 7.42227178769761587878e-19,
	GMP_RNDD, 53, 53, 53, 0.0);
  check(5.82106394662028628236e+234, -5.21514064202368477230e+89,
	GMP_RNDD, 53, 53, 53, 0.0);
  check(5.72931679569871602371e+122, -5.72886070363264321230e+122,
	GMP_RNDN, 53, 53, 53, 0.0);
  check(-5.09937369394650450820e+238, 2.70203299854862982387e+250,
	GMP_RNDD, 53, 53, 53, 0.0);
  check(-2.96695924472363684394e+27, 1.22842938251111500000e+16,
	GMP_RNDD, 53, 53, 53, 0.0);
  check(1.74693641655743793422e-227, -7.71776956366861843469e-229,
	GMP_RNDN, 53, 53, 53, 0.0);
  check(-1.03432206392780011159e-125, 1.30127034799251347548e-133,
	GMP_RNDN, 53, 53, 53, 0.0);
  check(1.05824655795525779205e+71, -1.06022698059744327881e+71,
	GMP_RNDZ, 53, 53, 53, 0.0);
  check(-5.84204911040921732219e+240, 7.26658169050749590763e+240,
	GMP_RNDD, 53, 53, 53, 0.0);
  /* the following check double overflow */
  /*  check(6.27557402141211962228e+307, 1.32141396570101687757e+308,
      GMP_RNDZ, 53, 53, 53, 0.0); */
  check(1.00944884131046636376e+221, 2.33809162651471520268e+215,
	GMP_RNDN, 53, 53, 53, 0.0);
  check(4.29232078932667367325e-278, 1.07735250473897938332e-281,
	GMP_RNDU, 53, 53, 53, 0.0);
  check(5.27584773801377058681e-80, 8.91207657803547196421e-91,
	GMP_RNDN, 53, 53, 53, 0.0);
  check(2.99280481918991653800e+272, 5.34637717585790933424e+271,
	GMP_RNDN, 53, 53, 53, 0.0);
  check(4.67302514390488041733e-184, 2.18321376145645689945e-190,
	GMP_RNDN, 53, 53, 53, 0.0);
  check(5.57294120336300389254e+71, 2.60596167942024924040e+65,
	GMP_RNDZ, 53, 53, 53, 0.0);
  x=6151626677899716.0; for (i=0;i<30;i++) x = 2.0*x;
  check(x, 4938448004894539.0, GMP_RNDU, 53, 53, 53, 0.0);
  check(1.23056185051606761523e-190, 1.64589756643433857138e-181,
	GMP_RNDU, 53, 53, 53, 0.0);
  check(2.93231171510175981584e-280, 3.26266919161341483877e-273,
	GMP_RNDU, 53, 53, 53, 0.0);
  check(5.76707395945001907217e-58, 4.74752971449827687074e-51,
	GMP_RNDD, 53, 53, 53, 0.0);
  check(277363943109.0, 11.0, GMP_RNDN, 53, 53, 53, 0.0);
  /* test denormalized numbers too */
  check(8.06294740693074521573e-310, 6.95250701071929654575e-310,
	GMP_RNDU, 53, 53, 53, 0.0);
  /* compares to results with double precision using machine arithmetic */
  for (i=0;i<1000000;i++) {
    x = drand(); 
    y = drand();
    if (ABS(x)>2.2e-307 && ABS(y)>2.2e-307 && x+y<1.7e+308 && x+y>-1.7e308) 
      /* avoid denormalized numbers and overflows */
      check(x, y, (rnd_mode==-1) ? lrand48()%4 : rnd_mode, 
	    prec, prec, prec, 0.0);
  } 
}

