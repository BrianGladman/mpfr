/* mpfr_rec_sqrt -- inverse square root

Copyright 2008 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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

#include <stdio.h>
#include <stdlib.h>

#define MPFR_NEED_LONGLONG_H /* for umul_ppmm */
#include "mpfr-impl.h"

/* number of guard bits in the one-word case */
#define GUARD 1

#define ONE ((mp_limb_t) 1)
#define LIMB_SIZE(x) ((((x)-1)>>MPFR_LOG2_BITS_PER_MP_LIMB) + 1)

#define MPFR_COM_N(x,y,n)				\
  {							\
    mp_size_t i;					\
    for (i = 0; i < n; i++)				\
      *((x)+i) = ~*((y)+i);				\
  }

/* check that X is at distance at most 1 ulp of A^{-1/2},
   where X has p bits (msb aligned).
*/
static void
check (mp_srcptr x, mp_srcptr a, mp_prec_t p)
{
  mp_size_t n = 1 + (p - 1) / GMP_NUMB_BITS;
  mp_prec_t l = n * GMP_NUMB_BITS - p; /* low ignored bits */
  mp_ptr y, t, u;

  /* we should have X - 2^p < A^{-1/2} < X + 2^p */

  y = (mp_ptr) malloc (n * sizeof (mp_limb_t));
  t = (mp_ptr) malloc (2 * n * sizeof (mp_limb_t));
  u = (mp_ptr) malloc ((3 * n + 1) * sizeof (mp_limb_t));  

  /* check the low bits of x are zero */
  if (l)
    {
      if (x[0] & ((ONE << l) - ONE))
	{
	  fprintf (stderr, "Error, low bits from x are not zero: %lu\n",
		   x[0] & ((ONE << l) - ONE));
	  exit (1);
	}
    }

  MPN_COPY (y, x, n);
  if (mpn_sub_1 (y, y, n, ONE << l))
    {
      fprintf (stderr, "Error in check, borrow in X - ulp\n");
      exit (1);
    }
  mpn_mul_n (t, y, y, n);
  mpn_mul (u, t, 2 * n, a, n + 1);
  if (u[3 * n] != 0)
    {
      fprintf (stderr, "(X-ulp)^2*A >= B^(3n)\n");
      exit (1);
    }

  MPN_COPY (y, x, n);
  /* if there is a carry in x + ulp(x), then x + ulp(x) = 1,
     and surely A^{-1/2} < x + ulp(x) */
  if (mpn_add_1 (y, y, n, ONE << l) == 0)
    {
      mpn_mul_n (t, y, y, n);
      mpn_mul (u, t, 2 * n, a, n + 1);
      if (u[3 * n] == 0)
	{
	  fprintf (stderr, "(X+ulp)^2*A < B^(3n)\n");
	  exit (1);
	}
    }

  free (u);
  free (t);
  free (y);
}

/* Put in X a p-bit approximation of 1/sqrt(A),
   where X = {x, n}/B^n, A = {a, n+1}/B^n, and n = ceil(p/GMP_NUMB_BITS),
   where B = 2^GMP_NUMB_BITS.
   Note: x and a are left-aligned.
   The error in the approximate result with respect to the true value
   1/sqrt(A) is bounded by 1 ulp.

   If p is not a multiple of GMP_NUMB_BITS, the extra low bits of the input
   A should be 0. Similarly, the extra low bits of the output X are set to 0.

   Assumptions:
   (1) A should be normalized, i.e., 1 <= A[n] < 4, thus 1/2 <= X < 1.
   (2) p >= 11
   (3) {a, n+1} and {x, n} should not overlap
   (4) GMP_NUMB_BITS >= 12 and is even

   Reference: Modern Computer Algebra, Richard Brent and Paul Zimmermann,
   http://www.loria.fr/~zimmerma/mca/pub226.html
   {a, n+1} and {x, n} should not overlap.
*/
void
mpfr_mpn_rec_sqrt (mp_ptr x, mp_srcptr a, mp_prec_t p)
{
  /* the following T1 and T2 are bipartite tables giving initial
     approximation for the inverse square root, with 13-bit input split in
     5+4+4, and 11-bit output. More precisely, if 2048 <= i < 8192,
     with i = a*2^8 + b*2^4 + c, we use for approximation of
     2048/sqrt(i/2048) the value x = T1[16*(a-8)+b] + T2[16*(a-8)+c].
     The largest error is obtained for i = 3347, where x = 1601,
     and 2048/sqrt(i/2048) = 1602.0168227...
  */
  static mp_limb_t T1[384] = {
2040, 2032, 2024, 2016, 2009, 2001, 1994, 1986, 1979, 1972, /* 10 */
1965, 1958, 1951, 1944, 1937, 1930, 1924, 1917, 1911, 1904, /* 20 */
1898, 1891, 1885, 1879, 1873, 1867, 1861, 1855, 1849, 1843, /* 30 */
1837, 1831, 1826, 1820, 1814, 1809, 1803, 1798, 1792, 1787, /* 40 */
1782, 1777, 1771, 1766, 1761, 1756, 1751, 1746, 1741, 1736, /* 50 */
1731, 1727, 1722, 1717, 1712, 1708, 1703, 1698, 1694, 1689, /* 60 */
1685, 1680, 1676, 1672, 1667, 1663, 1659, 1655, 1650, 1646, /* 70 */
1642, 1638, 1634, 1630, 1626, 1622, 1618, 1614, 1610, 1606, /* 80 */
1602, 1598, 1595, 1591, 1587, 1583, 1580, 1576, 1572, 1569, /* 90 */
1565, 1562, 1558, 1555, 1551, 1548, 1544, 1541, 1537, 1534, /* 100 */
1531, 1527, 1524, 1521, 1517, 1514, 1511, 1508, 1505, 1501, /* 110 */
1498, 1495, 1492, 1489, 1486, 1483, 1480, 1477, 1474, 1471, /* 120 */
1468, 1465, 1462, 1459, 1456, 1453, 1450, 1448, 1445, 1442, /* 130 */
1439, 1436, 1434, 1431, 1428, 1426, 1423, 1420, 1418, 1415, /* 140 */
1412, 1410, 1407, 1404, 1402, 1399, 1397, 1394, 1392, 1389, /* 150 */
1387, 1384, 1382, 1379, 1377, 1374, 1372, 1370, 1367, 1365, /* 160 */
1362, 1360, 1358, 1355, 1353, 1351, 1349, 1346, 1344, 1342, /* 170 */
1339, 1337, 1335, 1333, 1331, 1328, 1326, 1324, 1322, 1320, /* 180 */
1318, 1315, 1313, 1311, 1309, 1307, 1305, 1303, 1301, 1299, /* 190 */
1297, 1295, 1293, 1291, 1289, 1287, 1285, 1283, 1281, 1279, /* 200 */
1277, 1275, 1273, 1271, 1269, 1267, 1265, 1264, 1262, 1260, /* 210 */
1258, 1256, 1254, 1252, 1251, 1249, 1247, 1245, 1243, 1242, /* 220 */
1240, 1238, 1236, 1234, 1233, 1231, 1229, 1228, 1226, 1224, /* 230 */
1222, 1221, 1219, 1217, 1216, 1214, 1212, 1211, 1209, 1207, /* 240 */
1206, 1204, 1202, 1201, 1199, 1198, 1196, 1194, 1193, 1191, /* 250 */
1190, 1188, 1187, 1185, 1183, 1182, 1180, 1179, 1177, 1176, /* 260 */
1174, 1173, 1171, 1170, 1168, 1167, 1165, 1164, 1162, 1161, /* 270 */
1159, 1158, 1157, 1155, 1154, 1152, 1151, 1149, 1148, 1147, /* 280 */
1145, 1144, 1142, 1141, 1140, 1138, 1137, 1136, 1134, 1133, /* 290 */
1131, 1130, 1129, 1127, 1126, 1125, 1123, 1122, 1121, 1119, /* 300 */
1118, 1117, 1116, 1114, 1113, 1112, 1110, 1109, 1108, 1107, /* 310 */
1105, 1104, 1103, 1102, 1100, 1099, 1098, 1097, 1095, 1094, /* 320 */
1093, 1092, 1091, 1089, 1088, 1087, 1086, 1085, 1083, 1082, /* 330 */
1081, 1080, 1079, 1077, 1076, 1075, 1074, 1073, 1072, 1071, /* 340 */
1069, 1068, 1067, 1066, 1065, 1064, 1063, 1062, 1060, 1059, /* 350 */
1058, 1057, 1056, 1055, 1054, 1053, 1052, 1051, 1049, 1048, /* 360 */
1047, 1046, 1045, 1044, 1043, 1042, 1041, 1040, 1039, 1038, /* 370 */
1037, 1036, 1035, 1034, 1033, 1032, 1031, 1030, 1029, 1028, /* 380 */
1027, 1026, 1025, 1024};
  static mp_limb_t T2[384] = {
    8, 7, 7, 6, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 7, 6, 6, 5, /* 20 */
    5, 5, 4, 4, 4, 3, 3, 2, 2, 2, 1, 1, 6, 6, 5, 5, 5, 4, 4, 4, /* 40 */
    3, 3, 3, 2, 2, 2, 1, 1, 5, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, /* 60 */
    2, 1, 1, 1, 5, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, /* 80 */
    4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4, 3, 3, /* 100 */
    3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, /* 120 */
    2, 2, 2, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, /* 140 */
    1, 1, 1, 1, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, /* 160 */
    3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 3, 3, 2, 2, /* 180 */
    2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, /* 200 */
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, /* 220 */
    1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, /* 240 */
    2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, /* 260 */
    2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, /* 280 */
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, /* 300 */
    1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, /* 320 */
    2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, /* 340 */
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, /* 360 */
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, /* 380 */
    0, 0, 0, 0};
  mp_size_t n = LIMB_SIZE(p);

  /* A should be normalized */
  MPFR_ASSERTN((1 <= a[n]) && (a[n] < 4));
  /* we should have enough bits in one limb */
  MPFR_ASSERTN(GMP_NUMB_BITS >= 11);
  /* {a, n+1} and {x, n} should not overlap */
  MPFR_ASSERTN((a + n + 1 <= x) || (x + n <= a));
  MPFR_ASSERTN(p >= 11);

  if (p == 11)
    {
      unsigned long i, ab, ac;
      mp_limb_t t;

      /* take the 13 (or 12) most significant bits of A */
      i = (a[1] << 11) | (a[0] >> (GMP_NUMB_BITS - 11));
      switch (i)
	{
	case 2048:
	  t = 2047; /* otherwise we get t=2048, which gives X=1 */
	  break;
	case 2387: /* table value is 1896, error <= 1.0061 ulp */
	  t = 1897;
	  break;
	case 3347: /* table value is 1601, error <= 1.0169 ulp */
	  t = 1602;
	  break;
	default:
	  ab = i >> 4;
	  ac = (ab & 0x3F0) | (i & 0x0F);
	  t = T1[ab - 0x80] + T2[ac - 0x80];
	}
      x[0] = t << (GMP_NUMB_BITS - p);
    }
  else /* p >= 12 */
    {
      mp_prec_t h = (p < 18) ? 11 : (p >> 1) + 2; /* max(11, ceil((p+3)/2)) */
      mp_ptr r, s, t, u;
      mp_size_t xn = LIMB_SIZE(h);
      mp_size_t rn = LIMB_SIZE(2 * h);
      /* the remaining value of t has at most h+3 bits, but since they are
	 in 2h bits aligned left, the number of limbs is the following */
      mp_size_t th = (h - 3) >> MPFR_LOG2_BITS_PER_MP_LIMB;
      mp_size_t ln = n - xn;
      mp_size_t sn, tn;
      mp_size_t un = xn + rn - th;
      int neg;
      mp_limb_t cy, cu;
      mp_prec_t pl = n * GMP_NUMB_BITS - p; /* low bits from x */
      MPFR_TMP_DECL(marker);
      
      mpfr_mpn_rec_sqrt (x + ln, a + ln, h);
      /* the most h significant bits of X are set, X has ceil(h/GMP_NUMB_BITS)
	 limbs, the low (-h) % GMP_NUMB_BITS bits are zero */
      // gmp_printf ("xh=%Nd\n", x+ln, xn);

      MPFR_TMP_MARK (marker);
      /* first step: square X in r, result is exact */
      u = r = (mp_ptr) MPFR_TMP_ALLOC (un * sizeof (mp_limb_t));
      if (2 * h <= GMP_NUMB_BITS) /* xn=rn=1 */
	{
	  mp_limb_t xx = x[ln] >> (GMP_NUMB_BITS >> 1);
	  r ++;
	  r[0] = xx * xx;
	}
      else if (xn == 1) /* xn=1, rn=2 */
	umul_ppmm(r[1], r[0], x[ln], x[ln]);
      else
	{
	  mpn_mul_n (r, x + ln, x + ln, xn);
	  if (rn < 2 * xn)
	    r ++;
	}
      /* now the 2h most significant bits of {r, rn} contains X^2, r has rn
	 limbs, and the low (-2h) % GMP_NUMB_BITS bits are zero */
      // gmp_printf ("r=%Nd\n", r, rn);

      /* Second step: s <- A * (r^2), and truncate the low p bits,
	 i.e., at weight 2^{-2h} (s is aligned to the low significant bits)
       */
      sn = n + rn;
      s = (mp_ptr) MPFR_TMP_ALLOC ((sn + 1) * sizeof (mp_limb_t));
      if (2 * h + p + 2 <= GMP_NUMB_BITS)
	{
	  /* we should have n=rn=1, thus sn=2 */
	  mp_limb_t ah = (a[1] << p) | (a[0] >> (GMP_NUMB_BITS - p));
	  s[1] = (r[0] >> (GMP_NUMB_BITS - 2 * h)) * ah;
	  /* s[1] should have at most 2h+p+2 bits */
	  s[2] = s[1] >> (2 * h + p);
	  s[1] = s[1] << (GMP_NUMB_BITS - (2 * h + p));
	}
      else if (rn == 1) /* rn=1 implies n=1, since rn*GMP_NUMB_BITS >= 2h,
			   and 2h >= p+3 */
	{
	  /* necessarily p <= GMP_NUMB_BITS-3: we can ignore the two low
	     bits from A */
	  umul_ppmm (s[1], s[0], r[0], a[1] << (GMP_NUMB_BITS - 2)
		     | a[0] >> 2);
	  s[2] = mpn_lshift (s, s, 2, 2);
	}
      else
	{
	  /* we have p <= n * GMP_NUMB_BITS
             2h <= rn * GMP_NUMB_BITS with p+3 <= 2h <= p+4
	     thus n <= rn <= n + 1 */
	  MPFR_ASSERTN(rn <= n + 1);
	  mpn_mul (s, a, n + 1, r, rn);
	  /* s should be near B^sn, thus either s[sn] is 1 and s[sn-1] is
	     near 0, or t[sn]=0 and s[sn-1] is near 111...111.
	     We ignore the bits of s after the first 2h ones.
	  */
	}
      // gmp_printf ("s=%Nd\n", s, n+1+rn);
      /* We ignore the bits of s after the first 2h ones. */
      t = s + n; /* pointer to low limb of t */
      tn = rn;
      /* t has in theory (n+rn) * GMP_NUMB_BITS - floor(p / GMP_NUMB_BITS) bits
	 but the upper h-3 bits of 1-t should be zero */

      /* compute t <- 1 - t, which is B^tn - {t, tn+1} */
      neg = t[tn] != 0;
      if (neg == 0) /* Ax^2 < 1 */
	{
	  MPFR_COM_N (t, t, tn);
	  mpn_add_1 (t, t, tn, 1); /* no carry here */
	}
      /* otherwise we do nothing: we already have t-1 in {t, tn} */

      tn -= th; /* we know at least th = floor((h-3)/GMP_NUMB_LIMBS) of the
		   high limbs of {t, tn} are zero */

      /* tn = rn - th, where rn * GMP_NUMB_BITS >= 2*h and 
	 th * GMP_NUMB_BITS <= h-3, thus tn > 0 */
      MPFR_ASSERTN(tn > 0);

      /* Since |x-a^{-1/2}| <= 1.5*2^{-h}, we have 
	 x = a^{-1/2} + e with |e| <= 1.5*2^{-h}, thus
	 x^2 = 1/a + f with |f| <= 3*2^{-h}*a^{-1/2}+2.25*2^{-2h}
	 ax^2 = 1 + g with |g| <= 3*2^{-h}*a^{1/2} + 2.25*2^{-2h}*a
                               <= 6*2^{-h} + 9*2^{-2h} since a < 4.
	 Since we truncated s at 2^{-2h}, we have:
         |s - 1| <= 6*2^{-h} + 10*2^{-2h} < 2^{3-h} since h >= 3.
	 Thus t should have at most h+3 bits, instead of 2h in theory.
      */

      /* u <- x * t
	 {t, tn} contains at least h+3 bits, and {x, xn} contains h bits,
	 thus tn >= xn */
      MPFR_ASSERTN(tn >= xn);
      if (tn == 1) /* necessarily xn=1 */
	umul_ppmm (u[1], u[0], t[0], x[ln]);
      else
	mpn_mul (u, t, tn, x + ln, xn);

      /* we have already discarded th high limbs of t, thus we only have to
	 consider the upper n - th limbs of u */
      sn = n - th; /* sn cannot be zero, since for n=1 we have th=0, and
		      for n>=2 we have p <= n*GMP_NUMB_BITS,
		      thus h=ceil((p+3)/2) <= (p+4)/2 and 
		      th*GMP_NUMB_BITS <= (h-3) <= p/2 <= n/2*GMP_NUMB_BITS */
      MPFR_ASSERTN(sn > 0);
      un -= sn; /* xn + rn - n */
      u += un;
      
      /* u will be shifted to the right, and after that we want that the low
	 pl bits are zero, thus we want now that the pl+1 bits are zero,
	 thus we round u to nearest at bit pl+1 of u[0] */
      cu = mpn_add_1 (u, u, sn, u[0] & (ONE << pl));
      /* mask bits 0..pl of u[0] */
      if (pl + 1 == GMP_NUMB_BITS)
	u[0] = 0;
      else
	u[0] &= ~((ONE << (pl + 1)) - ONE);

      /* We already have filled {x + ln, xn = n - ln}, and we want to add or
	 subtract cu*B^sn + {u, sn} divided by two at position x.
	 sn = n - th, where th contains <= h-3 bits
	 ln = n - xn, where xn contains >= h bits
	 thus sn > ln.
	 Warning: ln might be zero.
      */
      MPFR_ASSERTN(sn > ln);
      /* we can have sn = ln + 2, for example with GMP_NUMB_BITS=32 and
	 p=62, then h=33, n=2, th=0, xn=2, thus sn=2 and ln=0. */
      MPFR_ASSERTN(sn == ln + 1 || sn == ln + 2);
      /* the high sn-ln limbs of u will overlap the low part of {x+ln,xn} */
      if (ln > 0)
	{
	  mpn_rshift (x, u, ln, 1); /* no carry out */
	  x[ln - 1] |= mpn_rshift (u, u + ln, sn - ln, 1);
	}
      else /* ln = 0 */
	mpn_rshift (u, u, sn, 1); /* no carry out */
      /* incorporate possible carry out from rounding of u */
      u[sn - ln - 1] |= cu << (GMP_NUMB_BITS - 1);
      /* we need to add or subtract the overlapping part {u, sn - ln} */
      if (neg == 0)
	cy = mpn_add (x + ln, x + ln, xn, u, sn - ln);
      else /* negative case */
	{
	  cy = mpn_sub (x + ln, x + ln, xn, u, sn - ln);
	  cy = -mpn_sub_1 (x + sn, x + sn, th, cy); /* n - sn = th */
	  if (ln > 0)
	    {
	      MPFR_COM_N (x, x, ln);
	      cy = mpn_add_1 (x, x, n, ONE) - cy;
	      /* we also need to subtract 1 at x[ln] */
	      cy -= mpn_sub_1 (x + ln, x + ln, xn, ONE); /* n - ln = xn */ 
	    }
	}

#if 0
      /* if cu is non-zero (necessarily 1), add/subtract cu/2 at x[sn] */
      if (MPFR_UNLIKELY(cu != 0))
	{
	  abort();
	  if (neg == 0)
	    cy += mpn_add_1 (x + sn - 1, x + sn - 1, th+1, MPFR_LIMB_HIGHBIT);
	  else
	    cy -= mpn_sub_1 (x + sn - 1, x + sn - 1, th+1, MPFR_LIMB_HIGHBIT);
	}
#endif

      /* cy can be 1 when A=1, i.e., {a, n} = B^n. Setting X to
	 1-ulp(1) satisties the error bound of 1 ulp. */
      if (MPFR_UNLIKELY(cy != 0))
	{
	  cy -= mpn_sub_1 (x, x, n, ONE << pl);
	  MPFR_ASSERTN(cy == 0);
	}

      MPFR_TMP_FREE (marker);
    }
}

#ifdef MAIN
#include "cputime.h"

#define N 1000000

/* test for mpfr_mpn_rec_sqrt */
int
main (int argc, char *argv[])
{
  mp_ptr X, A;
  int i;
  mp_prec_t p = atoi (argv[1]);
  mp_size_t n;
  mp_prec_t pl;
  int st;

  if (argc != 2)
    {
      fprintf (stderr, "Usage: %s precision\n", argv[0]);
      exit (1);
    }

  n = 1 + (p - 1) / GMP_NUMB_BITS;
  pl = n * GMP_NUMB_BITS - p;

  X = malloc (n * sizeof (mp_limb_t));
  A = malloc ((n + 1) * sizeof (mp_limb_t));

#ifdef TIMING
  mpn_random2 (A, n);
  A[n] = 1 + (lrand48 () % 3);
  A[0] &= ~((ONE << pl) - ONE);
  st = cputime ();
  for (i = 0; i < N; i++)
    mpfr_mpn_rec_sqrt (X, A, p);
  st = cputime () - st;
  printf ("time per call: %e ms\n", (double) st / (double) N);
#else
  for (i = 0; i < N; i++)
    {
      mpn_random2 (A, n);
      A[n] = 1 + (lrand48 () % 3);
      A[0] &= ~((ONE << pl) - ONE);
      // gmp_printf ("\nA=%Nd\n", A, n + 1);
      mpfr_mpn_rec_sqrt (X, A, p);
      // gmp_printf ("A=%Nd X=%Nd\n", A, n + 1, X, n);
      check (X, A, p);
    }
#endif
  return 0;
}
#endif
