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

#define LIMB_SIZE(x) ((((x)-1)>>MPFR_LOG2_BITS_PER_MP_LIMB) + 1)

#define MPFR_COM_N(x,y,n)                               \
  {                                                     \
    mp_size_t i;                                        \
    for (i = 0; i < n; i++)                             \
      *((x)+i) = ~*((y)+i);                             \
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
      if (x[0] & MPFR_LIMB_MASK(l))
        {
          fprintf (stderr, "Error, low bits from x are not zero: %lu\n",
                   x[0] & MPFR_LIMB_MASK(l));
          exit (1);
        }
    }

  MPN_COPY (y, x, n);
  if (mpn_sub_1 (y, y, n, MPFR_LIMB_ONE << l))
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
  if (mpn_add_1 (y, y, n, MPFR_LIMB_ONE << l) == 0)
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
     The largest error is obtained for i = 2054, where x = 2044,
     and 2048/sqrt(i/2048) = 2045.006576...
  */
  static short int T1[384] = {
2040, 2033, 2025, 2017, 2009, 2002, 1994, 1987, 1980, 1972, 1965, 1958, 1951,
1944, 1938, 1931, /* a=8 */
1925, 1918, 1912, 1905, 1899, 1892, 1886, 1880, 1874, 1867, 1861, 1855, 1849,
1844, 1838, 1832, /* a=9 */
1827, 1821, 1815, 1810, 1804, 1799, 1793, 1788, 1783, 1777, 1772, 1767, 1762,
1757, 1752, 1747, /* a=10 */
1742, 1737, 1733, 1728, 1723, 1718, 1713, 1709, 1704, 1699, 1695, 1690, 1686,
1681, 1677, 1673, /* a=11 */
1669, 1664, 1660, 1656, 1652, 1647, 1643, 1639, 1635, 1631, 1627, 1623, 1619,
1615, 1611, 1607, /* a=12 */
1603, 1600, 1596, 1592, 1588, 1585, 1581, 1577, 1574, 1570, 1566, 1563, 1559,
1556, 1552, 1549, /* a=13 */
1545, 1542, 1538, 1535, 1532, 1528, 1525, 1522, 1518, 1515, 1512, 1509, 1505,
1502, 1499, 1496, /* a=14 */
1493, 1490, 1487, 1484, 1481, 1478, 1475, 1472, 1469, 1466, 1463, 1460, 1457,
1454, 1451, 1449, /* a=15 */
1446, 1443, 1440, 1438, 1435, 1432, 1429, 1427, 1424, 1421, 1419, 1416, 1413,
1411, 1408, 1405, /* a=16 */
1403, 1400, 1398, 1395, 1393, 1390, 1388, 1385, 1383, 1380, 1378, 1375, 1373,
1371, 1368, 1366, /* a=17 */
1363, 1360, 1358, 1356, 1353, 1351, 1349, 1346, 1344, 1342, 1340, 1337, 1335,
1333, 1331, 1329, /* a=18 */
1327, 1325, 1323, 1321, 1319, 1316, 1314, 1312, 1310, 1308, 1306, 1304, 1302,
1300, 1298, 1296, /* a=19 */
1294, 1292, 1290, 1288, 1286, 1284, 1282, 1280, 1278, 1276, 1274, 1272, 1270,
1268, 1266, 1265, /* a=20 */
1262, 1260, 1258, 1256, 1254, 1253, 1251, 1249, 1247, 1245, 1244, 1242, 1240,
1238, 1236, 1235, /* a=21 */
1234, 1232, 1230, 1229, 1227, 1225, 1223, 1222, 1220, 1218, 1217, 1215, 1213,
1212, 1210, 1208, /* a=22 */
1206, 1204, 1203, 1201, 1199, 1198, 1196, 1195, 1193, 1191, 1190, 1188, 1187,
1185, 1184, 1182, /* a=23 */
1181, 1180, 1178, 1177, 1175, 1174, 1172, 1171, 1169, 1168, 1166, 1165, 1163,
1162, 1160, 1159, /* a=24 */
1157, 1156, 1154, 1153, 1151, 1150, 1149, 1147, 1146, 1144, 1143, 1142, 1140,
1139, 1137, 1136, /* a=25 */
1135, 1133, 1132, 1131, 1129, 1128, 1127, 1125, 1124, 1123, 1121, 1120, 1119,
1117, 1116, 1115, /* a=26 */
1114, 1113, 1111, 1110, 1109, 1108, 1106, 1105, 1104, 1103, 1101, 1100, 1099,
1098, 1096, 1095, /* a=27 */
1093, 1092, 1091, 1090, 1089, 1087, 1086, 1085, 1084, 1083, 1081, 1080, 1079,
1078, 1077, 1076, /* a=28 */
1075, 1073, 1072, 1071, 1070, 1069, 1068, 1067, 1065, 1064, 1063, 1062, 1061,
1060, 1059, 1058, /* a=29 */
1057, 1056, 1055, 1054, 1052, 1051, 1050, 1049, 1048, 1047, 1046, 1045, 1044,
1043, 1042, 1041, /* a=30 */
1040, 1039, 1038, 1037, 1036, 1035, 1034, 1033, 1032, 1031, 1030, 1029, 1028,
1027, 1026, 1025 /* a=31 */
};
  static unsigned char T2[384] = {
    7, 7, 6, 6, 5, 5, 4, 4, 4, 3, 3, 2, 2, 1, 1, 0, /* a=8 */
    6, 5, 5, 5, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 0, 0, /* a=9 */
    5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, /* a=10 */
    4, 4, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, /* a=11 */
    3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, /* a=12 */
    3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, /* a=13 */
    3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, /* a=14 */
    2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, /* a=15 */
    2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, /* a=16 */
    2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, /* a=17 */
    3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, /* a=18 */
    2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, /* a=19 */
    1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, /* a=20 */
    2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, /* a=21 */
    1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* a=22 */
    2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, /* a=23 */
    1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* a=24 */
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, /* a=25 */
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, /* a=26 */
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* a=27 */
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, /* a=28 */
    1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, /* a=29 */
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* a=30 */
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  /* a=31 */
};
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
      /* if one wants faithful rounding for p=11, replace #if 0 by #if 1 */
#if 0
      /* Note: we could have a table giving an error less than 1 ulp
         everywhere, but with t=2048 for i=2048, which would require
         a special treatment, since 2048 does not fit in 11 bits.
         We prefer to always require t <= 2047, with some slightly larger
         error for p=11. */
      switch (i)
        {
        case 2050:
          t = 2047; /* table gives 2046, error 1.000732 */
          break;
        case 2052:
          t = 2046; /* table gives 2045, error 1.002925 */
          break;
        case 2054:
          t = 2045; /* table gives 2044, error 1.006576 */
          break;
        case 2273:
          t = 1944; /* table gives 1945, error 1.004700 */
          break;
        default:
#else
        {
#endif
          ab = i >> 4;
          ac = (ab & 0x3F0) | (i & 0x0F);
          t = (mp_limb_t) T1[ab - 0x80] + (mp_limb_t) T2[ac - 0x80];
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
      cu = mpn_add_1 (u, u, sn, u[0] & (MPFR_LIMB_ONE << pl));
      /* mask bits 0..pl of u[0] */
      if (pl + 1 == GMP_NUMB_BITS)
        u[0] = 0;
      else
        u[0] &= ~MPFR_LIMB_MASK(pl + 1);

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
              cy = mpn_add_1 (x, x, n, MPFR_LIMB_ONE) - cy;
              /* we also need to subtract 1 at x[ln] */
              cy -= mpn_sub_1 (x + ln, x + ln, xn, MPFR_LIMB_ONE);
              /* n - ln = xn */
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
          cy -= mpn_sub_1 (x, x, n, MPFR_LIMB_ONE << pl);
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
  A[0] &= ~MPFR_LIMB_MASK(pl);
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
      A[0] &= ~MPFR_LIMB_MASK(pl);
      // gmp_printf ("\nA=%Nd\n", A, n + 1);
      mpfr_mpn_rec_sqrt (X, A, p);
      // gmp_printf ("A=%Nd X=%Nd\n", A, n + 1, X, n);
      check (X, A, p);
    }
#endif
  return 0;
}
#endif
