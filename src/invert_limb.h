/* __gmpfr_invert_limb -- implement GMP's invert_limb (which is not in GMP API)

Copyright 2016-2017 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* for now, we only provide __gmpfr_invert_limb for 64-bit limb */
#if GMP_NUMB_BITS == 64

/* for 256 <= d9 < 512, l[d9-256] = floor((2^19-3*2^8)/d9) */
static const mp_limb_t invert_limb_table[256] = {2045, 2037, 2029, 2021, 2013, 2005, 1998, 1990, 1983, 1975, 1968, 1960, 1953, 1946, 1938, 1931, 1924, 1917, 1910, 1903, 1896, 1889, 1883, 1876, 1869, 1863, 1856, 1849, 1843, 1836, 1830, 1824, 1817, 1811, 1805, 1799, 1792, 1786, 1780, 1774, 1768, 1762, 1756, 1750, 1745, 1739, 1733, 1727, 1722, 1716, 1710, 1705, 1699, 1694, 1688, 1683, 1677, 1672, 1667, 1661, 1656, 1651, 1646, 1641, 1636, 1630, 1625, 1620, 1615, 1610, 1605, 1600, 1596, 1591, 1586, 1581, 1576, 1572, 1567, 1562, 1558, 1553, 1548, 1544, 1539, 1535, 1530, 1526, 1521, 1517, 1513, 1508, 1504, 1500, 1495, 1491, 1487, 1483, 1478, 1474, 1470, 1466, 1462, 1458, 1454, 1450, 1446, 1442, 1438, 1434, 1430, 1426, 1422, 1418, 1414, 1411, 1407, 1403, 1399, 1396, 1392, 1388, 1384, 1381, 1377, 1374, 1370, 1366, 1363, 1359, 1356, 1352, 1349, 1345, 1342, 1338, 1335, 1332, 1328, 1325, 1322, 1318, 1315, 1312, 1308, 1305, 1302, 1299, 1295, 1292, 1289, 1286, 1283, 1280, 1276, 1273, 1270, 1267, 1264, 1261, 1258, 1255, 1252, 1249, 1246, 1243, 1240, 1237, 1234, 1231, 1228, 1226, 1223, 1220, 1217, 1214, 1211, 1209, 1206, 1203, 1200, 1197, 1195, 1192, 1189, 1187, 1184, 1181, 1179, 1176, 1173, 1171, 1168, 1165, 1163, 1160, 1158, 1155, 1153, 1150, 1148, 1145, 1143, 1140, 1138, 1135, 1133, 1130, 1128, 1125, 1123, 1121, 1118, 1116, 1113, 1111, 1109, 1106, 1104, 1102, 1099, 1097, 1095, 1092, 1090, 1088, 1086, 1083, 1081, 1079, 1077, 1074, 1072, 1070, 1068, 1066, 1064, 1061, 1059, 1057, 1055, 1053, 1051, 1049, 1047, 1044, 1042, 1040, 1038, 1036, 1034, 1032, 1030, 1028, 1026, 1024};

/* Implements Algorithm 2 from "Improved Division by Invariant Integers",
   Niels Möller and Torbjörn Granlund, IEEE Transactions on Computers,
   volume 60, number 2, pages 165-175, 2011. */
#ifdef HAVE_MULX_U64
#include <immintrin.h>
#define __gmpfr_invert_limb(r, d)                                       \
    do {                                                                \
      mp_limb_t _d, _d0, _d9, _d40, _d63, _v0, _v1, _v2, _e, _v3, _h, _l; \
      _d = (d);                                                         \
      _d9 = _d >> 55;                                                   \
      _v0 = invert_limb_table[_d9 - 256];                               \
      _d40 = (_d >> 24) + 1;                                            \
      _v1 = (_v0 << 11) - ((_v0 * _v0 * _d40) >> 40) - 1;               \
      _v2 = (_v1 << 13) + ((_v1 * (0x1000000000000000 - _v1 * _d40)) >> 47); \
      _d0 = _d & 1;                                                     \
      _d63 = ((_d - 1) >> 1) + 1;                                       \
      _e = - _v2 * _d63 + ((_v2 & -_d0) >> 1);                          \
      _mulx_u64 (_v2, _e, (unsigned long long *) &_h);                  \
      _v3 = (_v2 << 31) + (_h >> 1);                                    \
      umul_ppmm (_h, _l, _v3, d);                                       \
      /* v3 is too small iff (h+d)*2^64+l+d < 2^128 */                  \
      add_ssaaaa(_h, _l, _h, _l, _d, _d);                               \
      MPFR_ASSERTD(_h == MPFR_LIMB_ZERO || -_h == MPFR_LIMB_ONE);       \
      (r) = _v3 - _h;                                                   \
    } while (0)
#else
#define __gmpfr_invert_limb(r, d)                                       \
    do {                                                                \
      mp_limb_t _d, _d0, _d9, _d40, _d63, _v0, _v1, _v2, _e, _v3, _h, _l; \
      _d = (d);                                                         \
      _d9 = _d >> 55;                                                   \
      _v0 = invert_limb_table[_d9 - 256];                               \
      _d40 = (_d >> 24) + 1;                                            \
      _v1 = (_v0 << 11) - ((_v0 * _v0 * _d40) >> 40) - 1;               \
      _v2 = (_v1 << 13) + ((_v1 * (0x1000000000000000 - _v1 * _d40)) >> 47); \
      _d0 = _d & 1;                                                     \
      _d63 = ((_d - 1) >> 1) + 1;                                       \
      _e = - _v2 * _d63 + ((_v2 & -_d0) >> 1);                          \
      umul_ppmm (_h, _l, _v2, _e);                                      \
      _v3 = (_v2 << 31) + (_h >> 1);                                    \
      umul_ppmm (_h, _l, _v3, d);                                       \
      /* v3 is too small iff (h+d)*2^64+l+d < 2^128 */                  \
      add_ssaaaa(_h, _l, _h, _l, _d, _d);                               \
      MPFR_ASSERTD(_h == MPFR_LIMB_ZERO || -_h == MPFR_LIMB_ONE);       \
      (r) = _v3 - _h;                                                   \
    } while (0)
#endif /* HAVE_MULX_U64 */

#ifdef HAVE_MULX_U64
/* same algorithm, but return the value v3, which is such that
   v3 <= invert_limb (d) <= v3 + 1 */
#define __gmpfr_invert_limb_approx(r, d)                                \
    do {                                                                \
      mp_limb_t _d, _d0, _d9, _d40, _d63, _v0, _v1, _v2, _e, _h;        \
      _d = (d);                                                         \
      _d9 = _d >> 55;                                                   \
      _v0 = invert_limb_table[_d9 - 256];                               \
      _d40 = (_d >> 24) + 1;                                            \
      _v1 = (_v0 << 11) - ((_v0 * _v0 * _d40) >> 40) - 1;               \
      _v2 = (_v1 << 13) + ((_v1 * (0x1000000000000000 - _v1 * _d40)) >> 47); \
      _d0 = _d & 1;                                                     \
      _d63 = ((_d - 1) >> 1) + 1;                                       \
      _e = - _v2 * _d63 + ((_v2 & -_d0) >> 1);                          \
      _mulx_u64 (_v2, _e, (unsigned long long *) &_h);                  \
        (r) = (_v2 << 31) + (_h >> 1);                                  \
    } while (0)
#else
#define __gmpfr_invert_limb_approx(r, d)                                \
    do {                                                                \
      mp_limb_t _d, _d0, _d9, _d40, _d63, _v0, _v1, _v2, _e, _h, _l;    \
      _d = (d);                                                         \
      _d9 = _d >> 55;                                                   \
      _v0 = invert_limb_table[_d9 - 256];                               \
      _d40 = (_d >> 24) + 1;                                            \
      _v1 = (_v0 << 11) - ((_v0 * _v0 * _d40) >> 40) - 1;               \
      _v2 = (_v1 << 13) + ((_v1 * (0x1000000000000000 - _v1 * _d40)) >> 47); \
      _d0 = _d & 1;                                                     \
      _d63 = ((_d - 1) >> 1) + 1;                                       \
      _e = - _v2 * _d63 + ((_v2 & -_d0) >> 1);                          \
      umul_ppmm (_h, _l, _v2, _e);                                      \
      (r) = (_v2 << 31) + (_h >> 1);                                    \
    } while (0)
#endif /* HAVE_MULX_U64 */

#endif /* GMP_NUMB_BITS == 64 */
