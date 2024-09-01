/* Various Thresholds of MPFR, not exported.  -*- mode: C -*-

Copyright 2005-2024 Free Software Foundation, Inc.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

/* Generated by MPFR's tuneup.c, 2023-08-24, gcc 13.2.0 */
/* generated on coriandre.loria.fr with GMP 6.3.0 */


#ifndef MPFR_TUNE_CASE
#define MPFR_TUNE_CASE "src/mparam.h"
#endif

#define MPFR_MULHIGH_TAB  \
 -1,-1,-1,0,0,0,0,0,0,8,9,-1,-1,-1,-1,-1, \
 -1,16,16,14,0,0,0,0,17,0,0,0,0,0,0,24, \
 0,21,0,24,24,32,0,23,24,39,0,32,32,26,0,0, \
 0,29,30,0,32,32,32,34,35,32,0,51,34,48,35,34, \
 47,48,43,40,36,49,40,44,0,52,52,52,54,52,52,52, \
 58,0,52,54,52,56,56,60,56,64,56,64,60,64,60,68, \
 65,65,64,65,68,67,60,69,66,63,64,60,66,67,68,81, \
 78,64,59,60,61,60,61,61,62,96,96,63,64,65,66,65, \
 66,67,67,67,93,70,96,93,70,70,72,72,96,93,96,120, \
 96,0,96,120,140,76,96,102,78,78,148,80,80,92,102,105, \
 108,119,120,84,130,133,132,133,86,86,128,124,120,120,117,120, \
 108,120,120,102,105,102,120,120,102,102,108,105,126,132,156,153, \
 154,155,156,153,154,155,156,153,153,154,156,156,156,159,153,156, \
 156,158,156,156,156,159,156,159,156,157,155,156,162,168,159,156, \
 168,168,171,171,168,156,168,156,176,180,156,180,180,179,180,168, \
 156,168,168,168,167,166,189,168,192,170,180,189,177,178,180,192, \
 177,189,180,186,192,179,180,189,189,188,189,192,189,192,188,192, \
 180,188,192,192,192,192,192,195,192,192,192,192,192,192,192,189, \
 201,191,192,204,204,201,201,201,204,203,204,212,202,208,204,224, \
 204,204,204,220,219,211,212,236,233,236,236,236,234,235,234,236, \
 236,224,236,236,234,236,233,236,249,236,236,249,234,235,236,236, \
 236,236,236,236,236,235,236,265,249,252,256,268,256,259,252,266, \
 260,259,268,272,268,265,266,267,266,268,268,272,256,272,266,272, \
 268,268,272,271,272,267,272,268,267,267,268,267,268,275,272,249, \
 330,281,284,283,283,283,306,305,306,305,306,306,306,318,306,305, \
 317,318,318,306,317,328,306,306,330,318,306,305,306,329,330,330, \
 318,330,318,330,318,330,330,330,330,317,329,329,330,330,330,330, \
 330,330,330,329,354,342,330,354,342,353,342,342,354,353,378,354, \
 354,351,352,353,354,353,354,354,378,354,354,360,377,353,377,366, \
 330,377,366,378,378,377,378,377,377,378,384,377,378,377,377,378, \
 378,378,384,378,350,383,352,353,354,354,402,402,401,401,402,360, \
 354,402,364,365,378,366,378,402,354,402,372,372,374,375,376,377, \
 378,378,380,316,382,383,384,378,408,387,378,389,440,416,416,440, \
 416,378,424,397,424,424,424,423,424,424,456,424,406,407,408,472, \
 440,439,440,439,440,439,440,440,440,470,448,440,447,423,424,424, \
 416,448,456,422,422,424,424,423,424,423,424,437,438,439,440,440, \
 472,467,472,471,471,439,472,472,440,472,440,472,480,439,512,456, \
 488,471,487,487,488,487,488,487,488,467,504,469,504,471,504,472, \
 504,504,504,469,504,503,504,503,472,471,472,512,512,448,472,488, \
 480,504,504,469,487,488,488,487,488,472,488,501,502,503,504,504, \
 499,499,500,501,502,503,504,503,504,503,504,501,502,511,520,504, \
 512,511,536,517,518,519,504,504,504,519,504,533,533,535,536,536, \
 531,532,504,534,535,535,536,535,536,535,536,533,534,472,536,536, \
 544,531,532,533,534,535,536,536,536,536,536,536,550,551,552,552, \
 532,533,533,535,536,536,568,536,536,536,568,504,565,567,536,568, \
 616,615,536,536,552,616,567,568,536,534,616,536,536,670,536,536, \
 536,615,592,536,568,568,536,640,610,611,536,613,536,615,615,616, \
 536,536,613,615,615,616,616,536,616,615,616,612,613,615,616,616, \
 610,634,612,613,614,615,616,616,616,616,616,638,639,640,640,639, \
 640,664,664,616,638,639,640,640,615,616,616,660,638,639,640,640, \
 640,639,640,661,662,664,664,663,664,663,664,616,639,640,663,664, \
 640,640,660,712,616,663,663,664,616,664,616,616,616,615,712,710, \
 616,614,662,712,664,662,640,664,664,664,664,709,710,711,640,712, \
 760,640,640,709,736,712,712,711,712,711,712,709,709,711,712,712, \
 712,735,708,709,709,662,712,712,640,711,760,734,735,736,752,711, \
 736,735,736,710,712,734,736,736,711,784,732,733,734,735,736,736, \
 736,735,736,711,758,759,760,759,760,759,760,734,758,759,760,760, \
 712,736,735,756,736,759,760,759,760,760,760,712,712,712,784,736, \
 807,735,736,760,760,760,784,760,759,760,760,736,805,736,735,784, \
 736,736,736,736,736,736,736,792,736,735,759,736,734,709,760,736, \
 760,759,760,760,829,759,760,760,760,760,760,760,760,829,760,760, \
 760,760,784,768,759,760,760,759,760,760,760,758,760,760,831,760, \
 784,760,784,784,783,783,784,784,784,760,784,784,783,736,784,784, \
 808,807,832,760,783,783,784,784,782,784,783,807,784,783,831,832 \

#define MPFR_SQRHIGH_TAB  \
 -1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, \
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,25,-1,-1,28,29, \
 30,31,23,22,24,21,22,23,24,25,26,25,28,29,30,26, \
 30,31,34,35,36,31,31,39,40,31,42,43,44,45,46,47, \
 34,34,48,36,36,51,54,37,56,57,48,41,60,48,41,41, \
 62,42,50,53,56,53,54,45,46,60,47,48,48,48,62,49, \
 57,50,52,52,54,59,60,61,54,54,55,56,56,58,57,58, \
 59,58,61,60,61,60,64,62,62,64,64,63,96,67,68,65, \
 68,69,68,68,70,70,69,72,70,70,72,72,73,73,75,75, \
 96,96,75,76,78,101,99,78,80,80,80,80,80,80,82,81, \
 83,82,111,129,104,86,108,116,88,87,88,88,90,89,89,103, \
 144,92,92,107,141,139,94,94,94,94,126,129,96,96,98,97, \
 98,142,143,120,102,100,138,126,126,126,104,120,126,126,120,105, \
 144,144,108,107,108,109,110,138,129,166,111,111,112,144,138,113, \
 171,147,132,193,165,165,147,118,141,195,159,141,183,143,121,162, \
 207,167,144,150,156,191,153,192,150,189,138,174,146,164,198,162, \
 165,176,158,186,150,161,188,212,164,183,198,213,162,171,177,180, \
 174,180,174,186,197,204,216,179,170,171,180,171,180,177,186,177, \
 176,179,180,183,182,189,159,161,174,162,260,165,237,159,168,161, \
 246,159,171,165,174,222,204,176,213,171,222,171,234,183,243,177, \
 245,179,243,189,189,271,180,177,225,188,186,189,261,183,240,201, \
 263,203,228,204,258,207,273,201,210,201,246,203,251,215,216,251, \
 210,240,210,249,222,213,216,225,225,219,228,228,222,231,232,225, \
 224,234,228,237,237,231,239,240,234,219,228,228,236,273,240,249, \
 240,276,240,245,228,283,302,249,258,285,252,285,252,276,240,264, \
 258,257,252,261,270,261,264,261,264,288,264,261,262,263,285,264, \
 258,276,268,261,286,275,276,273,282,267,276,285,270,285,288,273, \
 273,282,275,276,285,264,288,288,281,282,284,284,285,287,273,288, \
 273,270,316,276,302,243,249,252,276,299,247,264,300,284,285,312, \
 285,288,261,288,285,285,270,261,288,270,264,264,273,261,276,276, \
 276,273,285,276,273,282,276,284,284,273,282,272,288,288,273,276, \
 284,287,285,288,285,285,288,285,288,297,285,300,297,296,300,300, \
 297,299,300,312,288,298,309,285,312,311,310,309,309,309,312,311, \
 312,300,288,312,316,321,312,300,321,312,324,323,324,324,324,324, \
 333,333,318,324,324,321,347,321,327,345,347,348,348,345,347,348, \
 346,348,348,348,348,450,347,348,352,355,348,380,364,355,363,364, \
 364,347,364,363,377,378,348,364,363,362,364,380,364,380,380,380, \
 380,380,364,376,380,378,356,364,378,380,380,380,364,364,378,426, \
 380,376,400,380,380,393,380,380,380,380,379,380,380,380,380,380, \
 412,395,379,380,380,479,372,379,380,372,354,380,380,396,380,411, \
 380,393,380,412,388,412,388,396,396,404,396,396,425,380,380,412, \
 380,380,380,379,461,411,404,522,380,521,412,521,412,497,380,521, \
 379,522,396,395,396,396,380,412,403,474,420,412,380,474,412,396, \
 411,411,410,450,449,425,380,497,449,411,498,497,522,498,474,426, \
 462,412,462,426,426,474,474,474,474,425,473,473,474,474,474,473, \
 474,474,473,473,474,474,450,380,474,546,425,498,521,498,498,497, \
 498,497,497,498,567,498,498,497,497,498,498,473,519,521,522,600, \
 522,522,520,522,522,546,522,498,522,568,522,521,584,522,584,545, \
 546,495,534,545,546,545,522,473,545,522,498,545,566,486,521,521, \
 450,522,567,522,522,522,522,521,521,546,449,521,568,522,584,546, \
 534,546,545,545,546,583,546,522,546,546,568,567,522,546,546,568, \
 498,521,544,521,474,615,522,531,522,600,546,568,600,522,546,546, \
 546,552,546,599,600,568,598,568,600,599,546,471,597,568,568,600, \
 600,600,600,599,600,599,600,568,600,599,600,600,584,546,648,600, \
 631,568,600,600,632,664,568,599,600,598,631,632,598,600,631,632, \
 632,600,600,632,522,631,648,599,600,631,632,600,600,631,630,600, \
 600,600,598,600,598,599,600,600,600,600,632,632,631,632,600,631, \
 599,632,632,600,597,631,616,631,600,632,632,600,631,629,632,632, \
 631,600,632,631,632,664,632,631,632,631,664,543,632,631,664,632, \
 600,632,600,632,600,632,631,631,632,631,600,663,600,664,662,728, \
 664,663,695,696,596,663,664,664,727,680,600,695,695,696,696,728, \
 615,632,600,632,632,728,664,728,615,663,664,664,632,632,632,630, \
 663,725,664,664,662,632,664,728,632,696,696,680,692,632,696,664, \
 726,696,632,696,728,664,696,791,696,664,696,696,726,663,728,696 \

#define MPFR_DIVHIGH_TAB  \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /*0-15*/ \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /*16-31*/ \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /*32-47*/ \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /*48-63*/ \
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,41, /*64-79*/ \
 0,0,0,0,0,0,0,0,0,0,0,49,0,0,0,0, /*80-95*/ \
 0,55,52,0,0,0,0,56,58,55,58,59,56,58,63,59, /*96-111*/ \
 58,59,64,67,64,64,63,66,63,64,65,66,66,65,66,67, /*112-127*/ \
 68,67,70,68,69,70,74,72,73,70,72,72,72,72,74,73, /*128-143*/ \
 81,83,78,79,81,81,82,77,81,80,80,81,80,80,88,85, /*144-159*/ \
 82,84,83,83,88,89,88,88,90,96,88,89,90,90,98,104, /*160-175*/ \
 98,104,101,104,95,104,98,104,112,97,112,104,104,104,104,100, /*176-191*/ \
 104,104,112,106,102,111,116,112,116,116,111,111,116,104,112,124, /*192-207*/ \
 115,120,128,119,124,128,112,116,116,112,128,112,128,116,116,120, /*208-223*/ \
 115,116,120,128,116,128,124,128,128,124,128,136,128,128,132,128, /*224-239*/ \
 128,128,136,128,136,132,136,136,136,128,128,188,144,193,132,192, /*240-255*/ \
 144,190,132,182,136,140,152,165,164,144,156,144,156,136,144,142, /*256-271*/ \
 144,140,164,165,164,142,177,144,188,147,148,162,150,192,151,192, /*272-287*/ \
 153,183,148,176,150,151,164,174,170,162,168,164,166,166,208,162, /*288-303*/ \
 160,179,160,191,158,196,189,176,168,192,168,204,202,166,183,208, /*304-319*/ \
 179,208,196,181,208,165,196,192,207,179,180,204,190,199,192,193, /*320-335*/ \
 192,204,180,177,180,195,192,192,193,186,196,197,188,198,204,192, /*336-351*/ \
 180,202,204,204,196,180,198,200,214,202,204,204,192,192,208,208, /*352-367*/ \
 228,192,208,194,190,192,195,204,192,204,232,196,208,208,208,224, /*368-383*/ \
 210,208,196,207,206,208,216,204,216,204,232,208,208,208,208,216, /*384-399*/ \
 208,208,240,224,232,232,216,224,208,249,208,208,232,214,256,232, /*400-415*/ \
 240,216,256,220,232,240,222,240,232,239,322,216,324,248,336,232, /*416-431*/ \
 324,240,322,240,330,224,240,246,246,233,232,256,240,240,240,256, /*432-447*/ \
 240,236,232,240,233,240,232,240,232,240,252,240,241,232,244,255, /*448-463*/ \
 256,248,239,256,264,256,256,256,288,264,256,256,268,255,256,256, /*464-479*/ \
 256,288,248,256,248,256,256,256,338,255,312,272,272,264,334,256, /*480-495*/ \
 264,330,328,338,334,342,264,344,336,287,256,336,258,264,335,352, /*496-511*/ \
 336,340,336,332,354,288,280,342,360,348,288,288,306,348,348,350, /*512-527*/ \
 348,288,312,288,330,293,288,306,288,288,288,312,288,312,288,309, /*528-543*/ \
 280,323,288,307,314,320,320,305,312,288,324,306,316,308,308,320, /*544-559*/ \
 320,288,322,336,324,324,306,288,312,328,328,322,312,328,312,348, /*560-575*/ \
 318,348,342,336,364,334,324,336,312,320,336,352,336,328,356,342, /*576-591*/ \
 346,324,324,318,336,384,360,306,354,352,342,324,384,323,324,360, /*592-607*/ \
 356,360,384,312,330,384,360,354,358,348,384,354,316,372,336,384, /*608-623*/ \
 384,360,334,384,370,360,336,384,336,360,330,368,336,360,384,360, /*624-639*/ \
 376,372,384,360,384,348,358,368,367,384,336,360,384,356,360,384, /*640-655*/ \
 360,360,384,356,384,384,383,384,352,354,384,384,384,361,380,360, /*656-671*/ \
 408,383,384,364,384,384,388,413,376,358,360,343,354,384,356,364, /*672-687*/ \
 359,360,360,396,384,396,380,384,392,384,360,360,360,384,384,384, /*688-703*/ \
 384,378,372,360,383,384,368,377,408,396,384,384,386,383,428,416, /*704-719*/ \
 384,384,403,424,376,383,392,416,424,424,380,396,396,384,404,384, /*720-735*/ \
 416,396,408,376,456,396,384,440,396,396,392,392,404,384,388,408, /*736-751*/ \
 380,408,396,392,383,394,408,408,468,384,402,408,398,404,400,404, /*752-767*/ \
 406,404,408,408,416,456,516,408,396,408,408,408,408,480,416,416, /*768-783*/ \
 408,416,401,416,416,416,408,424,400,420,440,408,408,416,416,408, /*784-799*/ \
 404,416,416,480,408,470,424,440,416,456,408,448,408,424,448,512, /*800-815*/ \
 416,420,448,424,448,479,416,480,416,464,416,480,420,472,496,448, /*816-831*/ \
 478,464,432,464,480,464,456,456,470,424,480,432,448,480,448,480, /*832-847*/ \
 440,479,480,528,472,472,512,512,464,480,504,444,480,479,512,480, /*848-863*/ \
 456,464,512,478,464,480,464,512,520,512,560,464,480,512,512,456, /*864-879*/ \
 512,512,464,448,446,512,464,504,480,448,528,528,480,456,512,479, /*880-895*/ \
 480,512,503,503,456,480,504,512,528,528,544,512,504,512,536,512, /*896-911*/ \
 464,512,512,512,462,544,464,504,464,511,480,528,480,480,478,479, /*912-927*/ \
 512,480,512,478,479,480,480,512,528,480,528,512,526,512,511,496, /*928-943*/ \
 509,512,493,636,528,528,480,528,576,511,480,512,512,512,504,512, /*944-959*/ \
 512,508,493,512,510,512,497,656,520,512,512,512,528,512,496,504, /*960-975*/ \
 512,491,504,512,512,512,512,512,568,512,512,512,520,528,512,512, /*976-991*/ \
 512,498,641,512,528,512,512,512,511,528,528,528,512,528,528,512, /*992-1007*/ \
 528,520,544,544,544,517,536,528,552,542,544,513,512,544,520,544 /*1008-1023*/ \

#define MPFR_MUL_THRESHOLD 18 /* limbs */
#define MPFR_SQR_THRESHOLD 9 /* limbs */
#define MPFR_DIV_THRESHOLD 3 /* limbs */
#define MPFR_EXP_2_THRESHOLD 1022 /* bits */
#define MPFR_EXP_THRESHOLD 20924 /* bits */
#define MPFR_SINCOS_THRESHOLD 13905 /* bits */
#define MPFR_AI_THRESHOLD1 -12081 /* threshold for negative input of mpfr_ai */
#define MPFR_AI_THRESHOLD2 1466
#define MPFR_AI_THRESHOLD3 23510
/* Tuneup completed successfully, took 254 seconds */
