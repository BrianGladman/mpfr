/* format for out_raw:
   A mpfr_t is represented by up to 3 fields, each one is represented by a
   sequence of 32-bit words. 32-bit words are stored as 4 bytes in little
   endian format:
   (a) a field for sign, precision and other bit-fields:
       - the sign (1 bit)
       - 2 bits for NaN, Inf, 0, normal numbers:
         00: 0
         01: normal
         10: Inf
         11: NaN
       - 1 bit for the precision encoding (prec_enc)
       - 1 bit for the exponent encoding (exp_enc) for normal numbers,
         otherwise 0 is stored here (reserved for future extensions)
       If prec_enc=0, the remaining 27 bits encode the precision (< 2^27)
       If prec_enc=1, the precision is stored in the following 27 bits
       (high part) and then 32 bits (low part). Thus the maximal precision
       is 59 bits.
   (b) (optional) a field for the exponent:
       - if the number is NaN, Inf, 0, this field is empty
       - if exp_enc=0, this field contains one 32-bit (signed) word encoding
         the exponent
       - if exp_enc=1, a first 32-bit word encodes a positive integer m,
         and the following m 32-bit words encode the exponent (in 2-complement
	 representation, with least significant words first)
   (c) (optional) a field for the significand:
       - if the number is NaN, Inf, 0, this field is empty
       - otherwise, let p = ceil(prec/32), the significand is represented
         by p consecutive 32-bit words (least significant words first).
	 Thus on a little-endian machine the significand can be directly
	 copied using memcopy.
   Examples:
   - a normal binary32 IEEE-754 number uses 96 bits: 32 for (a), 32 for (b),
     and 32 for (c);
   - a normal binary64 IEEE-754 number uses 128 bits: 32 for (a), 32 for (b),
     and 64 for (c) (idem for a significand of 64 bits, as in Intel x86
     double-extended format);
   - a normal binary128 IEEE-754 number uses 192 bits: 32 for (a), 32 for (b),
     and 128 for (c).
 */   
