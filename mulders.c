/* put in rp[n-1..2n-1] an approximation of the n+1 high limbs
   of {mp, n} * {np, n}.
   Assumes rp has 2n limbs.
*/
void
mpn_mul_hi_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  if (n < MPN_MUL_HI_THRESHOLD)
    mpn_mul_hi_basecase (rp, np, mp, n);
  else
    {
      mp_size_t k = MPN_MUL_HI_THRESHOLD_TABLE[n];
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_mul_n (rp + 2 * l, np + l, mp + l, k); /* fills rp[2l..2n-1] */
      mpn_mul_hi_n (rp, np + k, mp, l);          /* fills rp[l-1..2l-1] */
      cy = mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_mul_hi_n (rp, np, mp + k, l);          /* fills rp[l-1..2l-1] */
      cy += mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_add_1 (rp + n + l, rp + n + l, k, cy); /* propagate carry */
    }
}
