#include "smallcas.h"

size_t sc_tune_mul_ks_cutoff = 11;
size_t sc_tune_mul_toom3_cutoff = 24;
size_t sc_tune_mul_karatsuba_cutoff = 69;
size_t sc_tune_mul_karatsuba_low_bits = 128;
size_t sc_tune_mul_karatsuba_low_bits_cutoff = 91;
size_t sc_tune_mul_ntt_cutoff = 1748;
size_t sc_tune_mul_ssa_cutoff = 180;
unsigned sc_tune_ssa_mfa_cutoff_log = 15;

size_t sc_tune_mullow_dc_cutoff = 108;
size_t sc_tune_mullow_ntt_cutoff = 1170;
size_t sc_tune_mullow_ssa_cutoff = 151;
size_t sc_tune_mulhigh_ntt_cutoff = 983;
size_t sc_tune_mulhigh_ssa_cutoff = 135;
size_t sc_tune_mulmid_classical_cutoff = 33;
size_t sc_tune_mulmid_toom63_cutoff = 39;
size_t sc_tune_mulmid_ntt_cutoff = 651;
size_t sc_tune_mulmid_ssa_cutoff = 86;

size_t sc_tune_series_quo_dc_cutoff = 16;
size_t sc_tune_inv_series_newton_cutoff = 16;
size_t sc_tune_series_quo_newton_cutoff = 16;
size_t sc_tune_bidir_quo_cutoff = 16;
size_t sc_tune_mulders_quo_cutoff = 16;
size_t sc_tune_quo_mulders_cutoff = (size_t)-1;
size_t sc_tune_divrem_mulders_cutoff = (size_t)-1;
size_t sc_tune_divrem_dc_cutoff = 32;
size_t sc_tune_quo_newton_cutoff = 64;
size_t sc_tune_divrem_newton_cutoff = 64;
size_t sc_tune_quo_dc_cutoff = 32;
size_t sc_tune_divexact_bidir_cutoff = 32;
size_t sc_tune_pseudodiv_fast_cutoff = (size_t)-1;

size_t sc_tune_pseudorem_fast_cutoff = (size_t)-1;
size_t sc_tune_gcd_subresultant_cutoff = 5;
size_t sc_tune_resultant_subresultant_cutoff = 8;
size_t sc_tune_taylor_conv_cutoff = 3072;
