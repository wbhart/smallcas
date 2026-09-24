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
