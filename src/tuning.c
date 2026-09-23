#include "smallcas.h"

size_t sc_tune_mul_ks_cutoff = 16;
size_t sc_tune_mul_toom3_cutoff = 48;
size_t sc_tune_mul_karatsuba_cutoff = 12;
size_t sc_tune_mul_karatsuba_low_bits = 128;
size_t sc_tune_mul_karatsuba_low_bits_cutoff = 24;
size_t sc_tune_mul_ntt_cutoff = (size_t)-1;
size_t sc_tune_mul_ssa_cutoff = (size_t)-1;
