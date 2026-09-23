#ifndef SMALLCAS_TUNING_H
#define SMALLCAS_TUNING_H

#include <stddef.h>

#ifdef SC_TUNE
extern size_t sc_tune_mul_ks_cutoff;
extern size_t sc_tune_mul_toom3_cutoff;
extern size_t sc_tune_mul_karatsuba_cutoff;
extern size_t sc_tune_mul_karatsuba_low_bits;
extern size_t sc_tune_mul_karatsuba_low_bits_cutoff;
extern size_t sc_tune_mul_ntt_cutoff;
extern size_t sc_tune_mul_ssa_cutoff;
#define SC_MUL_KS_CUTOFF sc_tune_mul_ks_cutoff
#define SC_MUL_TOOM3_CUTOFF sc_tune_mul_toom3_cutoff
#define SC_MUL_KARATSUBA_CUTOFF sc_tune_mul_karatsuba_cutoff
#define SC_MUL_KARATSUBA_LOW_BITS sc_tune_mul_karatsuba_low_bits
#define SC_MUL_KARATSUBA_LOW_BITS_CUTOFF sc_tune_mul_karatsuba_low_bits_cutoff
#define SC_MUL_NTT_CUTOFF sc_tune_mul_ntt_cutoff
#define SC_MUL_SSA_CUTOFF sc_tune_mul_ssa_cutoff
#else
#define SC_MUL_KS_CUTOFF ((size_t)11)
#define SC_MUL_TOOM3_CUTOFF ((size_t)24)
#define SC_MUL_KARATSUBA_CUTOFF ((size_t)69)
#define SC_MUL_KARATSUBA_LOW_BITS ((size_t)128)
#define SC_MUL_KARATSUBA_LOW_BITS_CUTOFF ((size_t)91)
/* NTT uses shorter length / number of CRT primes. */
#define SC_MUL_NTT_CUTOFF ((size_t)1748)
#define SC_MUL_SSA_CUTOFF ((size_t)180)
#endif

#define SC_FFT_MFA_BASE_LOG ((unsigned)8)
#define SC_SSA_MFA_CUTOFF_LOG ((unsigned)15)
#define SC_MULLOW_DC_CUTOFF ((size_t)16)
#define SC_MULMID_CLASSICAL_CUTOFF ((size_t)16)
#define SC_MULMID_TOOM63_CUTOFF ((size_t)48)
#define SC_MULMID_GENERAL_CLASSICAL_CUTOFF ((size_t)24)
#define SC_SERIES_QUO_DC_CUTOFF ((size_t)16)
#define SC_INV_SERIES_NEWTON_CUTOFF ((size_t)16)
#define SC_SERIES_QUO_NEWTON_CUTOFF ((size_t)16)
#define SC_BIDIR_QUO_CUTOFF ((size_t)16)
#define SC_MULDERS_QUO_CUTOFF ((size_t)16)
#define SC_DIVREM_DC_CUTOFF ((size_t)32)
#define SC_QUO_NEWTON_CUTOFF ((size_t)64)
#define SC_QUO_DC_CUTOFF ((size_t)32)
#define SC_DIVEXACT_BIDIR_CUTOFF ((size_t)32)
#define SC_HGCD_BASE_CUTOFF ((size_t)12)
#define SC_EVALUATE_DC_CUTOFF ((size_t)16)
#define SC_COMPOSE_DC_CUTOFF ((size_t)8)
#define SC_TAYLOR_DC_CUTOFF ((size_t)4096)

#endif
