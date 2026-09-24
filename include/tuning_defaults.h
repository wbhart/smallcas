#ifndef SMALLCAS_TUNING_DEFAULTS_H
#define SMALLCAS_TUNING_DEFAULTS_H

#include <stddef.h>

#ifdef SC_TUNE
extern size_t sc_tune_mul_ks_cutoff;
extern size_t sc_tune_mul_toom3_cutoff;
extern size_t sc_tune_mul_karatsuba_cutoff;
extern size_t sc_tune_mul_karatsuba_low_bits;
extern size_t sc_tune_mul_karatsuba_low_bits_cutoff;
extern size_t sc_tune_mul_ntt_cutoff;
extern size_t sc_tune_mul_ssa_cutoff;
extern unsigned sc_tune_ssa_mfa_cutoff_log;
extern size_t sc_tune_mullow_dc_cutoff;
extern size_t sc_tune_mullow_ntt_cutoff;
extern size_t sc_tune_mullow_ssa_cutoff;
extern size_t sc_tune_mulhigh_ntt_cutoff;
extern size_t sc_tune_mulhigh_ssa_cutoff;
extern size_t sc_tune_mulmid_classical_cutoff;
extern size_t sc_tune_mulmid_toom63_cutoff;
extern size_t sc_tune_mulmid_ntt_cutoff;
extern size_t sc_tune_mulmid_ssa_cutoff;
extern size_t sc_tune_series_quo_dc_cutoff;
extern size_t sc_tune_inv_series_newton_cutoff;
extern size_t sc_tune_series_quo_newton_cutoff;
extern size_t sc_tune_bidir_quo_cutoff;
extern size_t sc_tune_mulders_quo_cutoff;
extern size_t sc_tune_quo_mulders_cutoff;
extern size_t sc_tune_divrem_mulders_cutoff;
extern size_t sc_tune_divrem_dc_cutoff;
extern size_t sc_tune_quo_newton_cutoff;
extern size_t sc_tune_divrem_newton_cutoff;
extern size_t sc_tune_quo_dc_cutoff;
extern size_t sc_tune_divexact_bidir_cutoff;
extern size_t sc_tune_pseudodiv_fast_cutoff;
extern size_t sc_tune_pseudorem_fast_cutoff;
extern size_t sc_tune_gcd_subresultant_cutoff;
#define SC_MUL_KS_CUTOFF sc_tune_mul_ks_cutoff
#define SC_MUL_TOOM3_CUTOFF sc_tune_mul_toom3_cutoff
#define SC_MUL_KARATSUBA_CUTOFF sc_tune_mul_karatsuba_cutoff
#define SC_MUL_KARATSUBA_LOW_BITS sc_tune_mul_karatsuba_low_bits
#define SC_MUL_KARATSUBA_LOW_BITS_CUTOFF sc_tune_mul_karatsuba_low_bits_cutoff
#define SC_MUL_NTT_CUTOFF sc_tune_mul_ntt_cutoff
#define SC_MUL_SSA_CUTOFF sc_tune_mul_ssa_cutoff
#define SC_SSA_MFA_CUTOFF_LOG sc_tune_ssa_mfa_cutoff_log
#define SC_MULLOW_DC_CUTOFF sc_tune_mullow_dc_cutoff
#define SC_MULLOW_NTT_CUTOFF sc_tune_mullow_ntt_cutoff
#define SC_MULLOW_SSA_CUTOFF sc_tune_mullow_ssa_cutoff
#define SC_MULHIGH_NTT_CUTOFF sc_tune_mulhigh_ntt_cutoff
#define SC_MULHIGH_SSA_CUTOFF sc_tune_mulhigh_ssa_cutoff
#define SC_MULMID_CLASSICAL_CUTOFF sc_tune_mulmid_classical_cutoff
#define SC_MULMID_TOOM63_CUTOFF sc_tune_mulmid_toom63_cutoff
#define SC_MULMID_NTT_CUTOFF sc_tune_mulmid_ntt_cutoff
#define SC_MULMID_SSA_CUTOFF sc_tune_mulmid_ssa_cutoff
#define SC_SERIES_QUO_DC_CUTOFF sc_tune_series_quo_dc_cutoff
#define SC_INV_SERIES_NEWTON_CUTOFF sc_tune_inv_series_newton_cutoff
#define SC_SERIES_QUO_NEWTON_CUTOFF sc_tune_series_quo_newton_cutoff
#define SC_BIDIR_QUO_CUTOFF sc_tune_bidir_quo_cutoff
#define SC_MULDERS_QUO_CUTOFF sc_tune_mulders_quo_cutoff
#define SC_QUO_MULDERS_CUTOFF sc_tune_quo_mulders_cutoff
#define SC_DIVREM_MULDERS_CUTOFF sc_tune_divrem_mulders_cutoff
#define SC_DIVREM_DC_CUTOFF sc_tune_divrem_dc_cutoff
#define SC_QUO_NEWTON_CUTOFF sc_tune_quo_newton_cutoff
#define SC_DIVREM_NEWTON_CUTOFF sc_tune_divrem_newton_cutoff
#define SC_QUO_DC_CUTOFF sc_tune_quo_dc_cutoff
#define SC_DIVEXACT_BIDIR_CUTOFF sc_tune_divexact_bidir_cutoff
#define SC_PSEUDODIV_FAST_CUTOFF sc_tune_pseudodiv_fast_cutoff
#define SC_PSEUDOREM_FAST_CUTOFF sc_tune_pseudorem_fast_cutoff
#define SC_GCD_SUBRESULTANT_CUTOFF sc_tune_gcd_subresultant_cutoff
#else
#ifndef SC_MUL_KS_CUTOFF
#define SC_MUL_KS_CUTOFF ((size_t)11)
#endif
#ifndef SC_MUL_TOOM3_CUTOFF
#define SC_MUL_TOOM3_CUTOFF ((size_t)24)
#endif
#ifndef SC_MUL_KARATSUBA_CUTOFF
#define SC_MUL_KARATSUBA_CUTOFF ((size_t)69)
#endif
#ifndef SC_MUL_KARATSUBA_LOW_BITS
#define SC_MUL_KARATSUBA_LOW_BITS ((size_t)128)
#endif
#ifndef SC_MUL_KARATSUBA_LOW_BITS_CUTOFF
#define SC_MUL_KARATSUBA_LOW_BITS_CUTOFF ((size_t)91)
#endif
#ifndef SC_MUL_NTT_CUTOFF
/* NTT uses shorter length / number of CRT primes. */
#define SC_MUL_NTT_CUTOFF ((size_t)1748)
#endif
#ifndef SC_MUL_SSA_CUTOFF
#define SC_MUL_SSA_CUTOFF ((size_t)180)
#endif
#ifndef SC_SSA_MFA_CUTOFF_LOG
#define SC_SSA_MFA_CUTOFF_LOG ((unsigned)15)
#endif
#endif

#ifndef SC_FFT_MFA_BASE_LOG
#define SC_FFT_MFA_BASE_LOG ((unsigned)8)
#endif
#ifndef SC_MULLOW_DC_CUTOFF
#define SC_MULLOW_DC_CUTOFF ((size_t)108)
#endif
#ifndef SC_MULLOW_NTT_CUTOFF
#define SC_MULLOW_NTT_CUTOFF ((size_t)1170)
#endif
#ifndef SC_MULLOW_SSA_CUTOFF
#define SC_MULLOW_SSA_CUTOFF ((size_t)151)
#endif
#ifndef SC_MULHIGH_NTT_CUTOFF
#define SC_MULHIGH_NTT_CUTOFF ((size_t)983)
#endif
#ifndef SC_MULHIGH_SSA_CUTOFF
#define SC_MULHIGH_SSA_CUTOFF ((size_t)135)
#endif
#ifndef SC_MULMID_CLASSICAL_CUTOFF
#define SC_MULMID_CLASSICAL_CUTOFF ((size_t)33)
#endif
#ifndef SC_MULMID_TOOM63_CUTOFF
#define SC_MULMID_TOOM63_CUTOFF ((size_t)39)
#endif
#ifndef SC_MULMID_NTT_CUTOFF
#define SC_MULMID_NTT_CUTOFF ((size_t)651)
#endif
#ifndef SC_MULMID_SSA_CUTOFF
#define SC_MULMID_SSA_CUTOFF ((size_t)86)
#endif
#ifndef SC_MULMID_GENERAL_CLASSICAL_CUTOFF
#define SC_MULMID_GENERAL_CLASSICAL_CUTOFF ((size_t)24)
#endif
#ifndef SC_SERIES_QUO_DC_CUTOFF
#define SC_SERIES_QUO_DC_CUTOFF ((size_t)16)
#endif
#ifndef SC_INV_SERIES_NEWTON_CUTOFF
#define SC_INV_SERIES_NEWTON_CUTOFF ((size_t)16)
#endif
#ifndef SC_SERIES_QUO_NEWTON_CUTOFF
#define SC_SERIES_QUO_NEWTON_CUTOFF ((size_t)16)
#endif
#ifndef SC_BIDIR_QUO_CUTOFF
#define SC_BIDIR_QUO_CUTOFF ((size_t)16)
#endif
#ifndef SC_MULDERS_QUO_CUTOFF
#define SC_MULDERS_QUO_CUTOFF ((size_t)16)
#endif
#ifndef SC_QUO_MULDERS_CUTOFF
#define SC_QUO_MULDERS_CUTOFF ((size_t)-1)
#endif
#ifndef SC_DIVREM_MULDERS_CUTOFF
#define SC_DIVREM_MULDERS_CUTOFF ((size_t)-1)
#endif
#ifndef SC_DIVREM_DC_CUTOFF
#define SC_DIVREM_DC_CUTOFF ((size_t)32)
#endif
#ifndef SC_QUO_NEWTON_CUTOFF
#define SC_QUO_NEWTON_CUTOFF ((size_t)64)
#endif
#ifndef SC_DIVREM_NEWTON_CUTOFF
#define SC_DIVREM_NEWTON_CUTOFF ((size_t)64)
#endif
#ifndef SC_QUO_DC_CUTOFF
#define SC_QUO_DC_CUTOFF ((size_t)32)
#endif
#ifndef SC_DIVEXACT_BIDIR_CUTOFF
#define SC_DIVEXACT_BIDIR_CUTOFF ((size_t)32)
#endif
#ifndef SC_PSEUDODIV_FAST_CUTOFF
#define SC_PSEUDODIV_FAST_CUTOFF ((size_t)-1)
#endif
#ifndef SC_PSEUDOREM_FAST_CUTOFF
#define SC_PSEUDOREM_FAST_CUTOFF ((size_t)-1)
#endif
#ifndef SC_GCD_SUBRESULTANT_CUTOFF
#define SC_GCD_SUBRESULTANT_CUTOFF ((size_t)5)
#endif
#ifndef SC_HGCD_BASE_CUTOFF
#define SC_HGCD_BASE_CUTOFF ((size_t)12)
#endif
#ifndef SC_EVALUATE_DC_CUTOFF
#define SC_EVALUATE_DC_CUTOFF ((size_t)16)
#endif
#ifndef SC_COMPOSE_DC_CUTOFF
#define SC_COMPOSE_DC_CUTOFF ((size_t)8)
#endif
#ifndef SC_TAYLOR_DC_CUTOFF
#define SC_TAYLOR_DC_CUTOFF ((size_t)4096)
#endif

#endif
