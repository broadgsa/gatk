#include "template.h"

#undef SIMD_TYPE
#undef SIMD_TYPE_AVX

#define SIMD_TYPE sse
#define SIMD_TYPE_SSE

#include "define-sse-float.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

#include "define-sse-double.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

template double compute_full_prob_ssed<double>(testcase* tc, double* nextlog);
template float compute_full_prob_sses<float>(testcase* tc, float* nextlog);
