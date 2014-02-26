#include "template.h"

#undef SIMD_ENGINE
#undef SIMD_ENGINE_AVX

#define SIMD_ENGINE sse
#define SIMD_ENGINE_SSE

#include "define-sse-float.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

#include "define-sse-double.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

template double compute_full_prob_ssed<double>(testcase* tc, double* nextlog);
template float compute_full_prob_sses<float>(testcase* tc, float* nextlog);
