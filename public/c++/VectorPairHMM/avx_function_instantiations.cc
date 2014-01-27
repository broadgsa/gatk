#include "template.h"

#undef SIMD_ENGINE
#undef SIMD_ENGINE_SSE

#define SIMD_ENGINE avx
#define SIMD_ENGINE_AVX

#include "define-float.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

#include "define-double.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

template double compute_full_prob_avxd<double>(testcase* tc, double* nextlog);
template float compute_full_prob_avxs<float>(testcase* tc, float* nextlog);

