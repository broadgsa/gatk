#include "template.h"

#undef SIMD_TYPE
#undef SIMD_TYPE_SSE

#define SIMD_TYPE avx
#define SIMD_TYPE_AVX

#include "define-float.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

#include "define-double.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

template double compute_full_prob_avxd<double>(testcase* tc, double* nextlog);
template float compute_full_prob_avxs<float>(testcase* tc, float* nextlog);
