/*Copyright (c) 2012 The Broad Institute

*Permission is hereby granted, free of charge, to any person
*obtaining a copy of this software and associated documentation
*files (the "Software"), to deal in the Software without
*restriction, including without limitation the rights to use,
*copy, modify, merge, publish, distribute, sublicense, and/or sell
*copies of the Software, and to permit persons to whom the
*Software is furnished to do so, subject to the following
*conditions:

*The above copyright notice and this permission notice shall be
*included in all copies or substantial portions of the Software.

*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
*EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
*OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
*NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
*HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
*WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
*FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
*THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/



#undef SIMD_ENGINE
#undef SIMD_ENGINE_AVX

#define SIMD_ENGINE sse
#define SIMD_ENGINE_SSE

#include "template.h"

#include "define-sse-float.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

#include "define-sse-double.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

template double compute_full_prob_ssed<double>(testcase* tc, double* nextlog);
template float compute_full_prob_sses<float>(testcase* tc, float* nextlog);
