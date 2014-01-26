#undef SIMD_ENGINE
#undef SIMD_ENGINE_AVX
#undef SIMD_ENGINE_SSE

#define SIMD_ENGINE avx
#define SIMD_ENGINE_AVX

#include "define-float.h"
#include "vector_function_prototypes.h"

#include "define-double.h"
#include "vector_function_prototypes.h"

#undef  SIMD_ENGINE
#undef  SIMD_ENGINE_AVX

#define SIMD_ENGINE sse
#define SIMD_ENGINE_SSE


#include "define-sse-float.h"
#include "vector_function_prototypes.h"

#include "define-sse-double.h"
#include "vector_function_prototypes.h"

#undef SIMD_ENGINE
#undef SIMD_ENGINE_AVX
#undef SIMD_ENGINE_SSE

