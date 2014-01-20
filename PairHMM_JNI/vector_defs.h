#undef SIMD_TYPE
#undef SIMD_TYPE_AVX
#undef SIMD_TYPE_SSE

#define SIMD_TYPE avx
#define SIMD_TYPE_AVX

#include "define-float.h"
#include "vector_function_prototypes.h"

#include "define-double.h"
#include "vector_function_prototypes.h"

#undef  SIMD_TYPE
#undef  SIMD_TYPE_AVX

#define SIMD_TYPE sse
#define SIMD_TYPE_SSE

#include "define-sse-float.h"
#include "vector_function_prototypes.h"

#include "define-sse-double.h"
#include "vector_function_prototypes.h"


