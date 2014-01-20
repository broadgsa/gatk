#include <iostream>

#ifdef PRECISION
	#undef PRECISION
        #undef MAIN_TYPE
        #undef MAIN_TYPE_SIZE
        #undef UNION_TYPE
        #undef IF_128
        #undef IF_MAIN_TYPE
        #undef SHIFT_CONST1
        #undef SHIFT_CONST2
        #undef SHIFT_CONST3
        #undef _128_TYPE
        #undef _256_TYPE
        #undef AVX_LENGTH
        #undef MAVX_COUNT
	#undef HAP_TYPE
        #undef MASK_TYPE
        #undef MASK_ALL_ONES

	#undef SET_VEC_ZERO(__vec)
	#undef VEC_OR(__v1, __v2)
	#undef VEC_ADD(__v1, __v2)
        #undef VEC_SUB(__v1, __v2)
	#undef VEC_MUL(__v1, __v2)
	#undef VEC_DIV(__v1, __v2)
	#undef VEC_BLEND(__v1, __v2, __mask)
	#undef VEC_BLENDV(__v1, __v2, __maskV)
        #undef VEC_CAST_256_128(__v1)
        #undef VEC_EXTRACT_128(__v1, __im)
        #undef VEC_EXTRACT_UNIT(__v1, __im)
        #undef VEC_SET1_VAL128(__val)
        #undef VEC_MOVE(__v1, __val)
        #undef VEC_CAST_128_256(__v1)
        #undef VEC_INSERT_VAL(__v1, __val, __pos)
	#undef VEC_CVT_128_256(__v1)	
	#undef VEC_SET1_VAL(__val)
	#undef VEC_POPCVT_CHAR(__ch)
	#undef VEC_LDPOPCVT_CHAR(__addr)
	#undef VEC_CMP_EQ(__v1, __v2)
	#undef VEC_SET_LSE(__val)
	#undef SHIFT_HAP(__v1, __val)
	#undef print256b(__v1)
        #undef MASK_VEC
        #undef VEC_SSE_TO_AVX(__vsLow, __vsHigh, __vdst)
        #undef VEC_SHIFT_LEFT_1BIT(__vs)
        #undef MASK_ALL_ONES
        #undef COMPARE_VECS(__v1, __v2)
        #undef _256_INT_TYPE
#endif

#define PRECISION s

#define MAIN_TYPE float
#define MAIN_TYPE_SIZE 32
#define UNION_TYPE mix_F
#define IF_128 IF_128f
#define IF_MAIN_TYPE IF_32
#define SHIFT_CONST1 12
#define SHIFT_CONST2 3
#define SHIFT_CONST3 4
#define _128_TYPE __m128
#define _256_TYPE __m256
#define _256_INT_TYPE __m256i
#define AVX_LENGTH 8
#define MAVX_COUNT  (MROWS+7)/AVX_LENGTH
#define HAP_TYPE UNION_TYPE
#define MASK_TYPE uint32_t
#define MASK_ALL_ONES 0xFFFFFFFF
#define MASK_VEC MaskVec_F

#define SET_VEC_ZERO(__vec)                     \
  __vec= _mm256_setzero_ps()

#define VEC_OR(__v1, __v2)                      \
  _mm256_or_ps(__v1, __v2)

#define VEC_ADD(__v1, __v2)                     \
  _mm256_add_ps(__v1, __v2)

#define VEC_SUB(__v1, __v2)                     \
  _mm256_sub_ps(__v1, __v2)

#define VEC_MUL(__v1, __v2)                     \
  _mm256_mul_ps(__v1, __v2)

#define VEC_DIV(__v1, __v2)			\
  _mm256_div_ps(__v1, __v2)

#define VEC_BLEND(__v1, __v2, __mask)           \
  _mm256_blend_ps(__v1, __v2, __mask)

#define VEC_BLENDV(__v1, __v2, __maskV)         \
  _mm256_blendv_ps(__v1, __v2, __maskV)

#define VEC_CAST_256_128(__v1)			\
  _mm256_castps256_ps128 (__v1)

#define VEC_EXTRACT_128(__v1, __im)             \
  _mm256_extractf128_ps (__v1, __im)

#define VEC_EXTRACT_UNIT(__v1, __im)		\
  _mm_extract_epi32(__v1, __im)

#define VEC_SET1_VAL128(__val)			\
  _mm_set1_ps(__val)

#define VEC_MOVE(__v1, __val)			\
  _mm_move_ss(__v1, __val)

#define VEC_CAST_128_256(__v1)                  \
  _mm256_castps128_ps256(__v1)

#define VEC_INSERT_VAL(__v1, __val, __pos)      \
  _mm256_insertf128_ps(__v1, __val, __pos)

#define VEC_CVT_128_256(__v1)                   \
  _mm256_cvtepi32_ps(__v1.i)

#define VEC_SET1_VAL(__val)                 	\
  _mm256_set1_ps(__val)

#define VEC_POPCVT_CHAR(__ch)                   \
  _mm256_cvtepi32_ps(_mm256_set1_epi32(__ch))

#define VEC_LDPOPCVT_CHAR(__addr)		\
   _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i const *)__addr))

#define VEC_CMP_EQ(__v1, __v2)                  \
  _mm256_cmp_ps(__v1, __v2, _CMP_EQ_OQ)

#define VEC_SET_LSE(__val)			\
  _mm256_set_ps(zero, zero, zero, zero, zero, zero, zero, __val);

#define SHIFT_HAP(__v1, __val)			\
  _vector_shift_lastavxs(__v1, __val.f);

#define print256b(__v1)                         \
  print256bFP(__v1)
  
#define VEC_SSE_TO_AVX(__vsLow, __vsHigh, __vdst)	\
  __vdst = _mm256_castps128_ps256(__vsLow) ;		\
  __vdst = _mm256_insertf128_ps(__vdst, __vsHigh, 1) ;

#define VEC_SHIFT_LEFT_1BIT(__vs)		\
  __vs = _mm_slli_epi32(__vs, 1) 

#define COMPARE_VECS(__v1, __v2, __first, __last) {			\
    float* ptr1 = (float*) (&__v1) ;					\
    float* ptr2 = (float*) (&__v2) ;					\
    for (int ei=__first; ei <= __last; ++ei) {				\
      if (ptr1[ei] != ptr2[ei]) {					\
	std::cout << "Float Mismatch at " << ei << ": "			\
		  << ptr1[ei] << " vs. " << ptr2[ei] << std::endl ;	\
	exit(0) ;							\
      }									\
    }									\
  }

