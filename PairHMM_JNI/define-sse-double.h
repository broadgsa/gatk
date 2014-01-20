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

	#undef VEC_EXTRACT_UNIT(__v1, __im)
	#undef VEC_INSERT_UNIT(__v1,__ins,__im)
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

#define SSE
#define PRECISION d

#define MAIN_TYPE double
#define MAIN_TYPE_SIZE 64
#define UNION_TYPE mix_D128
#define IF_128 IF_128d
#define IF_MAIN_TYPE IF_64
#define SHIFT_CONST1 1
#define SHIFT_CONST2 8
#define SHIFT_CONST3 0
#define _128_TYPE __m128d
#define _256_TYPE __m128d
#define _256_INT_TYPE __m128i
#define AVX_LENGTH 2
#define MAVX_COUNT  (MROWS+3)/AVX_LENGTH
#define HAP_TYPE __m128i
#define MASK_TYPE uint64_t
#define MASK_ALL_ONES 0xFFFFFFFFFFFFFFFFL
#define MASK_VEC MaskVec_D128

#define VEC_EXTRACT_UNIT(__v1, __im)            \
  _mm_extract_epi64(__v1, __im)

#define VEC_INSERT_UNIT(__v1,__ins,__im)	\
   _mm_insert_epi64(__v1,__ins,__im)

#define VEC_OR(__v1, __v2)                      \
  _mm_or_pd(__v1, __v2)

#define VEC_ADD(__v1, __v2)                     \
  _mm_add_pd(__v1, __v2)

#define VEC_SUB(__v1, __v2)                     \
  _mm_sub_pd(__v1, __v2)

#define VEC_MUL(__v1, __v2)                     \
  _mm_mul_pd(__v1, __v2)

#define VEC_DIV(__v1, __v2)                     \
  _mm_div_pd(__v1, __v2)

#define VEC_CMP_EQ(__v1, __v2)                  \
  _mm_cmpeq_pd(__v1, __v2)

#define VEC_BLEND(__v1, __v2, __mask)           \
  _mm_blend_pd(__v1, __v2, __mask)

#define VEC_BLENDV(__v1, __v2, __maskV)         \
  _mm_blendv_pd(__v1, __v2, __maskV)

#define SHIFT_HAP(__v1, __val)                  \
  __v1 = _mm_insert_epi32(_mm_slli_si128(__v1, 4), __val.i, 0)

#define VEC_CVT_128_256(__v1)                   \
  _mm_cvtepi32_pd(__v1)

#define VEC_SET1_VAL(__val)			\
  _mm_set1_pd(__val)
   
#define VEC_POPCVT_CHAR(__ch)                   \
  _mm_cvtepi32_pd(_mm_set1_epi32(__ch)) 

#define VEC_SET_LSE(__val)                      \
  _mm_set_pd(zero, __val);  

#define VEC_LDPOPCVT_CHAR(__addr)               \
   _mm_cvtepi32_pd(_mm_loadu_si128((__m128i const *)__addr)) 

#define VEC_SSE_TO_AVX(__vsLow, __vsHigh, __vdst)       \
   __vdst = _mm_castsi128_pd(_mm_set_epi64(__vsHigh, __vsLow))

#define VEC_SHIFT_LEFT_1BIT(__vs)               \
  __vs = _mm_slli_si64(__vs, 1)


