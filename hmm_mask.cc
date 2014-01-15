#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


// /usr/intel/pkgs/icc/13.0.0e/bin/icc -o ed -O3 ed.cpp -xAVX  -openmp -openmp-link static


#include <iostream>

#include <malloc.h>
#include <assert.h>
#include <sys/time.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>


#include <immintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>


#include <stdint.h>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <memory.h>
#include <string>
#include <sstream>
#include <fstream>

using namespace std ;

#define CHECK_MASK_CORRECTNESS

template <class T>
string getBinaryStr (T val, int numBitsToWrite) {
  
  ostringstream oss ;
  uint64_t mask = ((T) 0x1) << (numBitsToWrite-1) ;
  for (int i=numBitsToWrite-1; i >= 0; --i) {
    oss << ((val & mask) >> i) ;
    mask >>= 1 ;
  }
  return oss.str() ;
}


typedef struct 
{
	int rslen, haplen, *q, *i, *d, *c;
	char *hap, *rs;
} testcase;

int normalize(char c)
{
	return ((int) (c - 33));
}

class ConvertChar {

  static uint8_t conversionTable[255] ;


public:

  static void init() {
    conversionTable['A'] = 0 ;
    conversionTable['C'] = 1 ;
    conversionTable['T'] = 2 ;
    conversionTable['G'] = 3 ;
    conversionTable['N'] = 4 ;
  }

  static inline uint8_t get(uint8_t input) {
    return conversionTable[input] ;
  }

} ;

uint8_t ConvertChar::conversionTable[255] ;



int read_testcase(testcase *tc, FILE* ifp)
{
	char *q, *i, *d, *c, *line = NULL;
	int _q, _i, _d, _c;
	int x, size = 0;
	ssize_t read;

	read = getline(&line, (size_t *) &size, ifp);
	if (read == -1)
		return -1;


	tc->hap = (char *) malloc(size);
	tc->rs = (char *) malloc(size);
	q = (char *) malloc(size);
	i = (char *) malloc(size);
	d = (char *) malloc(size);
	c = (char *) malloc(size);

	if (sscanf(line, "%s %s %s %s %s %s\n", tc->hap, tc->rs, q, i, d, c) != 6)
		return -1;

	tc->haplen = strlen(tc->hap);
	tc->rslen = strlen(tc->rs);
	tc->q = (int *) malloc(sizeof(int) * tc->rslen);
	tc->i = (int *) malloc(sizeof(int) * tc->rslen);
	tc->d = (int *) malloc(sizeof(int) * tc->rslen);
	tc->c = (int *) malloc(sizeof(int) * tc->rslen);

	// Convert hap and rs to 3 bits
	
	/*
	for (int ci=0; ci < tc->haplen; ++ci)
	  tc->hap[ci] = ConvertChar::get(tc->hap[ci]) ;

	for (int ci=0; ci < tc->rslen; ++ci)
	  tc->rs[ci] = ConvertChar::get(tc->rs[ci]) ;
	*/

	for (x = 0; x < tc->rslen; x++)
	{
		_q = normalize(q[x]);
		_i = normalize(i[x]);
		_d = normalize(d[x]);
		_c = normalize(c[x]);
		tc->q[x] = (_q < 6) ? 6 : _q;
		tc->i[x] = _i;
		tc->d[x] = _d;
		tc->c[x] = _c;
	}
	

	free(q);
	free(i);
	free(d);
	free(c);
	free(line);

	return 0;
}




#define ALIGN __attribute__ ((aligned(32)))


typedef union ALIGN { 
  //__m256i vi; 
  __m128 vf ;
  __m128i vi ;
  __m128d vd; 

  uint8_t b[16];
  //uint16_t hw[8] ;
  uint32_t w[4] ;
  uint64_t dw[2] ;
  
  float f[4] ;
  //double d[2] ;

} v128 ;

typedef union ALIGN { 
  __m256i vi; 
  __m256 vf ;
  __m256d vd; 
  //__m128i vi128[2] ;

  uint8_t b[32];
  //uint16_t hw[16] ;
  uint32_t w[8] ;
  uint64_t dw[4] ;
  
  float f[8] ;
  double d[4] ;

} v256 ;


#define NUM_DISTINCT_CHARS 5
#define AMBIG_CHAR 4

#define VEC_ENTRY_CNT 8
#define VEC_LEN 256
#define VTYPE vf
#define VTYPEI vi
#define VENTRY f
#define VENTRYI w
#define MTYPE uint32_t
#define MASK_ALL_ONES 0xFFFFFFFF


#define VECTOR v256
#define VECTOR_SSE v128

#define SET_VEC_ZERO(__vec)			\
  __vec= _mm256_setzero_ps()

#define SET_VEC_ONES(__vec)			\
  __vec = _mm256_set1_epi32(0xFFFFFFFF)

#define VEC_OR(__v1, __v2)			\
  _mm256_or_ps(__v1, __v2)

#define VEC_ADD(__v1, __v2)			\
  _mm256_add_ps(__v1, __v2)

#define VEC_MUL(__v1, __v2)			\
  _mm256_mul_ps(__v1, __v2) 

#define VEC_BLENDV(__v1, __v2, __maskV)		\
  _mm256_blendv_ps(__v1, __v2, __maskV)

#define VEC_SSE_TO_AVX(__vsLow, __vsHigh, __vdst)	\
  __vdst = _mm256_castps128_ps256(__vsLow) ;		\
  __vdst = _mm256_insertf128_ps(__vdst, __vsHigh, 1) ;

#define VECS_SHIFT_LEFT_1BIT(__vs)		\
  __vs = _mm_slli_epi32(__vs, 1) 


#define VECS_SHIFT_RIGHT_WHOLE(__vs, __cnt) {				\
    uint64_t __mask = ;							\
    uint64_t __shiftWord = ((((uint64_t) 0x1) << __cnt) - 1 ) & __vs.dw[1] ; \
    __shiftWord <<= (64-cnt) ;						\
    __vs = _mm_slri_epi64(__vs, __cnt) ;				\
    __vs.dw[0] |= __shiftWord ;						\
  }


#define GET_MASK_WORD_OLD(__mask, __lastShiftOut, __shiftBy, __maskBitCnt){ \
    MTYPE __bitMask = (((MTYPE)0x1) << __shiftBy) - 1 ;			\
    MTYPE __nextShiftOut = (__mask & __bitMask) << (__maskBitCnt - __shiftBy) ; \
    __mask >>= __shiftBy ;						\
    __mask |= __lastShiftOut ;						\
    __lastShiftOut = __nextShiftOut ;					\
  }


#define SET_MASK_WORD(__dstMask, __srcMask, __lastShiftOut, __shiftBy, __maskBitCnt){ \
  MTYPE __bitMask = (((MTYPE)0x1) << __shiftBy) - 1 ;			\
  MTYPE __nextShiftOut = (__srcMask & __bitMask) << (__maskBitCnt - __shiftBy) ; \
  __dstMask = (__srcMask >> __shiftBy) | __lastShiftOut ;		\
  __lastShiftOut = __nextShiftOut ;					\ 
}


void precompute_masks_avx(const testcase& tc, int COLS, int numMaskVecs,
			  MTYPE (*maskArr)[NUM_DISTINCT_CHARS]) {

  const int maskBitCnt = VEC_LEN / VEC_ENTRY_CNT ;

  for (int vi=0; vi < numMaskVecs; ++vi) {
    for (int rs=0; rs < NUM_DISTINCT_CHARS; ++rs) {
      maskArr[vi][rs] = 0 ;
    }
    maskArr[vi][AMBIG_CHAR] = MASK_ALL_ONES ;
  }
 
  for (int col=1; col < COLS; ++col) {
    int mIndex = (col-1) / maskBitCnt ;
    int mOffset = (col-1) % maskBitCnt ;
    MTYPE bitMask = ((MTYPE)0x1) << (maskBitCnt-1-mOffset) ;

    char hapChar = tc.hap[col-1] ;

    if (hapChar == AMBIG_CHAR) {
      maskArr[mIndex][0] |= bitMask ;
      maskArr[mIndex][1] |= bitMask ;
      maskArr[mIndex][2] |= bitMask ;
      maskArr[mIndex][3] |= bitMask ;
    } 

    //cout << hapChar << " " << mIndex << " " << getBinaryStr<MTYPE>(bitMask, 32)
    // << endl ;
    //cout << getBinaryStr<MTYPE>(maskArr[0][hapChar],32) << endl ;
    //exit(0) ;
		     
    maskArr[mIndex][ConvertChar::get(hapChar)] |= bitMask ;


    // bit corresponding to col 1 will be the MSB of the mask 0
    // bit corresponding to col 2 will be the MSB-1 of the mask 0
    // ...
    // bit corresponding to col 32 will be the LSB of the mask 0
    // bit corresponding to col 33 will be the MSB of the mask 1
    // ...
  }

}


void test_mask_computations (testcase& tc, int tcID, bool printDebug=false) {

  int ROWS = tc.rslen + 1 ;
  int COLS = tc.haplen + 1 ;

  // only for testing
  VECTOR mismatchData, matchData ;
  //SET_VEC_ZERO(mismatchData.VTYPE) ;
  //SET_VEC_ONES(matchData.VTYPEI) ;
  for (int ei=0; ei < VEC_ENTRY_CNT; ++ei) {
    matchData.VENTRY[ei] = 1.0 ;
    mismatchData.VENTRY[ei] = 0.0 ;
  }
  
  const int maskBitCnt = VEC_LEN / VEC_ENTRY_CNT ;
  const int numMaskVecs = (COLS+ROWS+maskBitCnt-1)/maskBitCnt ; // ceil function

  MTYPE maskArr[numMaskVecs][NUM_DISTINCT_CHARS] ;
  precompute_masks_avx(tc, COLS, numMaskVecs, maskArr) ;

#ifdef DEBUG 
  if (printDebug) {
    cout << "The first 32 hap chars are: " ;
    for (int i=0; i < 32; ++i) {
      cout << tc.hap[i] ;
    }
    cout << endl ;

    cout << "Masks computed for A, C, T, G, N are: "  << endl ;
    cout << getBinaryStr<MTYPE>(maskArr[0][0], 32)  << endl ;
    cout << getBinaryStr<MTYPE>(maskArr[0][1], 32)  << endl ;
    cout << getBinaryStr<MTYPE>(maskArr[0][2], 32)  << endl ;
    cout << getBinaryStr<MTYPE>(maskArr[0][3], 32)  << endl ;
    cout << getBinaryStr<MTYPE>(maskArr[0][4], 32)  << endl ;
  }
#endif // #ifdef DEBUG 

  int beginRowIndex = 1 ;
  while (beginRowIndex < ROWS) {

    int numRowsToProcess = min(VEC_ENTRY_CNT, ROWS - beginRowIndex) ;

    char rsArr[VEC_ENTRY_CNT] ;
    for (int ri=0; ri < numRowsToProcess; ++ri) {
      rsArr[ri] = ConvertChar::get(tc.rs[ri+beginRowIndex-1]) ;
    }

    // Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors
    VECTOR_SSE currMaskVecLow ; // corresponding to entries 0-3
    VECTOR_SSE currMaskVecHigh ; // corresponding to entries 4-7

    MTYPE lastMaskShiftOut[VEC_ENTRY_CNT] ;
    for (int ei=0; ei < VEC_ENTRY_CNT; ++ei) 
      lastMaskShiftOut[ei] = 0 ;

    int col = 1 ;
    int diag = 1 ;
    for (int maskIndex=0; maskIndex < numMaskVecs; ++maskIndex) {
      // set up the mask vectors for the next maskBitCnt columns

      // For AVX, maskBitCnt = 32 (so, the operation below is amortized over 32 cols)
      for (int ei=0; ei < VEC_ENTRY_CNT/2; ++ei) {
	SET_MASK_WORD(currMaskVecLow.VENTRYI[ei], maskArr[maskIndex][rsArr[ei]], 
		      lastMaskShiftOut[ei], ei, maskBitCnt) ;

      	int ei2 = ei + VEC_ENTRY_CNT/2 ; // the second entry index
	SET_MASK_WORD(currMaskVecHigh.VENTRYI[ei], maskArr[maskIndex][rsArr[ei2]], 
		      lastMaskShiftOut[ei2], ei2, maskBitCnt) ;
      }

#ifdef DEBUG 
      if (printDebug && maskIndex == 0) {
	cout << "The masks for entry 1: " << endl  
	     << getBinaryStr<MTYPE>(maskArr[0][rsArr[1]], 32) << endl
	     << getBinaryStr<MTYPE>(currMaskVecLow.VENTRYI[1], 32) << endl ;

      }
#endif // #ifdef DEBUG 
      
      // iterate over mask bit indices and columns
      for (int mbi=0; mbi < maskBitCnt && diag < COLS + ROWS -2; ++mbi, ++diag) {

	VECTOR maskV ;
	VEC_SSE_TO_AVX(currMaskVecLow.VTYPE, currMaskVecHigh.VTYPE, maskV.VTYPE) ;

	VECTOR testData ;
	testData.VTYPE = VEC_BLENDV(mismatchData.VTYPE, matchData.VTYPE,
				    maskV.VTYPE) ;

	VECS_SHIFT_LEFT_1BIT(currMaskVecLow.VTYPEI) ;
	VECS_SHIFT_LEFT_1BIT(currMaskVecHigh.VTYPEI) ;

#ifdef DEBUG 
	if (printDebug && maskIndex == 0) {
	  cout << "The mask for entry 1, mbi=" << mbi << ": "
	       << getBinaryStr<MTYPE>(maskV.VENTRYI[1], 32) << endl ;

	}
#endif // #ifdef DEBUG 

#ifdef CHECK_MASK_CORRECTNESS

	int firstRowIndex = (diag < COLS) ? 0 : (diag - COLS) ;
	int lastRowIndex = min(col-1, numRowsToProcess-1) ;

	for (int ri=firstRowIndex; ri <= lastRowIndex; ++ri) {
	  int currRow = beginRowIndex + ri ;
	  int currCol = col - ri + firstRowIndex ;

	  char hapChar = tc.hap[currCol-1] ;
	  char rsChar = tc.rs[currRow-1] ;

	  bool match = (hapChar == rsChar || hapChar == 'N' || rsChar == 'N') ;

	  if ((bool) testData.VENTRYI[ri] != match) {
	    cout << "Error: Incorrect mask for tc " << tcID << ", diag = " << diag
		 << " (" << currRow  << ", " << currCol << ")" << endl ;

	    cout << "The chars are: " << hapChar << " and " << rsChar << endl ;
	    cout << "The selected value is: " << testData.VENTRYI[ri] << endl ;

	    exit(0) ;
	  }
	}
#endif // #ifdef CHECK_MASK_CORRECTNESS

	if (diag < COLS)
	  ++col ;

      } // mbi
    } // maskIndex

    beginRowIndex += VEC_ENTRY_CNT ;
  } // end of stripe

  //cout << "Finished validating entry " << endl ;
}


int main () {

  #define BATCH_SIZE 10000

  ConvertChar::init() ;
  
  testcase* tcBatch = new testcase[BATCH_SIZE] ;
  
  int numBatches = 0 ;
  int numRead = 0 ;

  const int DEBUG_TC = -1 ;

  FILE* ifp = stdin ;

  do {
    numRead = 0 ;

    while (numRead < BATCH_SIZE) {
      if (read_testcase(tcBatch+numRead, ifp) < 0)
	break ;
      
      ++numRead ;
    }

    if (numRead == 0)
      break ;

    for (int ti=0; ti < numRead; ++ti) {

      int tcID = numBatches * BATCH_SIZE + ti ;
      test_mask_computations(tcBatch[ti], tcID, tcID == DEBUG_TC) ;
    }

    ++numBatches ;
  } while (numRead == BATCH_SIZE) ;
  
  fclose(ifp) ;

  delete[] tcBatch ;

  return 0 ;
}
