#ifdef PRECISION

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

/*
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
*/

void GEN_INTRINSIC(precompute_masks_, PRECISION)(const testcase& tc, int COLS, int numMaskVecs, MASK_TYPE (*maskArr)[NUM_DISTINCT_CHARS]) {
  
  const int maskBitCnt = MAIN_TYPE_SIZE ;

  for (int vi=0; vi < numMaskVecs; ++vi) {
    for (int rs=0; rs < NUM_DISTINCT_CHARS; ++rs) {
      maskArr[vi][rs] = 0 ;
    }
    maskArr[vi][AMBIG_CHAR] = MASK_ALL_ONES ;
  }
 
  for (int col=1; col < COLS; ++col) {
    int mIndex = (col-1) / maskBitCnt ;
    int mOffset = (col-1) % maskBitCnt ;
    MASK_TYPE bitMask = ((MASK_TYPE)0x1) << (maskBitCnt-1-mOffset) ;

    char hapChar = ConvertChar::get(tc.hap[col-1]);

    if (hapChar == AMBIG_CHAR) {
      for (int ci=0; ci < NUM_DISTINCT_CHARS; ++ci) 
	maskArr[mIndex][ci] |= bitMask ;
    } 

    maskArr[mIndex][hapChar] |= bitMask ;
    // bit corresponding to col 1 will be the MSB of the mask 0
    // bit corresponding to col 2 will be the MSB-1 of the mask 0
    // ...
    // bit corresponding to col 32 will be the LSB of the mask 0
    // bit corresponding to col 33 will be the MSB of the mask 1
    // ...
  }

}

void GEN_INTRINSIC(init_masks_for_row_, PRECISION)(const testcase& tc, char* rsArr, MASK_TYPE* lastMaskShiftOut, int beginRowIndex, int numRowsToProcess) {

  for (int ri=0; ri < numRowsToProcess; ++ri) {
    rsArr[ri] = ConvertChar::get(tc.rs[ri+beginRowIndex-1]) ;
  }

  for (int ei=0; ei < AVX_LENGTH; ++ei) {
    lastMaskShiftOut[ei] = 0 ;
  }
}

#define SET_MASK_WORD(__dstMask, __srcMask, __lastShiftOut, __shiftBy, __maskBitCnt){ \
  MASK_TYPE __bitMask = (((MASK_TYPE)0x1) << __shiftBy) - 1 ;			\
  MASK_TYPE __nextShiftOut = (__srcMask & __bitMask) << (__maskBitCnt - __shiftBy) ; \
  __dstMask = (__srcMask >> __shiftBy) | __lastShiftOut ;		\
  __lastShiftOut = __nextShiftOut ;					\
}


void GEN_INTRINSIC(update_masks_for_cols_, PRECISION)(int maskIndex, MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, MASK_TYPE (*maskArr) [NUM_DISTINCT_CHARS], char* rsArr, MASK_TYPE* lastMaskShiftOut, MASK_TYPE maskBitCnt) {

  for (int ei=0; ei < AVX_LENGTH/2; ++ei) {
    SET_MASK_WORD(currMaskVecLow.masks[ei], maskArr[maskIndex][rsArr[ei]], 
	lastMaskShiftOut[ei], ei, maskBitCnt) ;

    int ei2 = ei + AVX_LENGTH/2 ; // the second entry index
    SET_MASK_WORD(currMaskVecHigh.masks[ei], maskArr[maskIndex][rsArr[ei2]], 
	lastMaskShiftOut[ei2], ei2, maskBitCnt) ;
  }

}

//void GEN_INTRINSIC(computeDistVec, PRECISION) (MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, _256_TYPE& distm, _256_TYPE& _1_distm, _256_TYPE& distmChosen, const _256_TYPE& distmSel, int firstRowIndex, int lastRowIndex) {

inline void GEN_INTRINSIC(computeDistVec, PRECISION) (MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, _256_TYPE& distm, _256_TYPE& _1_distm, _256_TYPE& distmChosen) {
  //#define computeDistVec() {					      

  _256_TYPE maskV ;
  VEC_SSE_TO_AVX(currMaskVecLow.vecf, currMaskVecHigh.vecf, maskV) ;

  distmChosen = VEC_BLENDV(distm, _1_distm, maskV) ;

  /*COMPARE_VECS(distmChosen, distmSel, firstRowIndex, lastRowIndex) ;*/

  VEC_SHIFT_LEFT_1BIT(currMaskVecLow.vec) ;
  VEC_SHIFT_LEFT_1BIT(currMaskVecHigh.vec) ;
}								     


/*
template<class NUMBER>
struct HmmData {
  int ROWS ;
  int COLS ;
  
  NUMBER shiftOutM[MROWS+MCOLS+AVX_LENGTH], shiftOutX[MROWS+MCOLS+AVX_LENGTH], shiftOutY[MROWS+MCOLS+AVX_LENGTH] ;
  Context<NUMBER> ctx ;
  testcase* tc ;
  _256_TYPE p_MM[MAVX_COUNT], p_GAPM[MAVX_COUNT], p_MX[MAVX_COUNT], p_XX[MAVX_COUNT], p_MY[MAVX_COUNT], p_YY[MAVX_COUNT], distm1D[MAVX_COUNT] ;
  _256_TYPE pGAPM, pMM, pMX, pXX, pMY, pYY ;

  UNION_TYPE M_t, M_t_1, M_t_2, X_t, X_t_1, X_t_2, Y_t, Y_t_1, Y_t_2, M_t_y, M_t_1_y ;
  UNION_TYPE rs , rsN ;
  _256_TYPE distmSel;
  _256_TYPE distm, _1_distm;

} ;
*/

template<class NUMBER> void GEN_INTRINSIC(initializeVectors, PRECISION)(int ROWS, int COLS, NUMBER* shiftOutM, NUMBER *shiftOutX, NUMBER *shiftOutY, Context<NUMBER> ctx, testcase *tc,  _256_TYPE *p_MM, _256_TYPE *p_GAPM, _256_TYPE *p_MX, _256_TYPE *p_XX, _256_TYPE *p_MY, _256_TYPE *p_YY, _256_TYPE *distm1D)
{
	NUMBER zero = ctx._(0.0);
        NUMBER init_Y = ctx.INITIAL_CONSTANT / (tc->haplen);
        for (int s=0;s<ROWS+COLS+AVX_LENGTH;s++)
        {
                shiftOutM[s] = zero;
                shiftOutX[s] = zero;
                shiftOutY[s] = init_Y;
        }

        NUMBER *ptr_p_MM = (NUMBER *)p_MM;
        NUMBER *ptr_p_XX = (NUMBER *)p_XX;
        NUMBER *ptr_p_YY = (NUMBER *)p_YY;
        NUMBER *ptr_p_MX = (NUMBER *)p_MX;
        NUMBER *ptr_p_MY = (NUMBER *)p_MY;
        NUMBER *ptr_p_GAPM = (NUMBER *)p_GAPM;

        *ptr_p_MM = ctx._(0.0);
        *ptr_p_XX = ctx._(0.0);
        *ptr_p_YY = ctx._(0.0);
        *ptr_p_MX = ctx._(0.0);
        *ptr_p_MY = ctx._(0.0);
        *ptr_p_GAPM = ctx._(0.0);


        for (int r = 1; r < ROWS; r++)
        {
                int _i = tc->i[r-1] & 127;
                int _d = tc->d[r-1] & 127;
                int _c = tc->c[r-1] & 127;

                *(ptr_p_MM+r-1) = ctx._(1.0) - ctx.ph2pr[(_i + _d) & 127];
                *(ptr_p_GAPM+r-1) = ctx._(1.0) - ctx.ph2pr[_c];
                *(ptr_p_MX+r-1) = ctx.ph2pr[_i];
                *(ptr_p_XX+r-1) = ctx.ph2pr[_c];
                *(ptr_p_MY+r-1) = ctx.ph2pr[_d];
                *(ptr_p_YY+r-1) = ctx.ph2pr[_c];
#ifdef DEBUG3
		debug_dump("transitions_jni.txt",to_string(*(ptr_p_MM+r-1)  ),true);
		debug_dump("transitions_jni.txt",to_string(*(ptr_p_GAPM+r-1)),true);
		debug_dump("transitions_jni.txt",to_string(*(ptr_p_MX+r-1)  ),true);
		debug_dump("transitions_jni.txt",to_string(*(ptr_p_XX+r-1)  ),true);
		debug_dump("transitions_jni.txt",to_string(*(ptr_p_MY+r-1)  ),true);
		debug_dump("transitions_jni.txt",to_string(*(ptr_p_YY+r-1)  ),true);
#endif
		//*(ptr_p_MY+r-1) = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_d];
		//*(ptr_p_YY+r-1) = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_c];
        }

        NUMBER *ptr_distm1D = (NUMBER *)distm1D;
        for (int r = 1; r < ROWS; r++)
        {
                int _q = tc->q[r-1] & 127;
                ptr_distm1D[r-1] = ctx.ph2pr[_q];
#ifdef DEBUG3
		debug_dump("priors_jni.txt",to_string(ptr_distm1D[r-1]),true);
#endif
        }
}


template<class NUMBER> inline void GEN_INTRINSIC(stripINITIALIZATION, PRECISION)(
			int stripIdx, Context<NUMBER> ctx, testcase *tc, _256_TYPE &pGAPM, _256_TYPE &pMM, _256_TYPE &pMX, _256_TYPE &pXX, _256_TYPE &pMY, _256_TYPE &pYY,
			_256_TYPE &rs, UNION_TYPE &rsN, _256_TYPE &distm, _256_TYPE &_1_distm,  _256_TYPE *distm1D, _256_TYPE N_packed256, _256_TYPE *p_MM , _256_TYPE *p_GAPM , 
			_256_TYPE *p_MX, _256_TYPE *p_XX , _256_TYPE *p_MY, _256_TYPE *p_YY, UNION_TYPE &M_t_2, UNION_TYPE &X_t_2, UNION_TYPE &M_t_1, UNION_TYPE &X_t_1, 
			UNION_TYPE &Y_t_2, UNION_TYPE &Y_t_1, UNION_TYPE &M_t_1_y, NUMBER* shiftOutX, NUMBER* shiftOutM)	 
{
	int i = stripIdx;
        pGAPM = p_GAPM[i];                                               
        pMM   = p_MM[i];                                                 
        pMX   = p_MX[i];                                                 
        pXX   = p_XX[i];                                                 
        pMY   = p_MY[i];                                                 
        pYY   = p_YY[i];               

	NUMBER zero = ctx._(0.0);        
	NUMBER init_Y = ctx.INITIAL_CONSTANT / (tc->haplen);
	UNION_TYPE packed1;  packed1.d = VEC_SET1_VAL(1.0);
#define TRISTATE_CORRECTION_FACTOR 3.0
	UNION_TYPE packed3;  packed3.d = VEC_SET1_VAL(TRISTATE_CORRECTION_FACTOR);
	
        /* compare rs and N */                                           
        //rs = VEC_LDPOPCVT_CHAR((tc->irs+i*AVX_LENGTH));        
        //rsN.d = VEC_CMP_EQ(N_packed256, rs);                             
                                                                         
        distm = distm1D[i];                                    
        _1_distm = VEC_SUB(packed1.d, distm);
#ifndef DO_NOT_USE_TRISTATE_CORRECTION
	distm = VEC_DIV(distm, packed3.d); 
#endif
                                                                         
        /* initialize M_t_2, M_t_1, X_t_2, X_t_1, Y_t_2, Y_t_1 */        
        M_t_2.d = VEC_SET1_VAL(zero);                                    
        X_t_2.d = VEC_SET1_VAL(zero);                                    
                                                                         
        if (i==0) {                                                      
                M_t_1.d = VEC_SET1_VAL(zero);                            
                X_t_1.d = VEC_SET1_VAL(zero);                            
                Y_t_2.d = VEC_SET_LSE(init_Y);                           
                Y_t_1.d = VEC_SET1_VAL(zero);                            
        }                                                                
        else {                                                           
                X_t_1.d = VEC_SET_LSE(shiftOutX[AVX_LENGTH]);            
                M_t_1.d = VEC_SET_LSE(shiftOutM[AVX_LENGTH]);            
                Y_t_2.d = VEC_SET1_VAL(zero);                            
                Y_t_1.d = VEC_SET1_VAL(zero);                            
        }                                                                
        M_t_1_y = M_t_1;
}        



inline _256_TYPE GEN_INTRINSIC(computeDISTM, PRECISION)(int d, int COLS, testcase * tc, HAP_TYPE &hap, _256_TYPE rs, UNION_TYPE rsN, _256_TYPE N_packed256, 
						_256_TYPE distm, _256_TYPE _1_distm)
{
        UNION_TYPE hapN, rshap;
        _256_TYPE  cond;
	IF_32 shiftInHap;

	int *hap_ptr = tc->ihap;

        shiftInHap.i = (d<COLS) ? hap_ptr[d-1] : hap_ptr[COLS-1];              

        /* shift hap */                                                        
        SHIFT_HAP(hap, shiftInHap);                                            
        _256_TYPE hapF = VEC_CVT_128_256(hap);                                 

        rshap.d = VEC_CMP_EQ(rs, hapF);                                        
        hapN.d  = VEC_CMP_EQ(N_packed256, hapF);                               

        /* OR rsN, rshap, hapN */                                              
        cond =  VEC_OR(rsN.d, rshap.d);                                        
        cond =  VEC_OR(cond, hapN.d);                                          
        
        /* distm1D = (cond) ? 1-distm1D : distm1D;  */                         
        _256_TYPE distmSel = VEC_BLENDV(distm, _1_distm, cond);
        
        return distmSel;
}        


inline void GEN_INTRINSIC(computeMXY, PRECISION)(UNION_TYPE &M_t, UNION_TYPE &X_t, UNION_TYPE &Y_t, UNION_TYPE &M_t_y, 
			UNION_TYPE M_t_2, UNION_TYPE X_t_2, UNION_TYPE Y_t_2, UNION_TYPE M_t_1, UNION_TYPE X_t_1, UNION_TYPE M_t_1_y, UNION_TYPE Y_t_1, 
			_256_TYPE pMM, _256_TYPE pGAPM, _256_TYPE pMX, _256_TYPE pXX, _256_TYPE pMY, _256_TYPE pYY, _256_TYPE distmSel)
{
	/* Compute M_t <= distm * (p_MM*M_t_2 + p_GAPM*X_t_2 + p_GAPM*Y_t_2) */
	M_t.d = VEC_MUL(VEC_ADD(VEC_ADD(VEC_MUL(M_t_2.d, pMM), VEC_MUL(X_t_2.d, pGAPM)), VEC_MUL(Y_t_2.d, pGAPM)), distmSel);
	M_t_y = M_t;

	/* Compute X_t */
	X_t.d = VEC_ADD(VEC_MUL(M_t_1.d, pMX) , VEC_MUL(X_t_1.d, pXX));

	/* Compute Y_t */
	Y_t.d = VEC_ADD(VEC_MUL(M_t_1_y.d, pMY) , VEC_MUL(Y_t_1.d, pYY));
}

template<class NUMBER> NUMBER GEN_INTRINSIC(compute_full_prob_avx, PRECISION) (testcase *tc, NUMBER *before_last_log = NULL)
{
        _256_TYPE p_MM   [MAVX_COUNT], p_GAPM [MAVX_COUNT], p_MX   [MAVX_COUNT];
	_256_TYPE p_XX   [MAVX_COUNT], p_MY   [MAVX_COUNT], p_YY   [MAVX_COUNT];
	_256_TYPE distm1D[MAVX_COUNT];
	NUMBER shiftOutM[MROWS+MCOLS+AVX_LENGTH], shiftOutX[MROWS+MCOLS+AVX_LENGTH], shiftOutY[MROWS+MCOLS+AVX_LENGTH];
	UNION_TYPE  M_t, M_t_1, M_t_2, X_t, X_t_1, X_t_2, Y_t, Y_t_1, Y_t_2, M_t_y, M_t_1_y;
	_256_TYPE pGAPM, pMM, pMX, pXX, pMY, pYY;

        struct timeval start, end;
        NUMBER result_avx2;
        Context<NUMBER> ctx;
        UNION_TYPE rs , rsN;
        HAP_TYPE hap;
        _256_TYPE distmSel, distmChosen ;
	_256_TYPE distm, _1_distm;

        int r, c;
        int ROWS = tc->rslen + 1;
        int COLS = tc->haplen + 1;
        int AVX_COUNT = (ROWS+7)/8;
        NUMBER zero = ctx._(0.0);
        UNION_TYPE packed1;  packed1.d = VEC_SET1_VAL(1.0);
	_256_TYPE N_packed256 = VEC_POPCVT_CHAR('N');
        NUMBER init_Y = ctx.INITIAL_CONSTANT / (tc->haplen);
        int remainingRows = (ROWS-1) % AVX_LENGTH;
        int strip_cnt = ((ROWS-1) / AVX_LENGTH) + (remainingRows!=0);

	const int maskBitCnt = MAIN_TYPE_SIZE ;
	const int numMaskVecs = (COLS+ROWS+maskBitCnt-1)/maskBitCnt ; // ceil function

	MASK_TYPE maskArr[numMaskVecs][NUM_DISTINCT_CHARS] ;
	GEN_INTRINSIC(precompute_masks_, PRECISION)(*tc, COLS, numMaskVecs, maskArr) ;

	char rsArr[AVX_LENGTH] ;
	MASK_TYPE lastMaskShiftOut[AVX_LENGTH] ;

	GEN_INTRINSIC(initializeVectors, PRECISION)<NUMBER>(ROWS, COLS, shiftOutM, shiftOutX, shiftOutY, 
							    ctx, tc, p_MM, p_GAPM, p_MX, p_XX, p_MY, p_YY, distm1D);

	//for (int __ii=0; __ii < 10; ++__ii)
        for (int i=0;i<strip_cnt-1;i++)
        {
		//STRIP_INITIALIZATION
		GEN_INTRINSIC(stripINITIALIZATION, PRECISION)(i, ctx, tc, pGAPM, pMM, pMX, pXX, pMY, pYY, rs.d, rsN, distm, _1_distm, distm1D, N_packed256, p_MM , p_GAPM ,
							      p_MX, p_XX , p_MY, p_YY, M_t_2, X_t_2, M_t_1, X_t_1, Y_t_2, Y_t_1, M_t_1_y, shiftOutX, shiftOutM);

		GEN_INTRINSIC(init_masks_for_row_, PRECISION)(*tc, rsArr, lastMaskShiftOut,
							      i*AVX_LENGTH+1, AVX_LENGTH) ;
		// Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors
		MASK_VEC currMaskVecLow ; // corresponding to lower half
		MASK_VEC currMaskVecHigh ; // corresponding to upper half

                for (int d=1;d<COLS+AVX_LENGTH;d++)
                {
		        if (d % MAIN_TYPE_SIZE == 1)
			  GEN_INTRINSIC(update_masks_for_cols_, PRECISION)((d-1)/MAIN_TYPE_SIZE, currMaskVecLow, currMaskVecHigh, maskArr, rsArr, lastMaskShiftOut, maskBitCnt) ;

			int ShiftIdx = d+AVX_LENGTH;
			//distmSel = GEN_INTRINSIC(computeDISTM, PRECISION)(d, COLS, tc, hap, rs.d, rsN, N_packed256, distm, _1_distm);

			//int firstRowIndex = (d < COLS) ? 0 : (d-COLS+1) ;
			//int lastRowIndex = std::min(d-1, AVX_LENGTH-1) ;
			//GEN_INTRINSIC(computeDistVec, PRECISION) (currMaskVecLow, currMaskVecHigh, distm, _1_distm, distmChosen, distmSel, firstRowIndex, lastRowIndex) ;

			GEN_INTRINSIC(computeDistVec, PRECISION) (currMaskVecLow, currMaskVecHigh, distm, _1_distm, distmChosen) ;
			

			GEN_INTRINSIC(computeMXY, PRECISION)(M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1, 
								pMM, pGAPM, pMX, pXX, pMY, pYY, distmChosen);

                        GEN_INTRINSIC(_vector_shift, PRECISION)(M_t, shiftOutM[ShiftIdx], shiftOutM[d]);

                        GEN_INTRINSIC(_vector_shift, PRECISION)(X_t, shiftOutX[ShiftIdx], shiftOutX[d]);

                        GEN_INTRINSIC(_vector_shift, PRECISION)(Y_t_1, shiftOutY[ShiftIdx], shiftOutY[d]);

                        M_t_2 = M_t_1; M_t_1 = M_t; X_t_2 = X_t_1; X_t_1 = X_t;
                        Y_t_2 = Y_t_1; Y_t_1 = Y_t; M_t_1_y = M_t_y;

		}
	}

        int i = strip_cnt-1;
        {
		//STRIP_INITIALIZATION
		GEN_INTRINSIC(stripINITIALIZATION, PRECISION)(i, ctx, tc, pGAPM, pMM, pMX, pXX, pMY, pYY, rs.d, rsN, distm, _1_distm, distm1D, N_packed256, p_MM , p_GAPM ,
                        p_MX, p_XX , p_MY, p_YY, M_t_2, X_t_2, M_t_1, X_t_1, Y_t_2, Y_t_1, M_t_1_y, shiftOutX, shiftOutM);

                if (remainingRows==0) remainingRows=AVX_LENGTH;

		GEN_INTRINSIC(init_masks_for_row_, PRECISION)(*tc, rsArr, lastMaskShiftOut,
							      i*AVX_LENGTH+1, remainingRows) ;

                _256_TYPE sumM, sumX;
                sumM = VEC_SET1_VAL(zero);
                sumX = VEC_SET1_VAL(zero);

		// Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors
		MASK_VEC currMaskVecLow ; // corresponding to lower half
		MASK_VEC currMaskVecHigh ; // corresponding to upper half

                for (int d=1;d<COLS+remainingRows-1;d++)
                {

  	 	        if (d % MAIN_TYPE_SIZE == 1)
			  GEN_INTRINSIC(update_masks_for_cols_, PRECISION)((d-1)/MAIN_TYPE_SIZE, currMaskVecLow, currMaskVecHigh, maskArr, rsArr, lastMaskShiftOut, maskBitCnt) ;

			int ShiftIdx = d+AVX_LENGTH;
			//distmSel = GEN_INTRINSIC(computeDISTM, PRECISION)(d, COLS, tc, hap, rs.d, rsN, N_packed256, distm, _1_distm);
			GEN_INTRINSIC(computeDistVec, PRECISION) (currMaskVecLow, currMaskVecHigh, distm, _1_distm, distmChosen) ;

			GEN_INTRINSIC(computeMXY, PRECISION)(M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
					                        pMM, pGAPM, pMX, pXX, pMY, pYY, distmChosen);

                        sumM  = VEC_ADD(sumM, M_t.d);
                        GEN_INTRINSIC(_vector_shift_last, PRECISION)(M_t, shiftOutM[ShiftIdx]);

                        sumX  = VEC_ADD(sumX, X_t.d);
                        GEN_INTRINSIC(_vector_shift_last, PRECISION)(X_t, shiftOutX[ShiftIdx]);

                        GEN_INTRINSIC(_vector_shift_last, PRECISION)(Y_t_1, shiftOutY[ShiftIdx]);

                        M_t_2 = M_t_1; M_t_1 = M_t; X_t_2 = X_t_1; X_t_1 = X_t;
                        Y_t_2 = Y_t_1; Y_t_1 = Y_t; M_t_1_y = M_t_y;

                }
                UNION_TYPE sumMX;
                sumMX.d = VEC_ADD(sumM, sumX);
                result_avx2 = sumMX.f[remainingRows-1];		
	}
	//printf("result_avx2: %f\n", result_avx2);
	return result_avx2;
}

#endif

