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


inline void CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn, MAIN_TYPE &shiftOut);
inline void CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn);
inline void CONCAT(CONCAT(precompute_masks_,SIMD_ENGINE), PRECISION)(const testcase& tc, int COLS, int numMaskVecs, MASK_TYPE (*maskArr)[NUM_DISTINCT_CHARS]);
inline void CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(const testcase& tc, char* rsArr, MASK_TYPE* lastMaskShiftOut, int beginRowIndex, int numRowsToProcess);
inline void CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE), PRECISION)(int maskIndex, MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, MASK_TYPE (*maskArr) [NUM_DISTINCT_CHARS], char* rsArr, MASK_TYPE* lastMaskShiftOut, MASK_TYPE maskBitCnt);
inline void CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, SIMD_TYPE& distm, SIMD_TYPE& _1_distm, SIMD_TYPE& distmChosen);
template<class NUMBER> inline void CONCAT(CONCAT(initializeVectors,SIMD_ENGINE), PRECISION)(int ROWS, int COLS, NUMBER* shiftOutM, NUMBER *shiftOutX, NUMBER *shiftOutY, Context<NUMBER> ctx, testcase *tc,  SIMD_TYPE *p_MM, SIMD_TYPE *p_GAPM, SIMD_TYPE *p_MX, SIMD_TYPE *p_XX, SIMD_TYPE *p_MY, SIMD_TYPE *p_YY, SIMD_TYPE *distm1D);
template<class NUMBER> inline void CONCAT(CONCAT(stripINITIALIZATION,SIMD_ENGINE), PRECISION)(
    int stripIdx, Context<NUMBER> ctx, testcase *tc, SIMD_TYPE &pGAPM, SIMD_TYPE &pMM, SIMD_TYPE &pMX, SIMD_TYPE &pXX, SIMD_TYPE &pMY, SIMD_TYPE &pYY,
    SIMD_TYPE &rs, UNION_TYPE &rsN, SIMD_TYPE &distm, SIMD_TYPE &_1_distm,  SIMD_TYPE *distm1D, SIMD_TYPE N_packed256, SIMD_TYPE *p_MM , SIMD_TYPE *p_GAPM , 
    SIMD_TYPE *p_MX, SIMD_TYPE *p_XX , SIMD_TYPE *p_MY, SIMD_TYPE *p_YY, UNION_TYPE &M_t_2, UNION_TYPE &X_t_2, UNION_TYPE &M_t_1, UNION_TYPE &X_t_1, 
    UNION_TYPE &Y_t_2, UNION_TYPE &Y_t_1, UNION_TYPE &M_t_1_y, NUMBER* shiftOutX, NUMBER* shiftOutM);
inline SIMD_TYPE CONCAT(CONCAT(computeDISTM,SIMD_ENGINE), PRECISION)(int d, int COLS, testcase * tc, HAP_TYPE &hap, SIMD_TYPE rs, UNION_TYPE rsN, SIMD_TYPE N_packed256, 
    SIMD_TYPE distm, SIMD_TYPE _1_distm);
inline void CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(UNION_TYPE &M_t, UNION_TYPE &X_t, UNION_TYPE &Y_t, UNION_TYPE &M_t_y, 
    UNION_TYPE M_t_2, UNION_TYPE X_t_2, UNION_TYPE Y_t_2, UNION_TYPE M_t_1, UNION_TYPE X_t_1, UNION_TYPE M_t_1_y, UNION_TYPE Y_t_1, 
    SIMD_TYPE pMM, SIMD_TYPE pGAPM, SIMD_TYPE pMX, SIMD_TYPE pXX, SIMD_TYPE pMY, SIMD_TYPE pYY, SIMD_TYPE distmSel);
template<class NUMBER> NUMBER CONCAT(CONCAT(compute_full_prob_,SIMD_ENGINE), PRECISION) (testcase *tc, NUMBER *before_last_log = NULL);

