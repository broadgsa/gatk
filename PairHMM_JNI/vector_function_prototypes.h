inline void GEN_INTRINSIC(GEN_INTRINSIC(_vector_shift,SIMD_TYPE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn, MAIN_TYPE &shiftOut);
inline void GEN_INTRINSIC(GEN_INTRINSIC(_vector_shift_last,SIMD_TYPE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn);
inline void GEN_INTRINSIC(GEN_INTRINSIC(precompute_masks_,SIMD_TYPE), PRECISION)(const testcase& tc, int COLS, int numMaskVecs, MASK_TYPE (*maskArr)[NUM_DISTINCT_CHARS]);
inline void GEN_INTRINSIC(GEN_INTRINSIC(init_masks_for_row_,SIMD_TYPE), PRECISION)(const testcase& tc, char* rsArr, MASK_TYPE* lastMaskShiftOut, int beginRowIndex, int numRowsToProcess);
inline void GEN_INTRINSIC(GEN_INTRINSIC(update_masks_for_cols_,SIMD_TYPE), PRECISION)(int maskIndex, MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, MASK_TYPE (*maskArr) [NUM_DISTINCT_CHARS], char* rsArr, MASK_TYPE* lastMaskShiftOut, MASK_TYPE maskBitCnt);
inline void GEN_INTRINSIC(GEN_INTRINSIC(computeDistVec,SIMD_TYPE), PRECISION) (MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, _256_TYPE& distm, _256_TYPE& _1_distm, _256_TYPE& distmChosen);
template<class NUMBER> inline void GEN_INTRINSIC(GEN_INTRINSIC(initializeVectors,SIMD_TYPE), PRECISION)(int ROWS, int COLS, NUMBER* shiftOutM, NUMBER *shiftOutX, NUMBER *shiftOutY, Context<NUMBER> ctx, testcase *tc,  _256_TYPE *p_MM, _256_TYPE *p_GAPM, _256_TYPE *p_MX, _256_TYPE *p_XX, _256_TYPE *p_MY, _256_TYPE *p_YY, _256_TYPE *distm1D);
template<class NUMBER> inline void GEN_INTRINSIC(GEN_INTRINSIC(stripINITIALIZATION,SIMD_TYPE), PRECISION)(
    int stripIdx, Context<NUMBER> ctx, testcase *tc, _256_TYPE &pGAPM, _256_TYPE &pMM, _256_TYPE &pMX, _256_TYPE &pXX, _256_TYPE &pMY, _256_TYPE &pYY,
    _256_TYPE &rs, UNION_TYPE &rsN, _256_TYPE &distm, _256_TYPE &_1_distm,  _256_TYPE *distm1D, _256_TYPE N_packed256, _256_TYPE *p_MM , _256_TYPE *p_GAPM , 
    _256_TYPE *p_MX, _256_TYPE *p_XX , _256_TYPE *p_MY, _256_TYPE *p_YY, UNION_TYPE &M_t_2, UNION_TYPE &X_t_2, UNION_TYPE &M_t_1, UNION_TYPE &X_t_1, 
    UNION_TYPE &Y_t_2, UNION_TYPE &Y_t_1, UNION_TYPE &M_t_1_y, NUMBER* shiftOutX, NUMBER* shiftOutM);
inline _256_TYPE GEN_INTRINSIC(GEN_INTRINSIC(computeDISTM,SIMD_TYPE), PRECISION)(int d, int COLS, testcase * tc, HAP_TYPE &hap, _256_TYPE rs, UNION_TYPE rsN, _256_TYPE N_packed256, 
    _256_TYPE distm, _256_TYPE _1_distm);
inline void GEN_INTRINSIC(GEN_INTRINSIC(computeMXY,SIMD_TYPE), PRECISION)(UNION_TYPE &M_t, UNION_TYPE &X_t, UNION_TYPE &Y_t, UNION_TYPE &M_t_y, 
    UNION_TYPE M_t_2, UNION_TYPE X_t_2, UNION_TYPE Y_t_2, UNION_TYPE M_t_1, UNION_TYPE X_t_1, UNION_TYPE M_t_1_y, UNION_TYPE Y_t_1, 
    _256_TYPE pMM, _256_TYPE pGAPM, _256_TYPE pMX, _256_TYPE pXX, _256_TYPE pMY, _256_TYPE pYY, _256_TYPE distmSel);
template<class NUMBER> NUMBER GEN_INTRINSIC(GEN_INTRINSIC(compute_full_prob_,SIMD_TYPE), PRECISION) (testcase *tc, NUMBER *before_last_log = NULL);

