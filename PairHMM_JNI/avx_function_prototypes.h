void GEN_INTRINSIC(_vector_shift, PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn, MAIN_TYPE &shiftOut);
void GEN_INTRINSIC(_vector_shift_last, PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn);
void GEN_INTRINSIC(precompute_masks_, PRECISION)(const testcase& tc, int COLS, int numMaskVecs, MASK_TYPE (*maskArr)[NUM_DISTINCT_CHARS]);
void GEN_INTRINSIC(init_masks_for_row_, PRECISION)(const testcase& tc, char* rsArr, MASK_TYPE* lastMaskShiftOut, int beginRowIndex, int numRowsToProcess);
void GEN_INTRINSIC(update_masks_for_cols_, PRECISION)(int maskIndex, MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, MASK_TYPE (*maskArr) [NUM_DISTINCT_CHARS], char* rsArr, MASK_TYPE* lastMaskShiftOut, MASK_TYPE maskBitCnt);
void GEN_INTRINSIC(computeDistVec, PRECISION) (MASK_VEC& currMaskVecLow, MASK_VEC& currMaskVecHigh, _256_TYPE& distm, _256_TYPE& _1_distm, _256_TYPE& distmChosen);
template<class NUMBER> void GEN_INTRINSIC(initializeVectors, PRECISION)(int ROWS, int COLS, NUMBER* shiftOutM, NUMBER *shiftOutX, NUMBER *shiftOutY, Context<NUMBER> ctx, testcase *tc,  _256_TYPE *p_MM, _256_TYPE *p_GAPM, _256_TYPE *p_MX, _256_TYPE *p_XX, _256_TYPE *p_MY, _256_TYPE *p_YY, _256_TYPE *distm1D);
template<class NUMBER> void GEN_INTRINSIC(stripINITIALIZATION, PRECISION)(
    int stripIdx, Context<NUMBER> ctx, testcase *tc, _256_TYPE &pGAPM, _256_TYPE &pMM, _256_TYPE &pMX, _256_TYPE &pXX, _256_TYPE &pMY, _256_TYPE &pYY,
    _256_TYPE &rs, UNION_TYPE &rsN, _256_TYPE &distm, _256_TYPE &_1_distm,  _256_TYPE *distm1D, _256_TYPE N_packed256, _256_TYPE *p_MM , _256_TYPE *p_GAPM , 
    _256_TYPE *p_MX, _256_TYPE *p_XX , _256_TYPE *p_MY, _256_TYPE *p_YY, UNION_TYPE &M_t_2, UNION_TYPE &X_t_2, UNION_TYPE &M_t_1, UNION_TYPE &X_t_1, 
    UNION_TYPE &Y_t_2, UNION_TYPE &Y_t_1, UNION_TYPE &M_t_1_y, NUMBER* shiftOutX, NUMBER* shiftOutM);
_256_TYPE GEN_INTRINSIC(computeDISTM, PRECISION)(int d, int COLS, testcase * tc, HAP_TYPE &hap, _256_TYPE rs, UNION_TYPE rsN, _256_TYPE N_packed256, 
    _256_TYPE distm, _256_TYPE _1_distm);
void GEN_INTRINSIC(computeMXY, PRECISION)(UNION_TYPE &M_t, UNION_TYPE &X_t, UNION_TYPE &Y_t, UNION_TYPE &M_t_y, 
    UNION_TYPE M_t_2, UNION_TYPE X_t_2, UNION_TYPE Y_t_2, UNION_TYPE M_t_1, UNION_TYPE X_t_1, UNION_TYPE M_t_1_y, UNION_TYPE Y_t_1, 
    _256_TYPE pMM, _256_TYPE pGAPM, _256_TYPE pMX, _256_TYPE pXX, _256_TYPE pMY, _256_TYPE pYY, _256_TYPE distmSel);
template<class NUMBER> NUMBER GEN_INTRINSIC(compute_full_prob_avx, PRECISION) (testcase *tc, NUMBER *before_last_log = NULL);

