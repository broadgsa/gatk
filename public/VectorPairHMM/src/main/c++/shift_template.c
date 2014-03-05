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


#ifdef PRECISION

#ifdef SIMD_ENGINE_AVX

inline void CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn, MAIN_TYPE &shiftOut)
{
    IF_128 xlow , xhigh;
    /* cast x to xlow */
    xlow.f = VEC_CAST_256_128(x.d);
    /* extract x,1 to xhigh */
    xhigh.f = VEC_EXTRACT_128(x.d, 1);
    /* extract xlow[3] */
    IF_128 shiftOutL128;
    shiftOutL128.i = _mm_srli_si128(xlow.i, SHIFT_CONST1);
    /* extract xhigh[3] */
    IF_MAIN_TYPE shiftOutH;
    shiftOutH.i = VEC_EXTRACT_UNIT(xhigh.i, SHIFT_CONST2);
    shiftOut = shiftOutH.f;
    /* shift xlow */
    xlow.i = _mm_slli_si128 (xlow.i, SHIFT_CONST3);
    /* shift xhigh */
    xhigh.i = _mm_slli_si128 (xhigh.i, SHIFT_CONST3);
    /*movss shiftIn to xlow[0] */
    _128_TYPE shiftIn128 = VEC_SET1_VAL128(shiftIn);
    xlow.f = VEC_MOVE(xlow.f , shiftIn128);
    /*movss xlow[3] to xhigh[0] */
    xhigh.f = VEC_MOVE(xhigh.f, shiftOutL128.f);
    /* cast xlow to x */
    x.d = VEC_CAST_128_256(xlow.f);
    /* insert xhigh to x,1 */
    x.d = VEC_INSERT_VAL(x.d, xhigh.f, 1);
}


inline void CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn)
{
    IF_128 xlow , xhigh;
    /* cast x to xlow */
    xlow.f = VEC_CAST_256_128(x.d);
    /* extract x,1 to xhigh */
    xhigh.f = VEC_EXTRACT_128(x.d, 1);
    /* extract xlow[3] */
    IF_128 shiftOutL128;
    shiftOutL128.i = _mm_srli_si128(xlow.i, SHIFT_CONST1);
    /* shift xlow */
    xlow.i = _mm_slli_si128 (xlow.i, SHIFT_CONST3);
    /* shift xhigh */
    xhigh.i = _mm_slli_si128 (xhigh.i, SHIFT_CONST3);
    /*movss shiftIn to xlow[0] */
    _128_TYPE shiftIn128 = VEC_SET1_VAL128(shiftIn);
    xlow.f = VEC_MOVE(xlow.f , shiftIn128);
    /*movss xlow[3] to xhigh[0] */
    xhigh.f = VEC_MOVE(xhigh.f, shiftOutL128.f);
    /* cast xlow to x */
    x.d = VEC_CAST_128_256(xlow.f);
    /* insert xhigh to x,1 */
    x.d = VEC_INSERT_VAL(x.d, xhigh.f, 1);
}

#endif

#ifdef SIMD_ENGINE_SSE

inline void CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn, MAIN_TYPE &shiftOut)
{
    IF_MAIN_TYPE tempIn, tempOut;
    tempIn.f = shiftIn;
    /* extratc H */
    tempOut.i = VEC_EXTRACT_UNIT(x.i, SHIFT_CONST1);
    shiftOut = tempOut.f;
    /* shift     */
    x.i = _mm_slli_si128(x.i, SHIFT_CONST2);
    /* insert  L */
    x.i = VEC_INSERT_UNIT(x.i , tempIn.i, SHIFT_CONST3);
}

inline void CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn)
{
    IF_MAIN_TYPE temp; temp.f = shiftIn;
    /* shift     */
    x.i = _mm_slli_si128(x.i, SHIFT_CONST2);
    /* insert  L */
    x.i = VEC_INSERT_UNIT(x.i , temp.i, SHIFT_CONST3);
}

#endif

#endif
