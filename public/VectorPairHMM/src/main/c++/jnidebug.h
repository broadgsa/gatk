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


#ifndef JNI_DEBUG_H
#define JNI_DEBUG_H

template<class NUMBER>
class DataHolder
{
#define INIT_MATRIX(X)								\
  X = new NUMBER*[m_paddedMaxReadLength];					\
  for(int i=0;i<m_paddedMaxReadLength;++i)					\
  {										\
    X[i] = new NUMBER[m_paddedMaxHaplotypeLength];					\
    for(int j=0;j<m_paddedMaxHaplotypeLength;++j)				\
      X[i][j] = (NUMBER)0; 							\
  }

#define FREE_MATRIX(X)								\
  for(int i=0;i<m_paddedMaxReadLength;++i)					\
    delete[] X[i];								\
  delete[] X;

  public:
    DataHolder() { m_is_initialized = false; }
    void initialize(int readMaxLength, int haplotypeMaxLength)
    {
      if(m_is_initialized)
      {
	FREE_MATRIX(m_matchMatrix); 
	FREE_MATRIX(m_insertionMatrix); 
	FREE_MATRIX(m_deletionMatrix); 
	FREE_MATRIX(m_prior);
	delete[] m_transition;
      }
      
      m_readMaxLength = readMaxLength;
      m_haplotypeMaxLength = haplotypeMaxLength;
      m_paddedMaxReadLength = readMaxLength + 1;
      m_paddedMaxHaplotypeLength = haplotypeMaxLength + 1;
           
      INIT_MATRIX(m_matchMatrix); 
      INIT_MATRIX(m_insertionMatrix); 
      INIT_MATRIX(m_deletionMatrix); 
      INIT_MATRIX(m_prior); 
      m_transition = new NUMBER[m_paddedMaxReadLength][6];
      for(int i=0;i<m_paddedMaxReadLength;++i)
	for(int j=0;j<6;++j)
	  m_transition[i][j] = (NUMBER)0;
      m_is_initialized = true;
    }

    //Corresponds to initializeProbabilities
    void initializeProbabilities(jint length, jbyte* insertionGOP, jbyte* deletionGOP, jbyte* overallGCP)
    {
      static unsigned g_num_prob_init = 0;
      Context<NUMBER> ctx;
      for (int r = 1; r <= length;r++)	//in original code, r < ROWS (where ROWS = paddedReadLength)
      {
	int _i = insertionGOP[r-1];	//insertionGOP
	int _d = deletionGOP[r-1];	//deletionGOP
	int _c = overallGCP[r-1];	//overallGCP
	m_transition[r][MM] = ctx._(1.0) - ctx.ph2pr[(_i + _d) & 127];	//lines 161-162
	m_transition[r][GapM] = ctx._(1.0) - ctx.ph2pr[_c];	//line 163
	m_transition[r][MX] = ctx.ph2pr[_i];	//164
	m_transition[r][XX] = ctx.ph2pr[_c];	//165
	m_transition[r][MY] =  ctx.ph2pr[_d];//last row seems different, compared to line 166
	m_transition[r][YY] = ctx.ph2pr[_c];//same as above for line 167
	//m_transition[r][MY] = (r == length) ? ctx._(1.0) : ctx.ph2pr[_d];//last row seems different, compared to line 166
	//m_transition[r][YY] = (r == length) ? ctx._(1.0) : ctx.ph2pr[_c];//same as above for line 167
#ifdef DEBUG3
	for(int j=0;j<6;++j)
	  debug_dump("transitions_jni.txt", to_string(m_transition[r][j]),true);
#endif
      }
      ++g_num_prob_init;
    }
    bool m_is_initialized;
    int m_readMaxLength;
    int m_haplotypeMaxLength;
    int m_paddedMaxReadLength;
    int m_paddedMaxHaplotypeLength;
    NUMBER** m_matchMatrix;
    NUMBER** m_insertionMatrix;
    NUMBER** m_deletionMatrix;
    NUMBER** m_prior; 
    NUMBER (*m_transition)[6];
};
extern DataHolder<double> g_double_dataholder;

template<class NUMBER>
NUMBER compute_full_prob(testcase *tc, NUMBER** M, NUMBER** X, NUMBER** Y, NUMBER (*p)[6], 
    bool do_initialization, jint hapStartIndex, NUMBER *before_last_log = NULL)
{
	int r, c;
	int ROWS = tc->rslen + 1;	//ROWS = paddedReadLength
	int COLS = tc->haplen + 1;	//COLS = paddedHaplotypeLength

	Context<NUMBER> ctx;
	//////NOTES
	////ctx.ph2pr[quality]; 	//This quantity is QualityUtils.qualToErrorProb(quality)
	////1-ctx.ph2pr[quality]; //This corresponds to QualityUtils.qualToProb(quality);

	//Initialization
	if(do_initialization)
	{
	  for (c = 0; c < COLS; c++)
	  {
	    M[0][c] = ctx._(0.0);
	    X[0][c] = ctx._(0.0);
	    Y[0][c] = ctx.INITIAL_CONSTANT / (tc->haplen);	//code from 87-90 in LoglessPairHMM
	  }

	  for (r = 1; r < ROWS; r++)
	  {
	    M[r][0] = ctx._(0.0);
	    //deletionMatrix row 0 in above nest is initialized in the Java code
	    //However, insertionMatrix column 0 is not initialized in Java code, could it be that
	    //values are re-used from a previous iteration?
	    //Why even do this, X[0][0] = 0 from above loop nest, X[idx][0] = 0 from this computation
	    X[r][0] = X[r-1][0] * p[r][XX];
	    Y[r][0] = ctx._(0.0);
	  }
	}

	for (r = 1; r < ROWS; r++)
		for (c = hapStartIndex+1; c < COLS; c++)
		{
		  	//The following lines correspond to initializePriors()
			char _rs = tc->rs[r-1];		//line 137
			char _hap = tc->hap[c-1];	//line 140
			//int _q = tc->q[r-1] & 127;	//line 138 - q is the quality (qual), should be byte hence int ANDed with 0xFF
			int _q = tc->q[r-1];	//line 138 - q is the quality (qual), should be byte hence int ANDed with 0xFF
			NUMBER distm = ctx.ph2pr[_q]; 	//This quantity is QualityUtils.qualToErrorProb(_q)
			//The assumption here is that doNotUseTristateCorrection is true	
			//TOASK
			if (_rs == _hap || _rs == 'N' || _hap == 'N')
				distm = ctx._(1.0) - distm; //This is the quantity QualityUtils.qualToProb(qual)
			else
			  distm = distm/3;
#ifdef DEBUG3
			debug_dump("priors_jni.txt",to_string(distm),true);
#endif

			//Computation inside updateCell
			M[r][c] = distm * (M[r-1][c-1] * p[r][MM] + X[r-1][c-1] * p[r][GapM] + Y[r-1][c-1] * p[r][GapM]);
			X[r][c] = M[r-1][c] * p[r][MX] + X[r-1][c] * p[r][XX];
			Y[r][c] = M[r][c-1] * p[r][MY] + Y[r][c-1] * p[r][YY];
#ifdef DEBUG3
			debug_dump("matrices_jni.txt",to_string(M[r][c]),true);
			debug_dump("matrices_jni.txt",to_string(X[r][c]),true);
			debug_dump("matrices_jni.txt",to_string(Y[r][c]),true);
#endif
		}

	NUMBER result = ctx._(0.0);
	for (c = 0; c < COLS; c++)
		result += M[ROWS-1][c] + X[ROWS-1][c];

	if (before_last_log != NULL)
		*before_last_log = result;	

#ifdef DEBUG
	debug_dump("return_values_jni.txt",to_string(ctx.LOG10(result) - ctx.LOG10_INITIAL_CONSTANT),true);
#endif
	return ctx.LOG10(result) - ctx.LOG10_INITIAL_CONSTANT;
}

#endif
