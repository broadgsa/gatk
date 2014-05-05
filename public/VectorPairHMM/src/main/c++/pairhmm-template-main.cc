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


#include "headers.h"

#define SIMD_ENGINE avx
#define SIMD_ENGINE_AVX

#include "utils.h"

#define BATCH_SIZE  10000
#define RUN_HYBRID

double getCurrClk();
int thread_level_parallelism_enabled = false ;


int main()
{
    testcase* tc = new testcase[BATCH_SIZE];
    float result[BATCH_SIZE], result_avxf;
    double result_avxd;
    double lastClk = 0.0 ;
    double aggregateTimeRead = 0.0;
    double aggregateTimeCompute = 0.0;
    double aggregateTimeWrite = 0.0;

    // Need to call it once to initialize the static array
    ConvertChar::init() ;

    //      char* ompEnvVar = getenv("OMP_NUM_THREADS") ;
    //      if (ompEnvVar != NULL && ompEnvVar != "" && ompEnvVar != "1" ) {
    //        thread_level_parallelism_enabled = true ;
    //      }

    bool noMoreData = false;
    int count =0;
    while (!noMoreData)
    {
        int read_count = BATCH_SIZE;

        lastClk = getCurrClk() ;
        for (int b=0;b<BATCH_SIZE;b++)
            if (read_testcase(&tc[b])==-1)
            {
                read_count = b;
                noMoreData = true;
                break;
            }
        aggregateTimeRead += (getCurrClk() - lastClk) ;
        lastClk = getCurrClk() ;

        //#pragma omp parallel for schedule(dynamic) if(thread_level_parallelism_enabled)
        for (int b=0;b<read_count;b++)
        {
            result_avxf = CONCAT(CONCAT(compute_full_prob_,SIMD_ENGINE), s)<float>(&tc[b]);

#ifdef RUN_HYBRID
#define MIN_ACCEPTED 1e-28f
            if (result_avxf < MIN_ACCEPTED) {
                count++;
                result_avxd = CONCAT(CONCAT(compute_full_prob_,SIMD_ENGINE), d)<double>(&tc[b]);
                result[b] = log10(result_avxd) - log10(ldexp(1.0, 1020.f));
            }
            else
                result[b] = log10f(result_avxf) - log10f(ldexpf(1.f, 120.f));
#endif

#ifndef RUN_HYBRID
            result[b] = log10f(result_avxf) - log10f(ldexpf(1.f, 120.f));
#endif
        }
        aggregateTimeCompute += (getCurrClk() - lastClk) ;
        lastClk = getCurrClk() ;
        for (int b=0;b<read_count;b++)
            printf("%E\n", result[b]);
        aggregateTimeWrite += (getCurrClk() - lastClk) ;
    }

    delete[] tc;
    printf("AVX Read Time: %.2f\n", aggregateTimeRead);
    printf("AVX Compute Time: %.2f\n", aggregateTimeCompute);
    printf("AVX Write Time: %.2f\n", aggregateTimeWrite);
    printf("AVX Total Time: %.2f\n", aggregateTimeRead + aggregateTimeCompute + aggregateTimeWrite);
    printf("# Double called: %d\n", count);

    return 0;
}



