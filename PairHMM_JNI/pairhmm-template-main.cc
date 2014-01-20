#include "headers.h"

#include "template.h"
#include "vector_defs.h"


#define BATCH_SIZE  10
#define RUN_HYBRID

int thread_level_parallelism_enabled = false ;

double getCurrClk() {
  struct timeval tv ;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

void print128b_F(__m128 x)
{	
	float *p = (float *)(&x);
	for (int i=3;i>=0;i--)
		printf("%f ", *(p+i));
	printf("\n");
}

void print128b_D(__m128d x)
{
        double *p = (double *)(&x);
        for (int i=1;i>=0;i--)
                printf("%f ", *(p+i));
        printf("\n");
}

int main()
{
/*
	IF_128f x;
	x.f = _mm_set_ps(1.0, 2.0, 3.0, 4.0);
	IF_32 shiftIn, shiftOut;
	shiftIn.f = 5.0f;
	print128b_F(x.f);
	GEN_INTRINSIC(_vector_shift, s)(x, shiftIn, shiftOut);
	print128b_F(x.f);

        IF_128d y;
        y.f = _mm_set_pd(10.0, 11.0);
        IF_64 shiftInd, shiftOutd;
        shiftInd.f = 12.0;
        print128b_D(y.f);
        GEN_INTRINSIC(_vector_shift, d)(y, shiftInd, shiftOutd);
        print128b_D(y.f);

	exit(0);
*/

        testcase tc[BATCH_SIZE];
        float result[BATCH_SIZE], result_avxf;
        double result_avxd;
        //struct timeval start, end;
	double lastClk = 0.0 ;
        double aggregateTimeRead = 0.0;
        double aggregateTimeCompute = 0.0;
        double aggregateTimeWrite = 0.0;

	// Need to call it once to initialize the static array
	ConvertChar::init() ;
	

//	char* ompEnvVar = getenv("OMP_NUM_THREADS") ;
//	if (ompEnvVar != NULL && ompEnvVar != "" && ompEnvVar != "1" ) {
//	  thread_level_parallelism_enabled = true ;
//	}

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
                //gettimeofday(&end, NULL);
                aggregateTimeRead += (getCurrClk() - lastClk) ;
		//((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

                //gettimeofday(&start, NULL);
		lastClk = getCurrClk() ;

//#pragma omp parallel for schedule(dynamic) if(thread_level_parallelism_enabled)
                for (int b=0;b<read_count;b++)
                {
                        result_avxf = GEN_INTRINSIC(GEN_INTRINSIC(compute_full_prob_, SIMD_TYPE), s)<float>(&tc[b]);

                        #ifdef RUN_HYBRID
                                #define MIN_ACCEPTED 1e-28f
                                if (result_avxf < MIN_ACCEPTED) {
				      //printf("**************** RUNNING DOUBLE ******************\n");
                                      count++;
                                      result_avxd = GEN_INTRINSIC(GEN_INTRINSIC(compute_full_prob_, SIMD_TYPE), d)<double>(&tc[b]);
                                      result[b] = log10(result_avxd) - log10(ldexp(1.0, 1020.f));
                                }
                                else
                                      result[b] = log10f(result_avxf) - log10f(ldexpf(1.f, 120.f));
                        #endif

                        #ifndef RUN_HYBRID
                                result[b] = log10f(result_avxf) - log10f(ldexpf(1.f, 120.f));
                        #endif

                }
                //gettimeofday(&end, NULL);
                aggregateTimeCompute += (getCurrClk() - lastClk) ;
		//((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

                //gettimeofday(&start, NULL);
		lastClk = getCurrClk() ;
                for (int b=0;b<read_count;b++)
                        printf("%E\n", result[b]);
                //gettimeofday(&end, NULL);
                aggregateTimeWrite += (getCurrClk() - lastClk) ;
		//((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

        }

        printf("AVX Read Time: %.2f\n", aggregateTimeRead);
        printf("AVX Compute Time: %.2f\n", aggregateTimeCompute);
        printf("AVX Write Time: %.2f\n", aggregateTimeWrite);
        printf("AVX Total Time: %.2f\n", aggregateTimeRead + aggregateTimeCompute + aggregateTimeWrite);
        printf("# Double called: %d\n", count);

        return 0;
}



