//#define DEBUG 1
//#define DEBUG0_1 1
//#define DEBUG3 1
#include "headers.h"
#include "template.h"
#include "utils.h"

#include "define-float.h"
#include "avx_function_prototypes.h"

#include "define-double.h"
#include "avx_function_prototypes.h"

using namespace std;
class LoadTimeInitializer
{
  public:
    LoadTimeInitializer()		//will be called when library is loaded
    {
      ConvertChar::init();
    }
};
LoadTimeInitializer g_load_time_initializer;

#define BATCH_SIZE  10000
#define RUN_HYBRID

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    cerr << "Needs path to input file as argument\n";
    exit(0);
  }
  bool use_old_read_testcase = false;
  if(argc >= 3 && string(argv[2]) == "1")
    use_old_read_testcase = true;

  if(true)
  {
    g_compute_full_prob_double = GEN_INTRINSIC(compute_full_prob_avx, d)<double>;
    g_compute_full_prob_float = GEN_INTRINSIC(compute_full_prob_avx, s)<float>;
  }
  else
  {
    g_compute_full_prob_double = compute_full_prob<double>;
    g_compute_full_prob_float = compute_full_prob<float>;
  }

  std::ifstream ifptr;
  FILE* fptr = 0;
  if(use_old_read_testcase)
  {
    fptr = fopen(argv[1],"r");
    assert(fptr);
  }
  else
  {
    ifptr.open(argv[1]);
    assert(ifptr.is_open());
  }

  testcase tc;
  while(1)
  {
    int break_value = use_old_read_testcase ? read_testcase(&tc, fptr) : read_mod_testcase(ifptr,&tc,true);
    if(break_value < 0)
      break;
    float result_avxf = g_compute_full_prob_float(&tc, 0);
    double result = 0;
    if (result_avxf < MIN_ACCEPTED) {
      double result_avxd = g_compute_full_prob_double(&tc, 0);
      result = log10(result_avxd) - log10(ldexp(1.0, 1020.0));
    }
    else
      result = (double)(log10f(result_avxf) - log10f(ldexpf(1.f, 120.f)));

    double baseline_result = compute_full_prob<double>(&tc);
    baseline_result = log10(baseline_result) - log10(ldexp(1.0, 1020.0)); 
    cout << std::scientific << baseline_result << " "<<result<<"\n";
    delete tc.rs;
    delete tc.hap;
  }
  if(use_old_read_testcase)
    fclose(fptr);
  else
    ifptr.close();
  return 0;  

#if 0
  float result[BATCH_SIZE], result_avxf;
  double result_avxd;
  struct timeval start, end;
  long long aggregateTimeRead = 0L;
  long long aggregateTimeCompute = 0L;
  long long aggregateTimeWrite = 0L;

  bool noMoreData = false;
  int count =0;
  while (!noMoreData)
  {
    int read_count = BATCH_SIZE;
    gettimeofday(&start, NULL);
    for (int b=0;b<BATCH_SIZE;b++)
      if (read_testcase(&tc[b], NULL)==-1)
      {
	read_count = b;
	noMoreData = true;
	break;
      }
    gettimeofday(&end, NULL);
    aggregateTimeRead += ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

    gettimeofday(&start, NULL);
    for (int b=0;b<read_count;b++)
    {
      result_avxf = compute_full_prob<float>(&tc[b]);

#ifdef RUN_HYBRID
#define MIN_ACCEPTED 1e-28f
      if (result_avxf < MIN_ACCEPTED) {
	count++;
	result_avxd = compute_full_prob<double>(&tc[b]);
	result[b] = log10(result_avxd) - log10(ldexp(1.0, 1020.f));
      }
      else
	result[b] = log10f(result_avxf) - log10f(ldexpf(1.f, 120.f));
#endif

#ifndef RUN_HYBRID
      result[b] = log10f(result_avxf) - log10f(ldexpf(1.f, 120.f));
#endif

    }
    gettimeofday(&end, NULL);
    aggregateTimeCompute += ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

    gettimeofday(&start, NULL);
    for (int b=0;b<read_count;b++)
      printf("%E\n", result[b]);
    gettimeofday(&end, NULL);
    aggregateTimeWrite += ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

  }

  printf("AVX Read Time: %ld\n", aggregateTimeRead);
  printf("AVX Compute Time: %ld\n", aggregateTimeCompute);
  printf("AVX Write Time: %ld\n", aggregateTimeWrite);
  printf("AVX Total Time: %ld\n", aggregateTimeRead + aggregateTimeCompute + aggregateTimeWrite);
  printf("# Double called: %d\n", count);

  return 0;
#endif
}

