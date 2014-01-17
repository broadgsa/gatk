//#define DEBUG 1
//#define DEBUG0_1 1
//#define DEBUG3 1
#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <emmintrin.h>
#include "template.h"

//#include "define-float.h"
//#include "shift_template.c"
//#include "pairhmm-template-kernel.cc"

#include "define-double.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"


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


template<class NUMBER>
NUMBER compute_full_prob(testcase *tc, NUMBER *before_last_log = NULL)
{
	int r, c;
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	Context<NUMBER> ctx;

	NUMBER M[MROWS][MCOLS];
	NUMBER X[MROWS][MCOLS];
	NUMBER Y[MROWS][MCOLS];
	NUMBER p[MROWS][6];

	p[0][MM] = ctx._(0.0);
	p[0][GapM] = ctx._(0.0);
	p[0][MX] = ctx._(0.0);
	p[0][XX] = ctx._(0.0);
	p[0][MY] = ctx._(0.0);
	p[0][YY] = ctx._(0.0);
	for (r = 1; r < ROWS; r++)
	{
		int _i = tc->i[r-1] & 127;
		int _d = tc->d[r-1] & 127;
		int _c = tc->c[r-1] & 127;
		p[r][MM] = ctx._(1.0) - ctx.ph2pr[(_i + _d) & 127];
		p[r][GapM] = ctx._(1.0) - ctx.ph2pr[_c];
		p[r][MX] = ctx.ph2pr[_i];
		p[r][XX] = ctx.ph2pr[_c];
		p[r][MY] = ctx.ph2pr[_d];
		p[r][YY] = ctx.ph2pr[_c];
		//p[r][MY] = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_d];
		//p[r][YY] = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_c];
	}

	for (c = 0; c < COLS; c++)
	{
		M[0][c] = ctx._(0.0);
		X[0][c] = ctx._(0.0);
		Y[0][c] = ctx.INITIAL_CONSTANT / (tc->haplen);
	}

	for (r = 1; r < ROWS; r++)
	{
		M[r][0] = ctx._(0.0);
		X[r][0] = X[r-1][0] * p[r][XX];
		Y[r][0] = ctx._(0.0);
	}

	 NUMBER result = ctx._(0.0);

	for (r = 1; r < ROWS; r++)
		for (c = 1; c < COLS; c++)
		{
			char _rs = tc->rs[r-1];
			char _hap = tc->hap[c-1];
			int _q = tc->q[r-1] & 127;
			NUMBER distm = ctx.ph2pr[_q];
			if (_rs == _hap || _rs == 'N' || _hap == 'N')
				distm = ctx._(1.0) - distm;
			else
			  distm = distm/3;
			M[r][c] = distm * (M[r-1][c-1] * p[r][MM] + X[r-1][c-1] * p[r][GapM] + Y[r-1][c-1] * p[r][GapM]);
			X[r][c] = M[r-1][c] * p[r][MX] + X[r-1][c] * p[r][XX];
			Y[r][c] = M[r][c-1] * p[r][MY] + Y[r][c-1] * p[r][YY];
		}

	for (c = 0; c < COLS; c++)
	{
		result += M[ROWS-1][c] + X[ROWS-1][c];
	}

	if (before_last_log != NULL)
		*before_last_log = result;

	return ctx.LOG10(result) - ctx.LOG10_INITIAL_CONSTANT;
}

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

  testcase tc;
  if(use_old_read_testcase)
  {
    FILE* fptr = fopen(argv[1],"r");
    while(!feof(fptr))
    {
      if(read_testcase(&tc, fptr) >= 0)
      {
	double result_avxd = GEN_INTRINSIC(compute_full_prob_avx, d)<double>(&tc);
	double result = log10(result_avxd) - log10(ldexp(1.0, 1020));

	cout << std::scientific << compute_full_prob<double>(&tc) << " "<<result<<"\n";
	delete tc.rs;
	delete tc.hap;
      }
    }
    fclose(fptr);
  }
  else
  {
    std::ifstream ifptr;
    std::vector<std::string> tokens;
    ifptr.open(argv[1]);
    assert(ifptr.is_open());
    while(1)
    {
      tokens.clear();
      if(read_mod_testcase(ifptr, &tc, false) < 0)
	break;
      //double result = 0;
      double result_avxd = GEN_INTRINSIC(compute_full_prob_avx, d)<double>(&tc);
      double result = log10(result_avxd) - log10(ldexp(1.0, 1020));

      cout << std::scientific << compute_full_prob<double>(&tc) << " "<<result<<"\n";
      delete tc.rs;
      delete tc.hap;
    }
    ifptr.close();
  }
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

