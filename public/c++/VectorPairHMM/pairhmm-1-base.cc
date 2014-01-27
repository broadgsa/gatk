//#define DEBUG 1
//#define DEBUG0_1 1
//#define DEBUG3 1
#include "headers.h"
#include "template.h"
#include "utils.h"
#include "LoadTimeInitializer.h"

using namespace std;


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
  unsigned chunk_size = 100;
  if(argc >= 4)
    chunk_size = strtol(argv[3],0,10);

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

  vector<testcase> tc_vector;
  tc_vector.clear();
  testcase tc;
  uint64_t total_time = 0;
  while(1)
  {
    int break_value = use_old_read_testcase ? read_testcase(&tc, fptr) : read_mod_testcase(ifptr,&tc,true);
    if(break_value >= 0)
      tc_vector.push_back(tc);
    if(tc_vector.size() == BATCH_SIZE || (break_value < 0 && tc_vector.size() > 0))
    {
      vector<double> results_vec;
      results_vec.clear();
      results_vec.resize(tc_vector.size());
      get_time();
#pragma omp parallel for schedule(dynamic,chunk_size)  num_threads(12)
      for(unsigned i=0;i<tc_vector.size();++i)
      {
	testcase& tc = tc_vector[i];
	float result_avxf = g_compute_full_prob_float(&tc, 0);
	double result = 0;
	if (result_avxf < MIN_ACCEPTED) {
	  double result_avxd = g_compute_full_prob_double(&tc, 0);
	  result = log10(result_avxd) - log10(ldexp(1.0, 1020.0));
	}
	else
	  result = (double)(log10f(result_avxf) - log10f(ldexpf(1.f, 120.f)));

	results_vec[i] = result;
      }
      total_time +=  get_time();
#pragma omp parallel for schedule(dynamic,chunk_size)
      for(unsigned i=0;i<tc_vector.size();++i)
      {
	testcase& tc = tc_vector[i];
	double baseline_result = compute_full_prob<double>(&tc);
	baseline_result = log10(baseline_result) - log10(ldexp(1.0, 1020.0)); 
	double abs_error = fabs(baseline_result-results_vec[i]);
	double rel_error = (baseline_result != 0) ? fabs(abs_error/baseline_result) : 0;
	if(abs_error > 1e-5 && rel_error > 1e-5)
	  cout << std::scientific << baseline_result << " "<<results_vec[i]<<"\n";
	delete tc_vector[i].rs;
	delete tc_vector[i].hap;
	delete tc_vector[i].q;
	delete tc_vector[i].i;
	delete tc_vector[i].d;
	delete tc_vector[i].c;
      }
      results_vec.clear();
      tc_vector.clear();
    }
    if(break_value < 0)
      break;
  }
  cout << "Total time "<< ((double)total_time)/1e9 << "\n";
  if(use_old_read_testcase)
    fclose(fptr);
  else
    ifptr.close();
  return 0;  
}

