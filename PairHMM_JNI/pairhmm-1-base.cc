//#define DEBUG 1
//#define DEBUG0_1 1
//#define DEBUG3 1
#include "headers.h"
#include "template.h"
#include "utils.h"

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

  initialize_function_pointers(); 

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
}

