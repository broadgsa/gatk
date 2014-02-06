#include "headers.h"
#include "template.h"
#include "utils.h"
#include "vector_defs.h"

uint8_t ConvertChar::conversionTable[255];
float (*g_compute_full_prob_float)(testcase *tc, float* before_last_log) = 0;
double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log) = 0;

using namespace std;

bool is_avx_supported()
{
  int ecx = 0, edx = 0, ebx = 0;
  __asm__("cpuid"
      : "=b" (ebx),
      "=c" (ecx),
      "=d" (edx)
      : "a" (1)
      );
  return ((ecx >> 28)&1) == 1;
}

bool is_sse42_supported()
{
  int ecx = 0, edx = 0, ebx = 0;
  __asm__("cpuid"
      : "=b" (ebx),
      "=c" (ecx),
      "=d" (edx)
      : "a" (1)
      );
  return ((ecx >> 20)&1) == 1;
}

uint64_t get_machine_capabilities()
{
  uint64_t machine_mask = 0ull;
  if(is_avx_supported())
    machine_mask |= (1 << AVX_CUSTOM_IDX);
  if(is_sse42_supported())
    machine_mask |= (1 << SSE42_CUSTOM_IDX);
  return machine_mask;
}

void initialize_function_pointers(uint64_t mask)
{
  //mask = 0ull;
  if(is_avx_supported() && (mask & (1<< AVX_CUSTOM_IDX)))
  {
    cout << "Using AVX accelerated implementation of PairHMM\n";
    g_compute_full_prob_float = compute_full_prob_avxs<float>;
    g_compute_full_prob_double = compute_full_prob_avxd<double>;
  }
  else
    if(is_sse42_supported() && (mask & (1<< SSE42_CUSTOM_IDX)))
    {
      cout << "Using SSE4.2 accelerated implementation of PairHMM\n";
      g_compute_full_prob_float = compute_full_prob_sses<float>;
      g_compute_full_prob_double = compute_full_prob_ssed<double>;
    }
    else
    {
      cout << "Using un-vectorized C++ implementation of PairHMM\n";
      g_compute_full_prob_float = compute_full_prob<float>;
      g_compute_full_prob_double = compute_full_prob<double>;
    }
}

int normalize(char c)
{
	return ((int) (c - 33));
}

int read_testcase(testcase *tc, FILE* ifp)
{
	char *q, *i, *d, *c, *line = NULL;
	int _q, _i, _d, _c;
	int x, size = 0;
	ssize_t read;


        read = getline(&line, (size_t *) &size, ifp == 0 ? stdin : ifp);
	if (read == -1)
        {
          free(line);
          return -1;
        }


	tc->hap = (char *) malloc(size);
	tc->rs = (char *) malloc(size);
	q = (char *) malloc(size);
	i = (char *) malloc(size);
	d = (char *) malloc(size);
	c = (char *) malloc(size);

	if (sscanf(line, "%s %s %s %s %s %s\n", tc->hap, tc->rs, q, i, d, c) != 6)
		return -1;


	tc->haplen = strlen(tc->hap);
	tc->rslen = strlen(tc->rs);
        assert(strlen(q) == tc->rslen);
        assert(strlen(i) == tc->rslen);
        assert(strlen(d) == tc->rslen);
        assert(strlen(c) == tc->rslen);
	//assert(tc->rslen < MROWS);
        //tc->ihap = (int *) malloc(tc->haplen*sizeof(int));
        //tc->irs = (int *) malloc(tc->rslen*sizeof(int));

	tc->q = (char *) malloc(sizeof(char) * tc->rslen);
	tc->i = (char *) malloc(sizeof(char) * tc->rslen);
	tc->d = (char *) malloc(sizeof(char) * tc->rslen);
	tc->c = (char *) malloc(sizeof(char) * tc->rslen);

	for (x = 0; x < tc->rslen; x++)
	{
		_q = normalize(q[x]);
		_i = normalize(i[x]);
		_d = normalize(d[x]);
		_c = normalize(c[x]);
                tc->q[x] = (_q < 6) ? 6 : _q;
                //tc->q[x] = _q;
		tc->i[x] = _i;
		tc->d[x] = _d;
		tc->c[x] = _c;
                //tc->irs[x] = tc->rs[x];
	}
        //for (x = 0; x < tc->haplen; x++)
        //tc->ihap[x] = tc->hap[x];

        
	free(q);
	free(i);
	free(d);
	free(c);
	free(line);



	return 0;
}

unsigned MAX_LINE_LENGTH = 65536;
int convToInt(std::string s)
{
  int i;
  std::istringstream strin(s);
  strin >> i;
  return i;
}

void tokenize(std::ifstream& fptr, std::vector<std::string>& tokens)
{
  int i = 0;
  std::string tmp;
  std::vector<std::string> myVec;
  vector<char> line;
  line.clear();
  line.resize(MAX_LINE_LENGTH);
  vector<char> tmpline;
  tmpline.clear();
  tmpline.resize(MAX_LINE_LENGTH);
  myVec.clear();

  while(!fptr.eof())
  {
    i = 0;
    bool still_read_line = true;
    unsigned line_position = 0;
    while(still_read_line)
    {
      fptr.getline(&(tmpline[0]), MAX_LINE_LENGTH);
      if(line_position + MAX_LINE_LENGTH > line.size())
	line.resize(2*line.size());
      for(unsigned i=0;i<MAX_LINE_LENGTH && tmpline[i] != '\0';++i,++line_position)
	line[line_position] = tmpline[i];
      if(fptr.eof() || !fptr.fail()) 
      {
	still_read_line = false;
	line[line_position++] = '\0';
      }
    }
    std::istringstream kap(&(line[0]));

    while(!kap.eof())
    {
      kap >> std::skipws >> tmp;
      if(tmp != "")
      {
	myVec.push_back(tmp);
	++i;
	//std::cout <<tmp <<"#";
      }
      tmp = "";
    }
    //std::cout << "\n";
    if(myVec.size() > 0)
      break;
  }
  tokens.clear();
  //std::cout << "Why "<<myVec.size()<<"\n";
  tokens.resize(myVec.size());
  for(i=0;i<(int)myVec.size();++i)
    tokens[i] = myVec[i];
  line.clear();
  tmpline.clear();
}

int read_mod_testcase(ifstream& fptr, testcase* tc, bool reformat)
{
  static bool first_call = true;
  vector<string> tokens;
  tokens.clear();
  tokenize(fptr, tokens);
  if(tokens.size() == 0)
    return -1;
  tc->hap = new char[tokens[0].size()+2];
  tc->haplen = tokens[0].size();
  memcpy(tc->hap, tokens[0].c_str(), tokens[0].size());
  tc->rs = new char[tokens[1].size()+2];
  tc->rslen = tokens[1].size();
  tc->q = new char[tc->rslen];
  tc->i = new char[tc->rslen];
  tc->d = new char[tc->rslen];
  tc->c = new char[tc->rslen];
  //cout << "Lengths "<<tc->haplen <<" "<<tc->rslen<<"\n";
  memcpy(tc->rs, tokens[1].c_str(),tokens[1].size());
  assert(tokens.size() == 2 + 4*(tc->rslen));
  //assert(tc->rslen < MROWS);
  for(unsigned j=0;j<tc->rslen;++j)
    tc->q[j] = (char)convToInt(tokens[2+0*tc->rslen+j]);
  for(unsigned j=0;j<tc->rslen;++j)
    tc->i[j] = (char)convToInt(tokens[2+1*tc->rslen+j]);
  for(unsigned j=0;j<tc->rslen;++j)
    tc->d[j] = (char)convToInt(tokens[2+2*tc->rslen+j]);
  for(unsigned j=0;j<tc->rslen;++j)
    tc->c[j] = (char)convToInt(tokens[2+3*tc->rslen+j]);
 
  if(reformat)
  {
    ofstream ofptr;
    ofptr.open("reformat/debug_dump.txt",first_call ? ios::out : ios::app);
    assert(ofptr.is_open());
    ofptr << tokens[0] << " ";
    ofptr << tokens[1] << " ";
    for(unsigned j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->q[j]+33));
    ofptr << " ";
    for(unsigned j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->i[j]+33));
    ofptr << " ";
    for(unsigned j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->d[j]+33));
    ofptr << " ";
    for(unsigned j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->c[j]+33));
    ofptr << " 0 false\n";

    ofptr.close();
    first_call = false;
  }


  return tokens.size();
}

double getCurrClk() {
  struct timeval tv ;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

uint64_t get_time(struct timespec* store_struct)
{
  static struct timespec start_time;
  struct timespec curr_time;
  struct timespec* ptr = (store_struct == 0) ? &curr_time : store_struct;
  clock_gettime(CLOCK_REALTIME, ptr);
  uint64_t diff_time = (ptr->tv_sec-start_time.tv_sec)*1000000000+(ptr->tv_nsec-start_time.tv_nsec);
  start_time = *ptr;
  return diff_time;
}

uint64_t diff_time(struct timespec& prev_time)
{
  struct timespec curr_time;
  clock_gettime(CLOCK_REALTIME, &curr_time);
  return (uint64_t)((curr_time.tv_sec-prev_time.tv_sec)*1000000000+(curr_time.tv_nsec-prev_time.tv_nsec));
}

//#define USE_PAPI
//#define COUNT_EXCEPTIONS
//#define CHECK_RESULTS
#define CHECK_UNDERFLOW 1
#ifdef USE_PAPI
#include "papi.h"
#define NUM_PAPI_COUNTERS 4
#endif

IF_32 g_converter;
FILE* g_debug_fptr = 0;
uint64_t exceptions_array[128];
void do_compute(char* filename)
{
  //g_debug_fptr = fopen("/mnt/app_hdd/scratch/karthikg/dump.log","w");
  //assert(g_debug_fptr);
  for(unsigned i=0;i<128;++i)
    exceptions_array[i] = 0ull;
  //assert(feenableexcept(FE_DIVBYZERO | FE_INVALID) >= 0);
#ifdef USE_PAPI
  PAPI_num_counters();
  //int events[NUM_PAPI_COUNTERS] = { PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L1_DCM, PAPI_L1_ICM, PAPI_L3_TCM, PAPI_TLB_DM, PAPI_TLB_IM };
  //char* eventnames[NUM_PAPI_COUNTERS]=  { "instructions", "cycles", "l1d_misses", "l1i_misses", "l3_misses", "dtlb_misses", "itlb_misses" };
  //long long values[NUM_PAPI_COUNTERS] = { 0, 0, 0, 0, 0, 0, 0 };
  //long long accum_values[NUM_PAPI_COUNTERS] = { 0, 0, 0, 0, 0, 0, 0 };
  //int events[NUM_PAPI_COUNTERS] = { PAPI_TOT_INS, PAPI_TOT_CYC, PAPI_L1_ICM };
  //char* eventnames[NUM_PAPI_COUNTERS]=  { "instructions", "cycles", "l1i_misses"};
  //assert(PAPI_event_name_to_code("PERF_COUNT_HW_STALLED_CYCLES_FRONTEND",&(events[2])) == PAPI_OK);
  int events[NUM_PAPI_COUNTERS] = { 0, 0, 0, 0 };
    //assert(PAPI_event_name_to_code("ICACHE:IFETCH_STALL",&(events[2])) == PAPI_OK);
  //assert(PAPI_event_name_to_code("MACHINE_CLEARS:e",&(events[3])) == PAPI_OK);
  char* eventnames[NUM_PAPI_COUNTERS]=  { "instructions", "cycles", "fp_assists", "idq_ms_cycles" };
  assert(PAPI_event_name_to_code("ix86arch::INSTRUCTION_RETIRED",&(events[0])) == PAPI_OK);
  assert(PAPI_event_name_to_code("UNHALTED_REFERENCE_CYCLES",&(events[1])) == PAPI_OK);
  assert(PAPI_event_name_to_code("FP_ASSIST:ANY",   &(events[2])) == PAPI_OK);
  assert(PAPI_event_name_to_code("IDQ:MS_UOPS_CYCLES",   &(events[3])) == PAPI_OK);
  long long values[NUM_PAPI_COUNTERS] = { 0, 0, 0, 0 };
  long long accum_values[NUM_PAPI_COUNTERS] = { 0, 0, 0, 0 };

#endif
#define BATCH_SIZE 10000
  bool use_old_read_testcase = true;
  unsigned chunk_size = 100;
  std::ifstream ifptr;
  FILE* fptr = 0;
  if(use_old_read_testcase)
  {
    fptr = fopen(filename,"r");
    if(fptr == 0)
      cerr << "Could not open file "<<filename<<"\n";
    assert(fptr);
  }
  else
  {
    ifptr.open(filename);
    assert(ifptr.is_open());
  }
  vector<testcase> tc_vector;
  tc_vector.clear();
  vector<double> results_vec;
  results_vec.clear();
  vector<double> baseline_results;
  baseline_results.clear();

  bool all_ok = true;
  uint64_t total_time = 0;
  uint64_t baseline_time = 0;
  unsigned total_count = 0;
  unsigned num_testcases = 0;
  //unsigned curr_batch_size = rand()%BATCH_SIZE + 4;     //min batch size
  unsigned curr_batch_size = BATCH_SIZE;

  testcase tc_in;
  int break_value = 0;
  uint64_t fp_single_exceptions_reexecute = 0;
  uint64_t fp_single_exceptions_continue = 0;
  uint64_t num_double_executions = 0;
  while(1)
  {
    break_value = use_old_read_testcase ? read_testcase(&tc_in, fptr) : 
      read_mod_testcase(ifptr, &tc_in, true);
    tc_vector.push_back(tc_in);
    if(break_value >= 0)
      ++num_testcases;
    if(num_testcases == curr_batch_size || (break_value < 0 && num_testcases > 0))
    {
      results_vec.resize(tc_vector.size());
      baseline_results.resize(tc_vector.size());

      get_time();
#ifdef USE_PAPI
      assert(PAPI_start_counters(events, NUM_PAPI_COUNTERS) == PAPI_OK);
#endif
#pragma omp parallel for schedule(dynamic,chunk_size)  num_threads(12)
      for(unsigned i=0;i<num_testcases;++i)
      {
        double result = 0;
#ifdef COUNT_EXCEPTIONS
        fexcept_t flagp = 0;
        feclearexcept(FE_ALL_EXCEPT | __FE_DENORM);                 
#endif
        float result_avxf = g_compute_full_prob_float(&(tc_vector[i]), 0);
        //CONVERT_AND_PRINT(result_avxf);
#ifdef COUNT_EXCEPTIONS
        STORE_FP_EXCEPTIONS(flagp, exceptions_array);
        bool fp_exception =  ((flagp & (FE_UNDERFLOW|FE_OVERFLOW|FE_INVALID)) != 0);
#endif
#ifdef CHECK_UNDERFLOW
        if (result_avxf < MIN_ACCEPTED)
#else
          if(false)
#endif
          {
#ifdef COUNT_EXCEPTIONS
            if(fp_exception)
              ++fp_single_exceptions_reexecute;
#endif
            double result_avxd = g_compute_full_prob_double(&(tc_vector[i]), 0);
            result = log10(result_avxd) - log10(ldexp(1.0, 1020.0));
            ++num_double_executions;
          }
          else
          {
#ifdef COUNT_EXCEPTIONS
            if(fp_exception)
              ++fp_single_exceptions_continue;
#endif
            result = (double)(log10f(result_avxf) - log10f(ldexpf(1.f, 120.f)));
          }
        results_vec[i] = result;
      }
#ifdef USE_PAPI
      //assert(PAPI_accum_counters(values, NUM_PAPI_COUNTERS) == PAPI_OK);
      assert(PAPI_stop_counters(values, NUM_PAPI_COUNTERS) == PAPI_OK);
#endif
      total_time +=  get_time();
#ifdef USE_PAPI
      for(unsigned k=0;k<NUM_PAPI_COUNTERS;++k)
        accum_values[k] += values[k];
#endif

#ifdef CHECK_RESULTS
#pragma omp parallel for schedule(dynamic,chunk_size)
      for(unsigned i=0;i<num_testcases;++i)
      {
        testcase& tc = tc_vector[i];
        float result_avxf = compute_full_prob<float>(&tc);
        double result = 0;
        if (result_avxf < MIN_ACCEPTED) {
          double result_avxd = compute_full_prob<double>(&tc);
          result = log10(result_avxd) - log10(ldexp(1.0, 1020.0));
        }
        else
          result = (double)(log10f(result_avxf) - log10f(ldexpf(1.f, 120.f)));
        baseline_results[i] = result;
      }
      baseline_time += get_time();
      for(unsigned i=0;i<num_testcases;++i)
      {
        double baseline_result = baseline_results[i];
        double abs_error = fabs(baseline_result-results_vec[i]);
        double rel_error = (baseline_result != 0) ? fabs(abs_error/baseline_result) : 0;
        if(abs_error > 1e-5 && rel_error > 1e-5)
        {
          cout << "Line "<<total_count+i<< " " << std::scientific << baseline_result << " "<<results_vec[i]<<"\n";
          all_ok = false;
        }
      }
#else
      all_ok = false;
#endif
      for(unsigned i=0;i<num_testcases;++i)
      {
        free(tc_vector[i].rs);
        free(tc_vector[i].hap);
        free(tc_vector[i].q);
        free(tc_vector[i].i);
        free(tc_vector[i].d);
        free(tc_vector[i].c);
      }
      total_count += num_testcases;
      num_testcases = 0;
      tc_vector.clear();
      baseline_results.clear();
      results_vec.clear();
      //curr_batch_size = rand()%BATCH_SIZE + 4;     //min batch size
      curr_batch_size = BATCH_SIZE;
    }
    if(break_value < 0)
      break;
  }

  baseline_results.clear();
  results_vec.clear();
  tc_vector.clear();
  if(all_ok)
    cout << "All outputs acceptable\n";
#ifdef USE_PAPI
  for(unsigned i=0;i<NUM_PAPI_COUNTERS;++i)
    cout << eventnames[i] << " : "<<accum_values[i]<<"\n";
#endif
  cout << "Total  vector time "<< (total_time*1e-9) << " baseline time "<<baseline_time*1e-9<<"\n";
  cout.flush();
  fflush(stdout);
  if(use_old_read_testcase)
    fclose(fptr);
  else
    ifptr.close();
#ifdef COUNT_EXCEPTIONS
  cout << "Exceptions "
    <<"invalid : "<<exceptions_array[FE_INVALID]<< " "
    <<"denormal : "<<exceptions_array[__FE_DENORM]<< " "
    <<"div_by_0 : "<<exceptions_array[FE_DIVBYZERO]<< " "
    <<"overflow : "<<exceptions_array[FE_OVERFLOW]<< " "
    <<"underflow : "<<exceptions_array[FE_UNDERFLOW]<< "\n";
  cout << "Single precision FP exceptions continuations "<<fp_single_exceptions_continue<<" re-executions "<<fp_single_exceptions_reexecute<<"\n";
#endif
  cout << "Num double executions "<<num_double_executions<<"\n";

  //fclose(g_debug_fptr);
}
