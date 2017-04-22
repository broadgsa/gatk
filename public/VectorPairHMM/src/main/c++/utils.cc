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
#include "utils.h"
#include "LoadTimeInitializer.h"
using namespace std;

//static members from ConvertChar
uint8_t ConvertChar::conversionTable[255];
//Global function pointers in utils.h
float (*g_compute_full_prob_float)(testcase *tc, float* before_last_log) = 0;
double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log) = 0;
//Static members in ContextBase
template<>
bool ContextBase<double>::staticMembersInitializedFlag = false;
template<>
double ContextBase<double>::jacobianLogTable[JACOBIAN_LOG_TABLE_SIZE] = { };
template<>
double ContextBase<double>::matchToMatchProb[((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] = { };
template<>
bool ContextBase<float>::staticMembersInitializedFlag = false;
template<>
float ContextBase<float>::jacobianLogTable[JACOBIAN_LOG_TABLE_SIZE] = { };
template<>
float ContextBase<float>::matchToMatchProb[((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] = { };


bool search_file_for_string(string filename, string search_string)
{
  ifstream fptr;
  fptr.open(filename.c_str(),ios::in);
  if(fptr.is_open())
  {
    string buffer;
    buffer.clear();
    buffer.resize(4096);
    bool retvalue = false;
    while(!fptr.eof())
    {
      fptr.getline(&(buffer[0]), 4096);
      if(buffer.find(search_string) != string::npos)    //found string
      {
        retvalue = true;
        break;
      }
    }
    buffer.clear();
    fptr.close();
    return retvalue;
  }
  else
    return false;
}

bool is_cpuid_ecx_bit_set(int eax, int bitidx)
{
  int ecx = 0, edx = 0, ebx = 0;
  __asm__ ("cpuid"
      :"=b" (ebx),
      "=c" (ecx),
      "=d" (edx)
      :"a" (eax)
      );
  return (((ecx >> bitidx)&1) == 1);
}

bool is_avx_supported()
{
#ifdef __INTEL_COMPILER
  bool use_avx = _may_i_use_cpu_feature(_FEATURE_AVX);
  if(use_avx)
    return true;
  else
  {
    //check if core supports AVX, but kernel does not and print info message
    if(!is_cpuid_ecx_bit_set(1, 28))  //core does not support AVX
      return false;
    //else fall through to end of function
  }
#else
  if(!__builtin_cpu_supports("avx"))  //core does not support AVX
    return false;
  else
  {
    //core supports AVX, check if kernel supports
    if(search_file_for_string("/proc/cpuinfo","avx"))
      return true;
    //else fall through to end of function

  }
#endif  //__INTEL_COMPILER
  clog << "INFO: Your CPU supports AVX vector instructions, but your kernel does not. Try upgrading to a kernel that supports AVX.\n";
  clog << "INFO: Your program will run correctly, but slower than the AVX version\n";
  return false;
}

bool is_sse41_supported()
{
#ifdef __INTEL_COMPILER
  return  (_may_i_use_cpu_feature(_FEATURE_SSE4_1) > 0);
#else
  return  __builtin_cpu_supports("sse4.1");
#endif
  //return is_cpuid_ecx_bit_set(1, 19);
}

bool is_sse42_supported()
{
#ifdef __INTEL_COMPILER
  return  (_may_i_use_cpu_feature(_FEATURE_SSE4_2) > 0);
#else
  return  __builtin_cpu_supports("sse4.2");
#endif
  //return is_cpuid_ecx_bit_set(1, 20);
}

uint64_t get_machine_capabilities()
{
  uint64_t machine_mask = 0ull;
  if(is_avx_supported())
    machine_mask |= (1 << AVX_CUSTOM_IDX);
  if(is_sse42_supported())
    machine_mask |= (1 << SSE42_CUSTOM_IDX);
  if(is_sse41_supported())
    machine_mask |= (1 << SSE41_CUSTOM_IDX);
  return machine_mask;
}

void initialize_function_pointers(uint64_t mask)
{
  //mask = 0ull;
  //mask = (1 << SSE41_CUSTOM_IDX);
  if(is_avx_supported() && (mask & (1<< AVX_CUSTOM_IDX)))
  {
    cerr << "Using AVX accelerated implementation of PairHMM\n";
    g_compute_full_prob_float = compute_full_prob_avxs<float>;
    g_compute_full_prob_double = compute_full_prob_avxd<double>;
  }
  else
    if(is_sse41_supported() && (mask & ((1<< SSE41_CUSTOM_IDX) | (1<<SSE42_CUSTOM_IDX))))
    {
      cerr << "Using SSE4.1 accelerated implementation of PairHMM\n";
      g_compute_full_prob_float = compute_full_prob_sses<float>;
      g_compute_full_prob_double = compute_full_prob_ssed<double>;
    }
    else
    {
      cerr << "Using un-vectorized C++ implementation of PairHMM\n";
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
        assert(strlen(q) == (size_t)tc->rslen);
        assert(strlen(i) == (size_t)tc->rslen);
        assert(strlen(d) == (size_t)tc->rslen);
        assert(strlen(c) == (size_t)tc->rslen);

        g_load_time_initializer.update_stat(READ_LENGTH_IDX, tc->rslen); 
        g_load_time_initializer.update_stat(HAPLOTYPE_LENGTH_IDX, tc->haplen);
        g_load_time_initializer.update_stat(PRODUCT_READ_LENGTH_HAPLOTYPE_LENGTH_IDX, tc->haplen*tc->rslen);
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
	//std::cerr <<tmp <<"#";
      }
      tmp = "";
    }
    //std::cerr << "\n";
    if(myVec.size() > 0)
      break;
  }
  tokens.clear();
  //std::cerr << "Why "<<myVec.size()<<"\n";
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
  //cerr << "Lengths "<<tc->haplen <<" "<<tc->rslen<<"\n";
  memcpy(tc->rs, tokens[1].c_str(),tokens[1].size());
  assert(tokens.size() == (size_t)(2 + 4*(tc->rslen)));
  //assert(tc->rslen < MROWS);
  for(int j=0;j<tc->rslen;++j)
    tc->q[j] = (char)convToInt(tokens[2+0*tc->rslen+j]);
  for(int j=0;j<tc->rslen;++j)
    tc->i[j] = (char)convToInt(tokens[2+1*tc->rslen+j]);
  for(int j=0;j<tc->rslen;++j)
    tc->d[j] = (char)convToInt(tokens[2+2*tc->rslen+j]);
  for(int j=0;j<tc->rslen;++j)
    tc->c[j] = (char)convToInt(tokens[2+3*tc->rslen+j]);
 
  if(reformat)
  {
    ofstream ofptr;
    ofptr.open("reformat/debug_dump.txt",first_call ? ios::out : ios::app);
    assert(ofptr.is_open());
    ofptr << tokens[0] << " ";
    ofptr << tokens[1] << " ";
    for(int j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->q[j]+33));
    ofptr << " ";
    for(int j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->i[j]+33));
    ofptr << " ";
    for(int j=0;j<tc->rslen;++j)
      ofptr << ((char)(tc->d[j]+33));
    ofptr << " ";
    for(int j=0;j<tc->rslen;++j)
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

inline unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

void get_time(struct timespec* store_struct)
{
  clock_gettime(CLOCK_REALTIME, store_struct);
}

uint64_t diff_time(struct timespec& prev_time)
{
  struct timespec curr_time;
  clock_gettime(CLOCK_REALTIME, &curr_time);
  return (uint64_t)((curr_time.tv_sec-prev_time.tv_sec)*1000000000+(curr_time.tv_nsec-prev_time.tv_nsec));
}


#ifdef USE_PAPI
#include "papi.h"
#define NUM_PAPI_COUNTERS 4
#endif

void do_compute(char* filename, bool use_old_read_testcase, unsigned chunk_size, bool do_check)
{
  FILE* fptr = 0;
  ifstream ifptr;
  if(use_old_read_testcase)
  {
    fptr = fopen(filename,"r");
    assert(fptr);
  }
  else
  {
    ifptr.open(filename);
    assert(ifptr.is_open());
  }
#ifdef PRINT_PER_INTERVAL_TIMINGS
  ofstream times_fptr;
  times_fptr.open("native_timed_intervals.csv",ios::out);
#endif
  vector<testcase> tc_vector;
  tc_vector.clear();
  testcase tc;
  uint64_t vector_compute_time = 0;
  uint64_t baseline_compute_time = 0;
  uint64_t num_double_calls = 0;
  unsigned num_testcases = 0;
  bool all_ok = do_check ? true : false;
#ifdef USE_PAPI
  uint32_t all_mask = (0);
  uint32_t no_usr_mask = (1 << 16);    //bit 16 user mode, bit 17 kernel mode
  uint32_t no_kernel_mask = (1 << 17);    //bit 16 user mode, bit 17 kernel mode
  PAPI_num_counters();
  int events[NUM_PAPI_COUNTERS] = { 0, 0, 0, 0 };
  char* eventnames[NUM_PAPI_COUNTERS]=  { "cycles", "l1_pending_miss", "lfb_hit", "l2_hit" };
  assert(PAPI_event_name_to_code("UNHALTED_REFERENCE_CYCLES:u=1:k=1",&(events[0])) == PAPI_OK);
  assert(PAPI_event_name_to_code("L1D_PEND_MISS:OCCURRENCES",   &(events[1])) == PAPI_OK);
  assert(PAPI_event_name_to_code("MEM_LOAD_UOPS_RETIRED:HIT_LFB",   &(events[2])) == PAPI_OK);
  assert(PAPI_event_name_to_code("MEM_LOAD_UOPS_RETIRED:L2_HIT",   &(events[3])) == PAPI_OK);
  long long values[NUM_PAPI_COUNTERS] = { 0, 0, 0, 0 };
  long long accum_values[NUM_PAPI_COUNTERS] = { 0, 0, 0, 0 };
#endif
  while(1)
  {
    int break_value = use_old_read_testcase ? read_testcase(&tc, fptr) : read_mod_testcase(ifptr,&tc,true);
    if(break_value >= 0)
      tc_vector.push_back(tc);
    if(tc_vector.size() == BATCH_SIZE || (break_value < 0 && tc_vector.size() > 0))
    {
      vector<double> results_vec;
      vector<double> baseline_results_vec;
      results_vec.clear();
      baseline_results_vec.clear();
      results_vec.resize(tc_vector.size());
      baseline_results_vec.resize(tc_vector.size());
      g_load_time_initializer.update_stat(NUM_TESTCASES_IDX, tc_vector.size()); 
      g_load_time_initializer.update_stat(NUM_READS_IDX, tc_vector.size()); 
      g_load_time_initializer.update_stat(NUM_HAPLOTYPES_IDX, tc_vector.size()); 
      struct timespec start_time;
#ifdef USE_PAPI
      assert(PAPI_start_counters(events, NUM_PAPI_COUNTERS) == PAPI_OK);
#endif
      get_time(&start_time);
#pragma omp parallel for schedule(dynamic,chunk_size)  num_threads(12)
#ifdef DO_REPEAT_PROFILING
      for(unsigned z=0;z<10;++z)
#endif
      {
        for(unsigned i=0;i<tc_vector.size();++i)
        {
          testcase& tc = tc_vector[i];
          float result_avxf = g_compute_full_prob_float(&tc, 0);
          double result = 0;
          if (result_avxf < MIN_ACCEPTED) {
            double result_avxd = g_compute_full_prob_double(&tc, 0);
            result = log10(result_avxd) - log10(ldexp(1.0, 1020.0));
            ++num_double_calls;
          }
          else
            result = (double)(log10f(result_avxf) - log10f(ldexpf(1.f, 120.f)));
#ifdef DUMP_COMPUTE_VALUES
          g_load_time_initializer.debug_dump("return_values_vector.txt",to_string(result),true);
#endif
          results_vec[i] = result;
        }
      }
#ifdef USE_PAPI
      assert(PAPI_stop_counters(values, NUM_PAPI_COUNTERS) == PAPI_OK);
#endif
      uint64_t curr_interval = diff_time(start_time);
#ifdef PRINT_PER_INTERVAL_TIMINGS
      times_fptr << curr_interval << "\n";
#endif
      vector_compute_time +=  curr_interval;
#ifdef USE_PAPI
      for(unsigned k=0;k<NUM_PAPI_COUNTERS;++k)
        accum_values[k] += values[k];
#endif
      num_testcases += tc_vector.size();
      if(do_check)
      {
        get_time(&start_time);
#pragma omp parallel for schedule(dynamic,chunk_size)
        for(unsigned i=0;i<tc_vector.size();++i)
        {
          testcase& tc = tc_vector[i];
          double baseline_result = compute_full_prob<double>(&tc);
          baseline_result = log10(baseline_result) - log10(ldexp(1.0, 1020.0));
          baseline_results_vec[i] = baseline_result;
        }
        baseline_compute_time += diff_time(start_time);
        for(unsigned i=0;i<tc_vector.size();++i)
        {
          double baseline_result = baseline_results_vec[i];
          double abs_error = fabs(baseline_result-results_vec[i]);
          double rel_error = (baseline_result != 0) ? fabs(abs_error/baseline_result) : 0;
          if(abs_error > 1e-5 && rel_error > 1e-5)
          {
            cerr << std::scientific << baseline_result << " "<<results_vec[i]<<"\n";
            all_ok = false;
          }
        }
      }
      for(unsigned i=0;i<tc_vector.size();++i)
      {
	delete[] tc_vector[i].rs;
	delete[] tc_vector[i].hap;
	delete[] tc_vector[i].q;
	delete[] tc_vector[i].i;
	delete[] tc_vector[i].d;
	delete[] tc_vector[i].c;
      }
      results_vec.clear();
      tc_vector.clear();
    }
    if(break_value < 0)
      break;
  }
#ifdef DUMP_COMPUTE_VALUES
  g_load_time_initializer.debug_close();
#endif
  if(all_ok)
  {
    cerr << "All output values within acceptable error\n";
    cerr << "Baseline double precision compute time "<<baseline_compute_time*1e-9<<"\n";
  }
  cerr << "Num testcase "<<num_testcases<< " num double invocations "<<num_double_calls<<"\n";
  cerr << "Vector compute time "<< vector_compute_time*1e-9 << "\n";
#ifdef USE_PAPI
  for(unsigned i=0;i<NUM_PAPI_COUNTERS;++i)
    cerr << eventnames[i] << " : "<<accum_values[i]<<"\n";
#endif
#ifdef PRINT_PER_INTERVAL_TIMINGS
  times_fptr.close();
#endif
  if(use_old_read_testcase)
    fclose(fptr);
  else
    ifptr.close();
  g_load_time_initializer.print_profiling();
}
