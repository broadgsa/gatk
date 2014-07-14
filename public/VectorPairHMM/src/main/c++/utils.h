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


#ifndef PAIRHMM_UTIL_H
#define PAIRHMM_UTIL_H

#include "common_data_structure.h"

template<class T>
std::string to_string(T obj)
{
  std::stringstream ss;
  std::string ret_string;
  ss.clear();
  ss << std::scientific << obj;
  ss >> ret_string;
  ss.clear();
  return ret_string;
}
void debug_dump(std::string filename, std::string s, bool to_append, bool add_newline=true);

int read_mod_testcase(std::ifstream& fptr, testcase* tc, bool reformat=false);

bool is_avx_supported();
bool is_sse42_supported();
extern float (*g_compute_full_prob_float)(testcase *tc, float *before_last_log);
extern double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log);
void debug_dump(std::string filename, std::string s, bool to_append, bool add_newline);
template<class NUMBER>
NUMBER compute_full_prob(testcase *tc, NUMBER *before_last_log=0);
template<class NUMBER>
NUMBER compute_full_prob_avxd(testcase *tc, NUMBER *before_last_log=0);
template<class NUMBER>
NUMBER compute_full_prob_avxs(testcase *tc, NUMBER *before_last_log=0);
template<class NUMBER>
NUMBER compute_full_prob_ssed(testcase *tc, NUMBER *before_last_log=0);
template<class NUMBER>
NUMBER compute_full_prob_sses(testcase *tc, NUMBER *before_last_log=0);

double getCurrClk();
void get_time(struct timespec* x);
uint64_t diff_time(struct timespec& prev_time);

//bit 0 is sse4.2, bit 1 is AVX
enum ProcessorCapabilitiesEnum
{
  SSE41_CUSTOM_IDX=0,
  SSE42_CUSTOM_IDX,
  AVX_CUSTOM_IDX
};
#define ENABLE_ALL_HARDWARE_FEATURES 0xFFFFFFFFFFFFFFFFull
uint64_t get_machine_capabilities();
void initialize_function_pointers(uint64_t mask=ENABLE_ALL_HARDWARE_FEATURES);
void do_compute(char* filename, bool use_old_read_testcase=true, unsigned chunk_size=10000, bool do_check=true);

//#define DO_WARMUP
//#define DO_REPEAT_PROFILING
/*#define DUMP_COMPUTE_VALUES 1*/
#define BATCH_SIZE  10000
#define RUN_HYBRID
/*#define PRINT_PER_INTERVAL_TIMINGS 1*/

#endif
