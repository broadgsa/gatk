#ifndef PAIRHMM_UTIL_H
#define PAIRHMM_UTIL_H

#include "template.h"

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
double getCurrClk();
uint64_t get_time(struct timespec* x=0);
uint64_t diff_time(struct timespec& prev_time);

//bit 0 is sse4.2, bit 1 is AVX
enum ProcessorCapabilitiesEnum
{
  SSE42_CUSTOM_IDX=0,
  AVX_CUSTOM_IDX
};
#define ENABLE_ALL_HARDWARE_FEATURES 0xFFFFFFFFFFFFFFFFull
uint64_t get_machine_capabilities();
void initialize_function_pointers(uint64_t mask=ENABLE_ALL_HARDWARE_FEATURES);
#endif
