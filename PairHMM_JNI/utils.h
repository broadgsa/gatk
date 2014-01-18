#ifndef PAIRHMM_UTIL_H
#define PAIRHMM_UTIL_H

#define MIN_ACCEPTED 1e-28f
bool is_avx_supported();
bool is_sse42_supported();
extern float (*g_compute_full_prob_float)(testcase *tc, float *before_last_log);
extern double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log);
void debug_dump(std::string filename, std::string s, bool to_append, bool add_newline);
template<class NUMBER>
NUMBER compute_full_prob(testcase *tc, NUMBER *before_last_log=0);

#endif
