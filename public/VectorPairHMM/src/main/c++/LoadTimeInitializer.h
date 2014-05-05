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


#ifndef LOAD_TIME_INITIALIZER_H
#define LOAD_TIME_INITIALIZER_H
#include "headers.h"
#include <jni.h>
/*#include "template.h"*/

enum LoadTimeInitializerStatsEnum
{
  NUM_REGIONS_IDX=0,
  NUM_READS_IDX,
  NUM_HAPLOTYPES_IDX,
  NUM_TESTCASES_IDX,
  NUM_DOUBLE_INVOCATIONS_IDX,
  HAPLOTYPE_LENGTH_IDX,
  READ_LENGTH_IDX,
  PRODUCT_READ_LENGTH_HAPLOTYPE_LENGTH_IDX,
  TOTAL_NUMBER_STATS
};
extern char* LoadTimeInitializerStatsNames[];

class LoadTimeInitializer
{
  public:
    LoadTimeInitializer();		//will be called when library is loaded
    void print_profiling();
    void debug_dump(std::string filename, std::string s, bool to_append, bool add_newline=true);
    void debug_close();
    
    void dump_sandbox(testcase& tc, unsigned tc_idx, unsigned numReads, unsigned numHaplotypes);
    void open_sandbox() { m_sandbox_fptr.open("sandbox.txt", std::ios::app); }
    void close_sandbox() { m_sandbox_fptr.close(); }
    
    jfieldID m_readBasesFID;
    jfieldID m_readQualsFID;
    jfieldID m_insertionGOPFID;
    jfieldID m_deletionGOPFID;
    jfieldID m_overallGCPFID;
    jfieldID m_haplotypeBasesFID;
    //profiling - update stats
    void update_stat(LoadTimeInitializerStatsEnum stat_idx, uint64_t value);
    //timing in nanoseconds
    uint64_t m_compute_time;
    uint64_t m_data_transfer_time;
    //bytes copied
    uint64_t m_bytes_copied;
  private:
    std::map<std::string, std::ofstream*> m_filename_to_fptr;
    std::set<std::string> m_written_files_set;
    std::ofstream m_sandbox_fptr;
    //used to compute various stats
    uint64_t m_sum_stats[TOTAL_NUMBER_STATS];
    double m_sum_square_stats[TOTAL_NUMBER_STATS];
    uint64_t m_min_stats[TOTAL_NUMBER_STATS];
    uint64_t m_max_stats[TOTAL_NUMBER_STATS];
};
extern LoadTimeInitializer g_load_time_initializer;

#define SIZE_PER_TESTCASE 6*10000
#define SIZE_PER_BUFFER 10000

#endif
