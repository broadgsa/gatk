#ifndef LOAD_TIME_INITIALIZER_H
#define LOAD_TIME_INITIALIZER_H
#include "headers.h"
#include <jni.h>
#include "template.h"
class LoadTimeInitializer
{
  public:
    LoadTimeInitializer();		//will be called when library is loaded
    void print_profiling();
    void debug_dump(std::string filename, std::string s, bool to_append, bool add_newline=true);
    void debug_close();
    
    void dump_sandbox(unsigned haplotypeLength, unsigned readLength, char* haplotypeBasesArray, testcase& tc);
    void open_sandbox() { m_sandbox_fptr.open("sandbox.txt", std::ios::app); }
    void close_sandbox() { m_sandbox_fptr.close(); }
    
    jfieldID m_readBasesFID;
    jfieldID m_readQualsFID;
    jfieldID m_insertionGOPFID;
    jfieldID m_deletionGOPFID;
    jfieldID m_overallGCPFID;
    jfieldID m_haplotypeBasesFID;
    //used to compute avg, variance of #testcases
    double m_sumNumReads;
    double m_sumSquareNumReads;
    double m_sumNumHaplotypes;
    double m_sumSquareNumHaplotypes;
    double m_sumNumTestcases;
    double m_sumSquareNumTestcases;
    uint64_t m_sumNumDoubleTestcases;
    double m_sumReadLengths;
    double m_sumHaplotypeLengths;
    double m_sumProductReadLengthHaplotypeLength;
    double m_sumSquareProductReadLengthHaplotypeLength;
    unsigned m_maxNumTestcases;
    unsigned m_num_invocations;
    //timing in nanoseconds
    uint64_t m_compute_time;
    uint64_t m_data_transfer_time;
    //bytes copied
    uint64_t m_bytes_copied;
  private:
    std::map<std::string, std::ofstream*> m_filename_to_fptr;
    std::ofstream m_sandbox_fptr;
};
extern LoadTimeInitializer g_load_time_initializer;


#endif
