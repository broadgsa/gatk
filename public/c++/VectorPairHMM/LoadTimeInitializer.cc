#include "LoadTimeInitializer.h"
#include "utils.h"
using namespace std;

LoadTimeInitializer g_load_time_initializer;

LoadTimeInitializer::LoadTimeInitializer()		//will be called when library is loaded
{
  ConvertChar::init();
  m_sumNumReads = 0;
  m_sumSquareNumReads = 0;
  m_sumNumHaplotypes = 0;
  m_sumSquareNumHaplotypes = 0;
  m_sumNumTestcases = 0;
  m_sumNumDoubleTestcases = 0;
  m_sumSquareNumTestcases = 0;
  m_sumReadLengths = 0;
  m_sumHaplotypeLengths = 0;
  m_sumProductReadLengthHaplotypeLength = 0;
  m_sumSquareProductReadLengthHaplotypeLength = 0;
  m_maxNumTestcases = 0;
  m_num_invocations = 0;

  m_compute_time = 0;
  m_data_transfer_time = 0;
  m_bytes_copied = 0;

  m_filename_to_fptr.clear();

  initialize_function_pointers();
  cout.flush();
}

void LoadTimeInitializer::print_profiling()
{
  double mean_val;
  cout << "Compute time "<<m_compute_time*1e-9<<"\n";
  cout << "Data initialization time "<<m_data_transfer_time*1e-9<<"\n";
  cout <<"Invocations : "<<m_num_invocations<<"\n";
  cout << "term\tsum\tsumSq\tmean\tvar\tmax\n";
  mean_val = m_sumNumReads/m_num_invocations;
  cout << "reads\t"<<m_sumNumReads<<"\t"<<m_sumSquareNumReads<<"\t"<<mean_val<<"\t"<<
    (m_sumSquareNumReads/m_num_invocations)-mean_val*mean_val<<"\n";
  mean_val = m_sumNumHaplotypes/m_num_invocations;
  cout << "haplotypes\t"<<m_sumNumHaplotypes<<"\t"<<m_sumSquareNumHaplotypes<<"\t"<<mean_val<<"\t"<<
    (m_sumSquareNumHaplotypes/m_num_invocations)-mean_val*mean_val<<"\n";
  mean_val = m_sumNumTestcases/m_num_invocations;
  cout << "numtestcases\t"<<m_sumNumTestcases<<"\t"<<m_sumSquareNumTestcases<<"\t"<<mean_val<<"\t"<<
    (m_sumSquareNumTestcases/m_num_invocations)-mean_val*mean_val<<"\t"<<m_maxNumTestcases<<"\n";
  mean_val = m_sumProductReadLengthHaplotypeLength/m_sumNumTestcases;
  cout <<"productReadLengthHaplotypeLength\t"<<m_sumProductReadLengthHaplotypeLength<<"\t"<<m_sumSquareProductReadLengthHaplotypeLength<<"\t"
    <<mean_val<<"\t"<<(m_sumSquareProductReadLengthHaplotypeLength/m_sumNumTestcases)-mean_val*mean_val<<"\n";
  cout <<"numDoubleTestcases\t"<<m_sumNumDoubleTestcases<<"\n";
  cout <<"numBytesCopied\t"<<m_bytes_copied<<"\n";
  cout.flush();
}

void LoadTimeInitializer::debug_dump(string filename, string s, bool to_append, bool add_newline)
{
  map<string, ofstream*>::iterator mI = m_filename_to_fptr.find(filename);
  ofstream* fptr = 0;
  if(mI == m_filename_to_fptr.end())
  {
    m_filename_to_fptr[filename] = new ofstream();
    fptr = m_filename_to_fptr[filename];
    fptr->open(filename.c_str(), to_append ? ios::app : ios::out);
    assert(fptr->is_open());
  }
  else
    fptr = (*mI).second;
  //ofstream fptr;
  //fptr.open(filename.c_str(), to_append ? ofstream::app : ofstream::out);
  (*fptr) << s;
  if(add_newline)
    (*fptr) << "\n";
  //fptr.close();
}
void LoadTimeInitializer::debug_close()
{
  for(map<string,ofstream*>::iterator mB = m_filename_to_fptr.begin(), mE = m_filename_to_fptr.end();
      mB != mE;mB++)
  {
    (*mB).second->close();
    delete (*mB).second;
  }
  m_filename_to_fptr.clear();
}

void LoadTimeInitializer::dump_sandbox(testcase& tc)
{
  unsigned haplotypeLength = tc.haplen;
  unsigned readLength = tc.rslen;
  ofstream& dumpFptr = m_sandbox_fptr;
  for(unsigned k=0;k<haplotypeLength;++k)
    dumpFptr<<(char)(tc.hap[k]);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.rs[k]);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.q[k]+33);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.i[k]+33);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.d[k]+33);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.c[k]+33);
  dumpFptr<<"\n";
}
