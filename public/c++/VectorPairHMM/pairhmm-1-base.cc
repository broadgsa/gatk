//#define DEBUG 1
//#define DEBUG0_1 1
//#define DEBUG3 1
#include "headers.h"
#include "utils.h"
#include "LoadTimeInitializer.h"
using namespace std;

int main(int argc, char** argv)
{
#define BATCH_SIZE 10000
  if(argc < 2)
  {
    cerr << "Needs path to input file as argument\n";
    exit(0);
  }
  bool use_old_read_testcase = false;
  if(argc >= 3 && string(argv[2]) == "1")
    use_old_read_testcase = true;
  unsigned chunk_size = 10000;
  bool do_check = true;
  uint64_t mask = ~(0ull);
  for(int i=3;i<argc;++i)
  {
    if(strncmp(argv[i], "-chunk_size", 15) == 0)
    {
      ++i;
      chunk_size = strtol(argv[i],0,10);
    }
    else
      if(strncmp(argv[i], "-mask", 15) == 0)
      {
        ++i;
        mask = strtoll(argv[i],0,16);
      }
      else
        if(strncmp(argv[i], "-no-check", 15) == 0)
          do_check = false;
  }
  if(mask != (~0ull))
    initialize_function_pointers(mask);
  do_compute(argv[1], use_old_read_testcase, chunk_size, do_check); 
  return 0;  
}

