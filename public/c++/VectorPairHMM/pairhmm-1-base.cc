//#define DEBUG 1
//#define DEBUG0_1 1
//#define DEBUG3 1
#include "headers.h"
#include "template.h"
#include "utils.h"
#include "LoadTimeInitializer.h"

using namespace std;

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
  unsigned chunk_size = 10000;
  if(argc >= 4)
    chunk_size = strtol(argv[3],0,10);

  do_compute(argv[1], use_old_read_testcase, chunk_size); 

  return 0;  
}

