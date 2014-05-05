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

