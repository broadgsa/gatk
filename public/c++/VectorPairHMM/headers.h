#ifndef COMMON_HEADERS_H
#define COMMON_HEADERS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <ctype.h>

#include <sys/time.h>

#include <immintrin.h>
#include <emmintrin.h>
#include <omp.h>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <fenv.h>

#define STORE_FP_EXCEPTIONS(flagp, exceptions_array)                                                  \
  fegetexceptflag(&flagp, FE_OVERFLOW | FE_UNDERFLOW | FE_DIVBYZERO | FE_INVALID | __FE_DENORM);      \
  exceptions_array[FE_INVALID] += ((flagp & FE_INVALID));                               \
  exceptions_array[__FE_DENORM] += ((flagp & __FE_DENORM) >> 1);                        \
  exceptions_array[FE_DIVBYZERO] += ((flagp & FE_DIVBYZERO) >> 2);                      \
  exceptions_array[FE_OVERFLOW] += ((flagp & FE_OVERFLOW) >> 3);                        \
  exceptions_array[FE_UNDERFLOW] += ((flagp & FE_UNDERFLOW) >> 4);                      \
  feclearexcept(FE_ALL_EXCEPT);                 


#endif
