#OMPCFLAGS=-fopenmp
#OMPLFLAGS=-fopenmp #-openmp-link static

#CFLAGS=-O2 -std=c++11 -W -Wall -march=corei7-avx -Wa,-q            -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas
#CFLAGS=-O2             -W -Wall -march=corei7 -mfpmath=sse -msse4.2 -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas

JAVA_ROOT=/opt/jdk1.7.0_25/
JNI_COMPILATION_FLAGS=-D_REENTRANT -fPIC -I${JAVA_ROOT}/include -I${JAVA_ROOT}/include/linux

COMMON_COMPILATION_FLAGS=$(JNI_COMPILATION_FLAGS) -O3 -W -Wall -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas
CC=icc
CXX=icc

LDFLAGS=-lm -lrt $(OMPLDFLAGS)

#USE_PAPI=1
PAPI_DIR=/home/karthikg/softwares/papi-5.3.0
ifdef USE_PAPI
  ifeq ($(USE_PAPI),1)
    COMMON_COMPILATION_FLAGS+=-I$(PAPI_DIR)/include
    LDFLAGS+=-L$(PAPI_DIR)/lib -lpapi
  endif
endif

BIN=libVectorLoglessPairHMM.so pairhmm-template-main checker
#BIN=checker

DEPDIR=.deps
DF=$(DEPDIR)/$(*).d

#Common across libJNI and sandbox
COMMON_SOURCES=utils.cc avx_function_instantiations.cc baseline.cc sse_function_instantiations.cc LoadTimeInitializer.cc
#Part of libJNI
LIBSOURCES=org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM.cc org_broadinstitute_sting_utils_pairhmm_DebugJNILoglessPairHMM.cc Sandbox.cc $(COMMON_SOURCES)
SOURCES=$(LIBSOURCES) pairhmm-template-main.cc pairhmm-1-base.cc
LIBOBJECTS=$(LIBSOURCES:.cc=.o)
COMMON_OBJECTS=$(COMMON_SOURCES:.cc=.o)


#No vectorization for these files
NO_VECTOR_SOURCES=org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM.cc org_broadinstitute_sting_utils_pairhmm_DebugJNILoglessPairHMM.cc pairhmm-template-main.cc pairhmm-1-base.cc utils.cc baseline.cc LoadTimeInitializer.cc Sandbox.cc
#Use -xAVX for these files
AVX_SOURCES=avx_function_instantiations.cc
#Use -xSSE4.2 for these files
SSE_SOURCES=sse_function_instantiations.cc

NO_VECTOR_OBJECTS=$(NO_VECTOR_SOURCES:.cc=.o)
AVX_OBJECTS=$(AVX_SOURCES:.cc=.o)
SSE_OBJECTS=$(SSE_SOURCES:.cc=.o)
$(NO_VECTOR_OBJECTS): CXXFLAGS=$(COMMON_COMPILATION_FLAGS)
$(AVX_OBJECTS): CXXFLAGS=$(COMMON_COMPILATION_FLAGS) -xAVX
$(SSE_OBJECTS): CXXFLAGS=$(COMMON_COMPILATION_FLAGS) -xSSE4.2
OBJECTS=$(NO_VECTOR_OBJECTS) $(AVX_OBJECTS) $(SSE_OBJECTS)

all: $(BIN) Sandbox.class

-include $(addprefix $(DEPDIR)/,$(SOURCES:.cc=.d))

checker: pairhmm-1-base.o $(COMMON_OBJECTS)
	$(CXX) $(OMPLFLAGS) -o $@ $^ $(LDFLAGS)

pairhmm-template-main:	pairhmm-template-main.o $(COMMON_OBJECTS)
	$(CXX) $(OMPLFLAGS) -o $@ $^ $(LDFLAGS)

libVectorLoglessPairHMM.so: $(LIBOBJECTS) 
	$(CXX) $(OMPLFLAGS) -shared -static-intel -o $@ $(LIBOBJECTS) ${LDFLAGS}


$(OBJECTS): %.o: %.cc
	@mkdir -p $(DEPDIR)
	$(CXX) -c -MMD -MF $(DF) $(CXXFLAGS) $(OUTPUT_OPTION) $<

Sandbox.class: Sandbox.java
	javac Sandbox.java

clean:
	rm -rf $(BIN) *.o $(DEPDIR) *.class
