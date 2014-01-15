OMPCFLAGS=-fopenmp
#OMPLDFLAGS=-lgomp

#CFLAGS=-O2 -std=c++11 -W -Wall -march=corei7-avx -Wa,-q            -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas
#CFLAGS=-O2             -W -Wall -march=corei7 -mfpmath=sse -msse4.2 -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas

JAVA_ROOT=/opt/jdk1.7.0_25/
JNI_COMPILATION_FLAGS=-D_REENTRANT -fPIC -I${JAVA_ROOT}/include -I${JAVA_ROOT}/include/linux

CFLAGS=-O3 -W -Wall -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas -xAVX

CXXFLAGS=$(CFLAGS)
CC=icc
CXX=icc

LDFLAGS=-lm $(OMPLDFLAGS)

#BIN:=pairhmm-1-base #pairhmm-2-omp pairhmm-3-hybrid-float-double pairhmm-4-hybrid-diagonal pairhmm-5-hybrid-diagonal-homogeneus pairhmm-6-onlythreediags pairhmm-7-presse pairhmm-8-sse #pairhmm-dev
BIN:=libJNILoglessPairHMM.so pairhmm-template-main

#SOURCES=pairhmm-1-base.cc input.cc 
LIBSOURCES=org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM.cc hmm_mask.cc
SOURCES=$(LIBSOURCES) pairhmm-template-main.cc
LIBOBJECTS=$(LIBSOURCES:.cc=.o)
DEPDIR=.deps
DF=$(DEPDIR)/$(*).d

all: $(BIN)

-include $(addprefix $(DEPDIR)/,$(SOURCES:.cc=.d))

pairhmm-template-main:	pairhmm-template-main.o hmm_mask.o
	$(CXX) -fopenmp -o $@ $^ $(LDFLAGS)

libJNILoglessPairHMM.so: $(LIBOBJECTS) 
	$(CXX) -shared -o $@ $(LIBOBJECTS) 

%.o: %.cc
	@mkdir -p $(DEPDIR)
	$(COMPILE.cpp) -MMD -MF $(DF) $(JNI_COMPILATION_FLAGS) $(CXXFLAGS) $(OUTPUT_OPTION) $<


clean:
	rm -f $(BIN) *.o
