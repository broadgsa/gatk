CXX=g++
CXXFLAGS=-g -Wall -O2 -m64 -fPIC

.cpp.o:
	$(CXX) -c $(CXXFLAGS) -I$(BWA_HOME) -I$(JAVA_INCLUDE) $< -o $@

all: init org_broadinstitute_sting_alignment_bwa_BWACAligner.o bwa_gateway.o 
	$(LIBTOOL_COMMAND) org_broadinstitute_sting_alignment_bwa_BWACAligner.o bwa_gateway.o -o $(TARGET_LIB) -L$(BWA_HOME) -lbwacore $(EXTRA_LIBS)

init:
	@echo Please make sure the following platforms are set correctly on your machine.
	@echo BWA_HOME=$(BWA_HOME)
	@echo JAVA_INCLUDE=$(JAVA_INCLUDE)
	@echo TARGET_LIB=$(TARGET_LIB)
	@echo EXTRA_LIBS=$(EXTRA_LIBS)
	@echo LIBTOOL_COMMAND=$(LIBTOOL_COMMAND)x

clean:
	rm *.o libbwa.*
