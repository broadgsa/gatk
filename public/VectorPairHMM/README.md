Implementation overview:
Created a new Java class called VectorLoglessPairHMM which extends LoglessPairHMM and 
overrides functions from both LoglessPairHMM and PairHMM.
1. Constructor: Call base class constructors. Then, load the native library located in this 
directory and call an init function (with suffix 'jniInitializeClassFieldsAndMachineMask') in the 
library to determine fields ids for the members of classes JNIReadDataHolder and 
JNIHaplotypeDataHolders. The native code stores the field ids (struct offsets) for the classes and 
re-uses them for subsequent computations. Optionally, the user can disable the vector 
implementation, by using the 'mask' argument (see comments for a more detailed explanation).
2. When the library is loaded, it invokes the constructor of the class LoadTimeInitializer (because 
a global variable g_load_time_initializer is declared in the library).  This constructor 
(LoadTimeInitializer.cc) can be used to perform various initializations.  Currently, it initializes 
two global function pointers to point to the function implementation that is supported on the 
machine (AVX/SSE/un-vectorized) on which the program is being run. The two pointers are for float 
and double respectively.  The global function pointers are declared in utils.cc and are assigned in 
the function initialize_function_pointers() defined in utils.cc and invoked from the constructor of 
LoadTimeInitializer.
Other initializations in LoadTimeInitializer:
* ConvertChar::init - sets some masks for the vector implementation
* FTZ for performance
* stat counters = 0
* debug structs (which are never used in non-debug mode)
This initialization is done only once for the whole program.
3. initialize(): To initialize the region for PairHMM. Pass haplotype bases to native code through 
the JNIHaplotypeDataHolder class.  Since the haplotype list is common across multiple samples in 
computeReadLikelihoods(), we can pass the haplotype bases to the native code once and re-use across 
multiple samples.
4. computeLikelihoods(): Copies array references for readBases/quals etc to array of 
JNIReadDataHolder objects.  Invokes the JNI function to perform the computation and updates the 
likelihoodMap.
The JNI function copies the byte array references into an array of testcase structs and invokes the 
compute_full_prob function through the function pointers initialized earlier.
The primary native function called is 
Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniComputeLikelihoods. It uses 
standard JNI calls to get and return data from/to the Java class VectorLoglessPairHMM.  The last 
argument to the function is the maximum number of OpenMP threads to use while computing PairHMM in 
C++. This option is set when the native function call is made from JNILoglessPairHMM 
computeLikelihoods - currently it is set to 12 (no logical reason).
Note: OpenMP has been disabled for now - insufficient #testcases per call to computeLikelihoods() to 
justify multi-threading.
5. finalizeRegion(): Releases the haplotype arrays initialized in step 3 - should be called at the 
end of every region (line 351 in PairHMMLikelihoodCalculationEngine).

Note: Debug code has been moved to a separate class DebugJNILoglessPairHMM.java.

Compiling:
The native library (called libVectorLoglessPairHMM.so) can be compiled with icc (Intel C compiler) 
or gcc versions >= 4.8.1 that support AVX intrinsics. By default, the make process tries to invoke 
icc. To use gcc, edit the file 'pom.xml' (in this directory) and enable the environment variables 
USE_GCC,C_COMPILER and CPP_COMPILER (edit and uncomment lines 60-62).
Using Maven:
Type 'mvn install' in this directory - this will build the library (by invoking 'make') and copy the 
native library to the directory 
${sting-utils.basedir}/src/main/resources/org/broadinstitute/sting/utils/pairhmm
The GATK maven build process (when run) will bundle the library into the StingUtils jar file from 
the copied directory.
Simple build:
cd src/main/c++
make  

Running:
The default implementation of PairHMM is now VECTOR_LOGLESS_CACHING in HaplotypeCaller.java. To use 
the Java version, use the command line argument "--pair_hmm_implementation LOGLESS_CACHING". (see 
run.sh in src/main/c++).
The native library is bundled with the StingUtils jar file. When HaplotypeCaller is invoked, then 
the library is unpacked from the jar file, copied to the /tmp directory (with a unique id) and 
loaded by the Java class VectorLoglessPairHMM in the constructor (if it has not been loaded 
already).
The default library can be overridden by using the -Djava.library.path argument (see 
src/main/c++/run.sh for an example) for the JVM to pass the path to the library. If the library 
libVectorLoglessPairHMM.so can be found in java.library.path, then it is loaded and the 'packed' 
library is not used.
