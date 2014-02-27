Implementation overview:
Created a new Java class called VectorLoglessPairHMM which extends LoglessPairHMM and 
overrides functions from both LoglessPairHMM and PairHMM.
1. Constructor: Call base class constructors. Then, load the native library located in this 
directory and call a global init function in the library to determine fields ids for the 
members of classes JNIReadDataHolder and JNIHaplotypeDataHolders.
2. When the library is loaded, it initializes two global function pointers to point to the
function implementation that is supported on the machine on which the program is being 
run. The two pointers are for float and double respectively. This initialization is done 
only once for the whole program.
3. initialize(): To initialized the region for PairHMM. Pass haplotype bases to native 
code through the JNIHaplotypeDataHolders class.  Since the haplotype list is common across multiple 
samples in computeReadLikelihoods(), we can store the haplotype bases to the native code once and 
re-use across multiple samples.
4. computeLikelihoods(): Copies array references for readBases/quals etc to array of 
JNIReadDataHolder objects.  Invokes the JNI function to perform the computation and 
updates the likelihoodMap.

Note: Debug code has been moved to a separate class DebugJNILoglessPairHMM.java.

On the C++ side, the primary function called is 
Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniComputeLikelihoods. It 
uses standard JNI calls to get and return data from/to the Java class 
VectorLoglessPairHMM.  The last argument to the function is the maximum number of OpenMP 
threads to use while computing PairHMM in C++. This option is set when the native function 
call is made from JNILoglessPairHMM computeLikelihoods - currently it is set to 12 (no 
logical reason).
Note: OpenMP has been disabled for now.

Compiling:
Make sure you have icc (Intel C compiler) available. Currently, gcc does not seem to 
support all AVX intrinsics.
Type 'make'. This should create a library called libVectorLoglessPairHMM.so 

Running:
The default implementation of PairHMM is still LOGLESS_CACHING in HaplotypeCaller.java. To use the 
native version, use the command line argument "--pair_hmm_implementation VECTOR_LOGLESS_CACHING"
(see run.sh in src/main/c++).
The native library is bundled with the StingUtils jar file. When HaplotypeCaller is invoked with the 
VectorLoglessPairHMM implementation (see run.sh in the directory src/main/c++), then the library is 
unpacked from the jar file, copied to the /tmp directory (with a unique id) and loaded by the Java class 
VectorLoglessPairHMM in the constructor (if it has not been loaded already).
The default library can be overridden by using the -Djava.library.path argument for the JVM to pass 
the path to the library. If the library libVectorLoglessPairHMM.so can be found in 
java.library.path, then it is loaded and the 'packed' library is not used.
See run.sh in this directory on how to invoke HaplotypeCaller with the vector implementation of 
PairHMM. The argument -Djava.library.path is needed if you wish to override the default packed 
library, else unnecessary.
