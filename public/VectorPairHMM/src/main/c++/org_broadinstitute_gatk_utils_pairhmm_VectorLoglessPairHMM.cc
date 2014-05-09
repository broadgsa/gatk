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
#include "jni_common.h"
#include "org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM.h"
//#include "template.h"
#include "utils.h"
#include "LoadTimeInitializer.h"

using namespace std;

JNIEXPORT jlong JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniGetMachineType
  (JNIEnv* env, jobject thisObject)
{
  return (jlong)get_machine_capabilities(); 
}

//Should be called only once for the whole Java process - initializes field ids for the classes JNIReadDataHolderClass
//and JNIHaplotypeDataHolderClass
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniInitializeClassFieldsAndMachineMask
  (JNIEnv* env, jobject thisObject, jclass readDataHolderClass, jclass haplotypeDataHolderClass, jlong mask)
{
  assert(readDataHolderClass);
  assert(haplotypeDataHolderClass);
  jfieldID fid;
  fid = env->GetFieldID(readDataHolderClass, "readBases", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for readBases");
  g_load_time_initializer.m_readBasesFID = fid;
  fid = env->GetFieldID(readDataHolderClass, "readQuals", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for readQuals");
  g_load_time_initializer.m_readQualsFID = fid;
  fid = env->GetFieldID(readDataHolderClass, "insertionGOP", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for insertionGOP");
  g_load_time_initializer.m_insertionGOPFID = fid;
  fid = env->GetFieldID(readDataHolderClass, "deletionGOP", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for deletionGOP");
  g_load_time_initializer.m_deletionGOPFID = fid;
  fid = env->GetFieldID(readDataHolderClass, "overallGCP", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for overallGCP");
  g_load_time_initializer.m_overallGCPFID = fid;

  fid = env->GetFieldID(haplotypeDataHolderClass, "haplotypeBases", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for haplotypeBases");
  g_load_time_initializer.m_haplotypeBasesFID = fid;
  if(mask != ENABLE_ALL_HARDWARE_FEATURES)
  {
    cout << "Using user supplied hardware mask to re-initialize function pointers for PairHMM\n";
    initialize_function_pointers((uint64_t)mask);
    cout.flush();
  }
}

JNIEXPORT void JNICALL initializeHaplotypes
  (JNIEnv * env, jobject& thisObject, jint numHaplotypes, jobjectArray& haplotypeDataArray,
   vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector, vector<unsigned>& haplotypeBasesLengths)
{
  jboolean is_copy = JNI_FALSE;
  haplotypeBasesArrayVector.clear();
  haplotypeBasesLengths.clear();
  haplotypeBasesArrayVector.resize(numHaplotypes);
  haplotypeBasesLengths.resize(numHaplotypes);
  jsize haplotypeBasesLength = 0;
  for(unsigned j=0;j<numHaplotypes;++j)
  {
    jobject haplotypeObject = env->GetObjectArrayElement(haplotypeDataArray, j);
    jbyteArray haplotypeBases = (jbyteArray)env->GetObjectField(haplotypeObject, g_load_time_initializer.m_haplotypeBasesFID);
#ifdef ENABLE_ASSERTIONS
    assert(haplotypeBases && ("haplotypeBases is NULL at index : "+to_string(j)+"\n").c_str());
#endif
    //Need a global reference as this will be accessed across multiple JNI calls to JNIComputeLikelihoods()
    jbyteArray haplotypeBasesGlobalRef = (jbyteArray)env->NewGlobalRef(haplotypeBases);
#ifdef ENABLE_ASSERTIONS
    assert(haplotypeBasesGlobalRef && ("Could not get global ref to haplotypeBases at index : "+to_string(j)+"\n").c_str());
#endif
    env->DeleteLocalRef(haplotypeBases);	//free the local reference
    jbyte* haplotypeBasesArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(haplotypeBasesGlobalRef, &is_copy);
    haplotypeBasesLength = env->GetArrayLength(haplotypeBasesGlobalRef);
#ifdef ENABLE_ASSERTIONS
    assert(haplotypeBasesArray && "haplotypeBasesArray not initialized in JNI"); 
    //assert(haplotypeBasesLength < MCOLS);
#endif
#ifdef DEBUG0_1
    cout << "JNI haplotype length "<<haplotypeBasesLength<<"\n";
#endif
    haplotypeBasesArrayVector[j] = make_pair(haplotypeBasesGlobalRef, haplotypeBasesArray);
    haplotypeBasesLengths[j] = haplotypeBasesLength;
#ifdef DEBUG3
    for(unsigned k=0;k<haplotypeBasesLength;++k)
      g_load_time_initializer.debug_dump("haplotype_bases_jni.txt",to_string((int)haplotypeBasesArray[k]),true);
#endif
#ifdef DO_PROFILING
    g_load_time_initializer.update_stat(HAPLOTYPE_LENGTH_IDX, haplotypeBasesLength);
    g_load_time_initializer.m_bytes_copied += (is_copy ? haplotypeBasesLength : 0);
#endif
  }
}

JNIEXPORT void JNICALL releaseHaplotypes(JNIEnv * env, jobject thisObject,
    vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector, vector<unsigned>& haplotypeBasesLengths
    )
{
  //Now release haplotype arrays
  for(int j=haplotypeBasesArrayVector.size()-1;j>=0;--j)	//note the order - reverse of GET
  {
    RELEASE_BYTE_ARRAY_ELEMENTS(haplotypeBasesArrayVector[j].first, haplotypeBasesArrayVector[j].second, JNI_RO_RELEASE_MODE);
    env->DeleteGlobalRef(haplotypeBasesArrayVector[j].first);	//free the global reference
  }
  haplotypeBasesArrayVector.clear();
  haplotypeBasesLengths.clear(); 
}


vector<pair<jbyteArray, jbyte*> > g_haplotypeBasesArrayVector;
vector<unsigned> g_haplotypeBasesLengths;
//Since the list of haplotypes against which the reads are evaluated in PairHMM is the same for a region,
//transfer the list only once
//Works only for ST case as the haplotype data is stored in global variables
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniInitializeHaplotypes
  (JNIEnv * env, jobject thisObject, jint numHaplotypes, jobjectArray haplotypeDataArray)
{
#ifdef SINGLE_THREADED_ONLY
  //To ensure, GET_BYTE_ARRAY_ELEMENTS is invoked only once for each haplotype, store bytearrays in a vector
  initializeHaplotypes(env, thisObject, numHaplotypes, haplotypeDataArray, g_haplotypeBasesArrayVector, g_haplotypeBasesLengths);
#endif
}


//Create a vector of testcases for computation - copy the references to bytearrays read/readQuals etc into the appropriate
//testcase struct
inline JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniInitializeTestcasesVector
  (JNIEnv* env, jint numReads, jint numHaplotypes, jobjectArray& readDataArray,
   vector<vector<pair<jbyteArray,jbyte*> > >& readBasesArrayVector,
   vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector, vector<unsigned>& haplotypeBasesLengths,
   vector<testcase>& tc_array)
{
  jboolean is_copy = JNI_FALSE;
  unsigned tc_idx = 0;
  for(unsigned i=0;i<numReads;++i)
  {
    //Get bytearray fields from read
    jobject readObject = env->GetObjectArrayElement(readDataArray, i);
    jbyteArray readBases = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_readBasesFID);
    jbyteArray insertionGOP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_insertionGOPFID);
    jbyteArray deletionGOP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_deletionGOPFID);
    jbyteArray overallGCP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_overallGCPFID);
    jbyteArray readQuals = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_readQualsFID);

#ifdef ENABLE_ASSERTIONS
    assert(readBases && ("readBases is NULL at index : "+to_string(i)+"\n").c_str());
    assert(insertionGOP && ("insertionGOP is NULL at index : "+to_string(i)+"\n").c_str());
    assert(deletionGOP && ("deletionGOP is NULL at index : "+to_string(i)+"\n").c_str());
    assert(overallGCP && ("overallGCP is NULL at index : "+to_string(i)+"\n").c_str());
    assert(readQuals && ("readQuals is NULL at index : "+to_string(i)+"\n").c_str());
#endif
    jsize readLength = env->GetArrayLength(readBases);

    jbyte* readBasesArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(readBases, &is_copy);	//order of GET-RELEASE is important
    jbyte* readQualsArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(readQuals, &is_copy);
    jbyte* insertionGOPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(insertionGOP, &is_copy);
    jbyte* deletionGOPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(deletionGOP, &is_copy);
    jbyte* overallGCPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(overallGCP, &is_copy);
#ifdef DO_PROFILING
    g_load_time_initializer.m_bytes_copied += (is_copy ? readLength*5 : 0);
    g_load_time_initializer.update_stat(READ_LENGTH_IDX, readLength);
#endif
#ifdef ENABLE_ASSERTIONS
    assert(readBasesArray && "readBasesArray not initialized in JNI"); 
    assert(readQualsArray && "readQualsArray not initialized in JNI"); 
    assert(insertionGOPArray && "insertionGOP array not initialized in JNI");
    assert(deletionGOPArray && "deletionGOP array not initialized in JNI");
    assert(overallGCPArray && "overallGCP array not initialized in JNI");
    //assert(readLength < MROWS); 
    assert(readLength == env->GetArrayLength(readQuals));
    assert(readLength == env->GetArrayLength(insertionGOP));
    assert(readLength == env->GetArrayLength(deletionGOP));
    assert(readLength == env->GetArrayLength(overallGCP));
#endif
#ifdef DEBUG0_1
    cout << "JNI read length "<<readLength<<"\n";
#endif
#ifdef DEBUG3
    for(unsigned j=0;j<readLength;++j)
    {
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)readBasesArray[j]),true);
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)readQualsArray[j]),true);
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)insertionGOPArray[j]),true);
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)deletionGOPArray[j]),true);
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)overallGCPArray[j]),true);
    }
#endif
    for(unsigned j=0;j<numHaplotypes;++j)
    {
      jsize haplotypeLength = (jsize)haplotypeBasesLengths[j];
      jbyte* haplotypeBasesArray = haplotypeBasesArrayVector[j].second;
      tc_array[tc_idx].rslen = (int)readLength;
      tc_array[tc_idx].haplen = (int)haplotypeLength;
      tc_array[tc_idx].hap = (char*)haplotypeBasesArray;
      tc_array[tc_idx].rs = (char*)readBasesArray;
      tc_array[tc_idx].q = (char*)readQualsArray;
      tc_array[tc_idx].i = (char*)insertionGOPArray;
      tc_array[tc_idx].d = (char*)deletionGOPArray;
      tc_array[tc_idx].c = (char*)overallGCPArray;
#ifdef DO_PROFILING
      g_load_time_initializer.update_stat(PRODUCT_READ_LENGTH_HAPLOTYPE_LENGTH_IDX, ((uint64_t)readLength)*((uint64_t)haplotypeLength));
#endif
#ifdef DUMP_TO_SANDBOX
      g_load_time_initializer.dump_sandbox(tc_array[tc_idx], tc_idx, numReads, numHaplotypes);
#endif
      ++tc_idx;  
    }
    //Store the read array references and release them at the end because they are used by compute_full_prob
    //Maintain order in which GET_BYTE_ARRAY_ELEMENTS called
    readBasesArrayVector[i].clear();
    readBasesArrayVector[i].resize(5);
    readBasesArrayVector[i][0] = make_pair(readBases, readBasesArray);
    readBasesArrayVector[i][1] = make_pair(readQuals, readQualsArray);
    readBasesArrayVector[i][2] = make_pair(insertionGOP, insertionGOPArray);
    readBasesArrayVector[i][3] = make_pair(deletionGOP, deletionGOPArray);
    readBasesArrayVector[i][4] = make_pair(overallGCP, overallGCPArray);
  }
}

//Do compute over vector of testcase structs
inline void compute_testcases(vector<testcase>& tc_array, unsigned numTestCases, double* likelihoodDoubleArray,
    unsigned maxNumThreadsToUse)
{
#ifdef DO_REPEAT_PROFILING
  for(unsigned i=0;i<10;++i)
#endif
  {
#pragma omp parallel for schedule (dynamic,10000) num_threads(maxNumThreadsToUse)
    for(unsigned tc_idx=0;tc_idx<numTestCases;++tc_idx)
    {
      float result_avxf = g_compute_full_prob_float(&(tc_array[tc_idx]), 0);
      double result = 0;
      if (result_avxf < MIN_ACCEPTED) {
        double result_avxd = g_compute_full_prob_double(&(tc_array[tc_idx]), 0);
        result = log10(result_avxd) - log10(ldexp(1.0, 1020.0));
#ifdef DO_PROFILING
        g_load_time_initializer.update_stat(NUM_DOUBLE_INVOCATIONS_IDX, 1);
#endif
      }
      else
        result = (double)(log10f(result_avxf) - log10f(ldexpf(1.f, 120.f)));
      likelihoodDoubleArray[tc_idx] = result;
    }
  }
}

//Inform the Java VM that we no longer need access to the read arrays (and free memory)
inline JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniReleaseReadArrays
  (JNIEnv* env, vector<vector<pair<jbyteArray,jbyte*> > >& readBasesArrayVector)
{
  //Release read arrays first
  for(int i=readBasesArrayVector.size()-1;i>=0;--i)//note the order - reverse of GET
  {
    for(int j=readBasesArrayVector[i].size()-1;j>=0;--j)
      RELEASE_BYTE_ARRAY_ELEMENTS(readBasesArrayVector[i][j].first, readBasesArrayVector[i][j].second, JNI_RO_RELEASE_MODE);
    readBasesArrayVector[i].clear();
  }
  readBasesArrayVector.clear();
}


#ifdef DO_WARMUP
uint64_t g_sum = 0;
#endif
//JNI function to invoke compute_full_prob_avx
//readDataArray - array of JNIReadDataHolderClass objects which contain the readBases, readQuals etc
//haplotypeDataArray - array of JNIHaplotypeDataHolderClass objects which contain the haplotypeBases
//likelihoodArray - array of doubles to return results back to Java. Memory allocated by Java prior to JNI call
//maxNumThreadsToUse - Max number of threads that OpenMP can use for the HMM computation
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniComputeLikelihoods
  (JNIEnv* env, jobject thisObject, jint numReads, jint numHaplotypes, 
   jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray, jint maxNumThreadsToUse)
{
#ifdef DEBUG0_1
  cout << "JNI numReads "<<numReads<<" numHaplotypes "<<numHaplotypes<<"\n";
#endif
  jboolean is_copy = JNI_FALSE;
  struct timespec start_time;
  unsigned numTestCases = numReads*numHaplotypes;
  //vector to store testcases
  vector<testcase> tc_array;
  tc_array.clear();
  tc_array.resize(numTestCases);
  //Store read arrays for release later
  vector<vector<pair<jbyteArray,jbyte*> > > readBasesArrayVector;
  readBasesArrayVector.clear();
  readBasesArrayVector.resize(numReads);
#ifdef DUMP_TO_SANDBOX
  g_load_time_initializer.open_sandbox();
#endif
#ifdef DO_PROFILING
  get_time(&start_time);
#endif

#ifdef SINGLE_THREADED_ONLY
  vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector = g_haplotypeBasesArrayVector;
  vector<unsigned>& haplotypeBasesLengths = g_haplotypeBasesLengths;
#else
  vector<pair<jbyteArray, jbyte*> > l_haplotypeBasesArrayVector;
  vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector = l_haplotypeBasesArrayVector;
  vector<unsigned> l_haplotypeBasesLengths;
  vector<unsigned>& haplotypeBasesLengths = l_haplotypeBasesLengths;
  initializeHaplotypes(env, thisObject, numHaplotypes, haplotypeDataArray, haplotypeBasesArrayVector, haplotypeBasesLengths);
#endif
  //Copy byte array references from Java memory into vector of testcase structs
  Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniInitializeTestcasesVector(env,
      numReads, numHaplotypes, readDataArray, readBasesArrayVector, haplotypeBasesArrayVector, haplotypeBasesLengths, tc_array);

#ifdef DO_PROFILING
  g_load_time_initializer.m_data_transfer_time += diff_time(start_time);
#endif

  //Get double array where results are stored (to pass back to java)
  jdouble* likelihoodDoubleArray = (jdouble*)GET_DOUBLE_ARRAY_ELEMENTS(likelihoodArray, &is_copy);
#ifdef ENABLE_ASSERTIONS
  assert(likelihoodDoubleArray && "likelihoodArray is NULL");
  assert(env->GetArrayLength(likelihoodArray) == numTestCases);
#endif
#ifdef DO_WARMUP        //ignore - only for crazy profiling
  for(unsigned i=0;i<haplotypeBasesArrayVector.size();++i)
  {
    unsigned curr_size = env->GetArrayLength(haplotypeBasesArrayVector[i].first);
    for(unsigned j=0;j<curr_size;++j)
      g_sum += ((uint64_t)((haplotypeBasesArrayVector[i].second)[j]));
  }
  for(unsigned i=0;i<readBasesArrayVector.size();++i)
  {
    for(unsigned j=0;j<readBasesArrayVector[i].size();++j)
    {
      unsigned curr_size = env->GetArrayLength(readBasesArrayVector[i][j].first);
      for(unsigned k=0;k<curr_size;++k)
        g_sum += ((uint64_t)((readBasesArrayVector[i][j].second)[k]));
    }
  }
#endif
#ifdef DO_PROFILING
  g_load_time_initializer.m_bytes_copied += (is_copy ? numTestCases*sizeof(double) : 0);
  get_time(&start_time);
#endif
  compute_testcases(tc_array, numTestCases, likelihoodDoubleArray, maxNumThreadsToUse); //actual computation
#ifdef DO_PROFILING
  g_load_time_initializer.m_compute_time += diff_time(start_time);
#endif
#ifdef DUMP_COMPUTE_VALUES
  for(unsigned tc_idx=0;tc_idx<numTestCases;++tc_idx)
    g_load_time_initializer.debug_dump("return_values_jni.txt",to_string(likelihoodDoubleArray[tc_idx]),true);
#endif
#ifdef DO_PROFILING
  get_time(&start_time);
#endif
  RELEASE_DOUBLE_ARRAY_ELEMENTS(likelihoodArray, likelihoodDoubleArray, 0); //release mode 0, copy back results to Java memory (if copy made)
  Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniReleaseReadArrays(env, readBasesArrayVector);
#ifndef SINGLE_THREADED_ONLY
  releaseHaplotypes(env, thisObject, haplotypeBasesArrayVector, haplotypeBasesLengths);
#endif

#ifdef DO_PROFILING
  g_load_time_initializer.m_data_transfer_time += diff_time(start_time);
  g_load_time_initializer.update_stat(NUM_REGIONS_IDX, 1);
  g_load_time_initializer.update_stat(NUM_READS_IDX, numReads);
  g_load_time_initializer.update_stat(NUM_HAPLOTYPES_IDX, numHaplotypes);
  g_load_time_initializer.update_stat(NUM_TESTCASES_IDX, numTestCases);
#endif
  tc_array.clear();
#ifdef DUMP_TO_SANDBOX
  g_load_time_initializer.close_sandbox();
#endif
}

//If single threaded, release haplotypes at the end of a region
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniFinalizeRegion
  (JNIEnv * env, jobject thisObject)
{
#ifdef SINGLE_THREADED_ONLY
  releaseHaplotypes(env, thisObject, g_haplotypeBasesArrayVector, g_haplotypeBasesLengths);
#endif
}


JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniClose
  (JNIEnv* env, jobject thisObject)
{
#ifdef DO_PROFILING
  g_load_time_initializer.print_profiling();
#endif
#ifdef DUMP_COMPUTE_VALUES
  g_load_time_initializer.debug_close();
#endif
}

