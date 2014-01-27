#include "headers.h"
#include "jni_common.h"
#include "org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM.h"
#include "template.h"
#include "utils.h"
#include "LoadTimeInitializer.h"

using namespace std;

//Should be called only once for the whole Java process - initializes field ids for the classes JNIReadDataHolderClass
//and JNIHaplotypeDataHolderClass
JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniGlobalInit
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
}

//Since the list of haplotypes against which the reads are evaluated in PairHMM is the same for a region,
//transfer the list only once
vector<pair<jbyteArray, jbyte*> > g_haplotypeBasesArrayVector;
vector<unsigned> g_haplotypeBasesLengths;
JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniInitializeHaplotypes
  (JNIEnv * env, jobject thisObject, jint numHaplotypes, jobjectArray haplotypeDataArray)
{
  jboolean is_copy = JNI_FALSE;
  //To ensure, GET_BYTE_ARRAY_ELEMENTS is invoked only once for each haplotype, store bytearrays in a vector
  vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector = g_haplotypeBasesArrayVector;
  haplotypeBasesArrayVector.clear();
  g_haplotypeBasesLengths.clear();
  haplotypeBasesArrayVector.resize(numHaplotypes);
  g_haplotypeBasesLengths.resize(numHaplotypes);
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
    assert(haplotypeBasesLength < MCOLS);
#endif
#ifdef DEBUG0_1
    cout << "JNI haplotype length "<<haplotypeBasesLength<<"\n";
#endif
    haplotypeBasesArrayVector[j] = make_pair(haplotypeBasesGlobalRef, haplotypeBasesArray);
    g_haplotypeBasesLengths[j] = haplotypeBasesLength;
#ifdef DEBUG3
    for(unsigned k=0;k<haplotypeBasesLength;++k)
      g_load_time_initializer.debug_dump("haplotype_bases_jni.txt",to_string((int)haplotypeBasesArray[k]),true);
#endif
#ifdef DO_PROFILING
    g_load_time_initializer.m_sumHaplotypeLengths += haplotypeBasesLength;
#endif
  }
}

//JNI function to invoke compute_full_prob_avx
//readDataArray - array of JNIReadDataHolderClass objects which contain the readBases, readQuals etc
//haplotypeDataArray - array of JNIHaplotypeDataHolderClass objects which contain the haplotypeBases
//likelihoodArray - array of doubles to return results back to Java. Memory allocated by Java prior to JNI call
//maxNumThreadsToUse - Max number of threads that OpenMP can use for the HMM computation
JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniComputeLikelihoods
  (JNIEnv* env, jobject thisObject, jint numReads, jint numHaplotypes, 
   jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray, jint maxNumThreadsToUse)
{
#ifdef DEBUG0_1
  cout << "JNI numReads "<<numReads<<" numHaplotypes "<<numHaplotypes<<"\n";
#endif
  double start_time = 0;
  //haplotype vector from earlier store - note the reference to vector, not copying
  vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector = g_haplotypeBasesArrayVector;
  jboolean is_copy = JNI_FALSE;

  unsigned numTestCases = numReads*numHaplotypes;
  //vector to store results
  vector<testcase> tc_array;
  tc_array.clear();
  tc_array.resize(numTestCases);
  unsigned tc_idx = 0;
  //Store arrays for release later
  vector<vector<pair<jbyteArray,jbyte*> > > readBasesArrayVector;
  readBasesArrayVector.clear();
  readBasesArrayVector.resize(numReads);
#ifdef DO_PROFILING
  start_time = get_time();
#endif
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
#ifdef ENABLE_ASSERTIONS
    assert(readBasesArray && "readBasesArray not initialized in JNI"); 
    assert(readQualsArray && "readQualsArray not initialized in JNI"); 
    assert(insertionGOPArray && "insertionGOP array not initialized in JNI");
    assert(deletionGOPArray && "deletionGOP array not initialized in JNI");
    assert(overallGCPArray && "overallGCP array not initialized in JNI");
    assert(readLength < MROWS); 
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
      jsize haplotypeLength = (jsize)g_haplotypeBasesLengths[j];
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
      g_load_time_initializer.m_sumProductReadLengthHaplotypeLength += (readLength*haplotypeLength);
      g_load_time_initializer.m_sumSquareProductReadLengthHaplotypeLength += ((readLength*haplotypeLength)*(readLength*haplotypeLength));
#endif
      ++tc_idx;  
    }
    //Release read arrays at end because they are used by compute_full_prob
    //Maintain order in which GET_BYTE_ARRAY_ELEMENTS called
    readBasesArrayVector[i].clear();
    readBasesArrayVector[i].resize(5);
    readBasesArrayVector[i][0] = make_pair(readBases, readBasesArray);
    readBasesArrayVector[i][1] = make_pair(readQuals, readQualsArray);
    readBasesArrayVector[i][2] = make_pair(insertionGOP, insertionGOPArray);
    readBasesArrayVector[i][3] = make_pair(deletionGOP, deletionGOPArray);
    readBasesArrayVector[i][4] = make_pair(overallGCP, overallGCPArray);
#ifdef DO_PROFILING
    g_load_time_initializer.m_sumReadLengths += readLength;
#endif
  }
#ifdef DO_PROFILING
  g_load_time_initializer.m_data_transfer_time += get_time();
#endif

  jdouble* likelihoodDoubleArray = (jdouble*)GET_DOUBLE_ARRAY_ELEMENTS(likelihoodArray, &is_copy);
#ifdef ENABLE_ASSERTIONS
  assert(likelihoodDoubleArray && "likelihoodArray is NULL");
  assert(env->GetArrayLength(likelihoodArray) == numTestCases);
#endif
#ifdef DO_PROFILING
  start_time = get_time();
#endif
#pragma omp parallel for schedule (dynamic,10) private(tc_idx) num_threads(maxNumThreadsToUse) 
  for(tc_idx=0;tc_idx<numTestCases;++tc_idx)
  {
    float result_avxf = g_compute_full_prob_float(&(tc_array[tc_idx]), 0);
    double result = 0;
    if (result_avxf < MIN_ACCEPTED) {
      double result_avxd = g_compute_full_prob_double(&(tc_array[tc_idx]), 0);
      result = log10(result_avxd) - log10(ldexp(1.0, 1020.0));
    }
    else
      result = (double)(log10f(result_avxf) - log10f(ldexpf(1.f, 120.f)));

    likelihoodDoubleArray[tc_idx] = result;
  }
#ifdef DO_PROFILING
  g_load_time_initializer.m_compute_time += get_time();
#endif
#ifdef DEBUG
  for(tc_idx=0;tc_idx<numTestCases;++tc_idx)
  {
    g_load_time_initializer.debug_dump("return_values_jni.txt",to_string(likelihoodDoubleArray[tc_idx]),true);
  }
#endif
#ifdef DO_PROFILING
  start_time = get_time();
#endif
  RELEASE_DOUBLE_ARRAY_ELEMENTS(likelihoodArray, likelihoodDoubleArray, 0); //release mode 0, copy back results to Java memory
  
  //Release read arrays first
  for(int i=readBasesArrayVector.size()-1;i>=0;--i)//note the order - reverse of GET
  {
    for(int j=readBasesArrayVector[i].size()-1;j>=0;--j)
      RELEASE_BYTE_ARRAY_ELEMENTS(readBasesArrayVector[i][j].first, readBasesArrayVector[i][j].second, JNI_RO_RELEASE_MODE);
    readBasesArrayVector[i].clear();
  }
  readBasesArrayVector.clear();
#ifdef DO_PROFILING
  g_load_time_initializer.m_data_transfer_time += get_time();
#endif
  tc_array.clear();
#ifdef DO_PROFILING
  g_load_time_initializer.m_sumNumReads += numReads;
  g_load_time_initializer.m_sumSquareNumReads += numReads*numReads;
  g_load_time_initializer.m_sumNumHaplotypes += numHaplotypes;
  g_load_time_initializer.m_sumSquareNumHaplotypes += numHaplotypes*numHaplotypes;
  g_load_time_initializer.m_sumNumTestcases += numTestCases;
  g_load_time_initializer.m_sumSquareNumTestcases += numTestCases*numTestCases;
  g_load_time_initializer.m_maxNumTestcases = numTestCases > g_load_time_initializer.m_maxNumTestcases ? numTestCases
    : g_load_time_initializer.m_maxNumTestcases;
  ++(g_load_time_initializer.m_num_invocations);
#endif
#ifdef DEBUG
  g_load_time_initializer.debug_close();
#endif
}

//Release haplotypes at the end of a region
JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniFinalizeRegion
  (JNIEnv * env, jobject thisObject)
{
  vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector = g_haplotypeBasesArrayVector;
  //Now release haplotype arrays
  for(int j=haplotypeBasesArrayVector.size()-1;j>=0;--j)	//note the order - reverse of GET
  {
    RELEASE_BYTE_ARRAY_ELEMENTS(haplotypeBasesArrayVector[j].first, haplotypeBasesArrayVector[j].second, JNI_RO_RELEASE_MODE);
    env->DeleteGlobalRef(haplotypeBasesArrayVector[j].first);	//free the global reference
  }
  haplotypeBasesArrayVector.clear();
  g_haplotypeBasesLengths.clear(); 
}


JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniClose
  (JNIEnv* env, jobject thisObject)
{
#ifdef DO_PROFILING
  g_load_time_initializer.print_profiling();
#endif
#ifdef DEBUG
  g_load_time_initializer.debug_close();
#endif
}

