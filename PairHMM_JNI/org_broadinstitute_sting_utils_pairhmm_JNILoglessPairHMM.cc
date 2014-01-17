#include <jni.h>
#include <assert.h>
#include <stdio.h>
#include "org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM.h"
#include <vector>

#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>


#include <immintrin.h>
#include <emmintrin.h>
#include <omp.h>
using namespace std;

#define ENABLE_ASSERTIONS 1
//#define DEBUG 1
//#define DEBUG0_1 1
//#define DEBUG3 1


#include "template.h"

#include "define-double.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"



#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5

class LoadTimeInitializer
{
  public:
    LoadTimeInitializer()		//will be called when library is loaded
    {
      ConvertChar::init();
    }
    jfieldID m_readBasesFID;
    jfieldID m_readQualsFID;
    jfieldID m_insertionGOPFID;
    jfieldID m_deletionGOPFID;
    jfieldID m_overallGCPFID;
    jfieldID m_haplotypeBasesFID;
};
LoadTimeInitializer g_load_time_initializer;

void debug_dump(string filename, string s, bool to_append, bool add_newline)
{
  ofstream fptr;
  fptr.open(filename.c_str(), to_append ? ofstream::app : ofstream::out);
  fptr << s;
  if(add_newline)
    fptr << "\n";
  fptr.close();
}

JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM_jniInitializeProbabilities
(JNIEnv* env, jclass thisObject,
 jobjectArray transition, jbyteArray insertionGOP, jbyteArray deletionGOP, jbyteArray overallGCP
 )
{}

  JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM_jniInitialize
(JNIEnv* env, jobject thisObject,
 jint readMaxLength, jint haplotypeMaxLength)
{}

JNIEXPORT jdouble JNICALL 
Java_org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM_jniInitializePriorsAndUpdateCells( 
    JNIEnv* env, jobject thisObject,
    jboolean doInitialization, jint paddedReadLength, jint paddedHaplotypeLength,
    jbyteArray readBases, jbyteArray haplotypeBases, jbyteArray readQuals,
    jint hapStartIndex
    )

{ return 0.0; }

#define DIRECT_ACCESS_TO_JAVA_HEAP_MEMORY 1

#ifdef DIRECT_ACCESS_TO_JAVA_HEAP_MEMORY
//Gets direct access to Java arrays
#define GET_BYTE_ARRAY_ELEMENTS env->GetPrimitiveArrayCritical
#define RELEASE_BYTE_ARRAY_ELEMENTS env->ReleasePrimitiveArrayCritical
#define JNI_RELEASE_MODE JNI_ABORT
#define GET_DOUBLE_ARRAY_ELEMENTS env->GetPrimitiveArrayCritical
#define RELEASE_DOUBLE_ARRAY_ELEMENTS env->ReleasePrimitiveArrayCritical

#else
//Likely makes copy of Java arrays to JNI C++ space
#define GET_BYTE_ARRAY_ELEMENTS env->GetByteArrayElements
#define RELEASE_BYTE_ARRAY_ELEMENTS env->ReleaseByteArrayElements
#define JNI_RELEASE_MODE JNI_ABORT
#define GET_DOUBLE_ARRAY_ELEMENTS env->GetDoubleArrayElements
#define RELEASE_DOUBLE_ARRAY_ELEMENTS env->ReleaseDoubleArrayElements

#endif		//ifdef DIRECT_ACCESS_TO_JAVA_HEAP_MEMORY

JNIEXPORT jdouble JNICALL 
Java_org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM_jniSubComputeReadLikelihoodGivenHaplotypeLog10( 
    JNIEnv* env, jobject thisObject,
    jint readLength, jint haplotypeLength,
    jbyteArray readBases, jbyteArray haplotypeBases, jbyteArray readQuals,
    jbyteArray insertionGOP, jbyteArray deletionGOP, jbyteArray overallGCP,
    jint hapStartIndex
    )
{
  jboolean is_copy = JNI_FALSE;
  jbyte* readBasesArray =   (jbyte*)GET_BYTE_ARRAY_ELEMENTS(readBases, &is_copy);
  jbyte* haplotypeBasesArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(haplotypeBases, &is_copy);
  jbyte* readQualsArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(readQuals, &is_copy);
  jbyte* insertionGOPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(insertionGOP, &is_copy);
  jbyte* deletionGOPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(deletionGOP, &is_copy);
  jbyte* overallGCPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(overallGCP, &is_copy);
#ifdef DEBUG
  assert(readBasesArray && "readBasesArray not initialized in JNI"); 
  assert(haplotypeBasesArray && "haplotypeBasesArray not initialized in JNI"); 
  assert(readQualsArray && "readQualsArray not initialized in JNI"); 
  assert(insertionGOPArray && "insertionGOP array not initialized in JNI");
  assert(deletionGOPArray && "deletionGOP array not initialized in JNI");
  assert(overallGCPArray && "OverallGCP array not initialized in JNI");
  assert(readLength < MROWS);
#endif
  testcase tc;
  tc.rslen = readLength;
  tc.haplen = haplotypeLength;
  tc.rs = (char*)readBasesArray;
  tc.hap = (char*)haplotypeBasesArray;
  for(unsigned i=0;i<readLength;++i)
  {
    tc.q[i] = (int)readQualsArray[i];
    tc.i[i] = (int)insertionGOPArray[i];
    tc.d[i] = (int)deletionGOPArray[i];
    tc.c[i] = (int)overallGCPArray[i];
  }

  double result_avxd = GEN_INTRINSIC(compute_full_prob_avx, d)<double>(&tc);
  double result = log10(result_avxd) - log10(ldexp(1.0, 1020));
#ifdef DEBUG
  debug_dump("return_values_jni.txt",to_string(result),true);
#endif


  RELEASE_BYTE_ARRAY_ELEMENTS(overallGCP, overallGCPArray, JNI_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(deletionGOP, deletionGOPArray, JNI_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(insertionGOP, insertionGOPArray, JNI_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(readQuals, readQualsArray, JNI_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(haplotypeBases, haplotypeBasesArray, JNI_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(readBases, readBasesArray, JNI_RELEASE_MODE);

  return 0.0;
}

//Should be called only once for the whole Java process - initializes field ids for the classes JNIReadDataHolderClass
//and JNIHaplotypeDataHolderClass
JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM_jniGlobalInit
  (JNIEnv* env, jobject thisObject, jclass readDataHolderClass, jclass haplotypeDataHolderClass)
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

//JNI function to invoke compute_full_prob_avx
//readDataArray - array of JNIReadDataHolderClass objects which contain the readBases, readQuals etc
//haplotypeDataArray - array of JNIHaplotypeDataHolderClass objects which contain the haplotypeBases
//likelihoodArray - array of doubles to return results back to Java. Memory allocated by Java prior to JNI call
//maxNumThreadsToUse - Max number of threads that OpenMP can use for the HMM computation
JNIEXPORT void JNICALL Java_org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM_jniComputeLikelihoods
  (JNIEnv* env, jobject thisObject, jint numReads, jint numHaplotypes, 
   jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray, jint maxNumThreadsToUse)
{
  jboolean is_copy = JNI_FALSE;
  //To ensure, GET_BYTE_ARRAY_ELEMENTS is invoked only once for each haplotype, store bytearrays in a vector
  vector<pair<jbyteArray, jbyte*> > haplotypeBasesArrayVector;
  haplotypeBasesArrayVector.clear();
  haplotypeBasesArrayVector.resize(numHaplotypes);
#ifdef DEBUG0_1
  cout << "JNI numReads "<<numReads<<" numHaplotypes "<<numHaplotypes<<"\n";
#endif
  for(unsigned j=0;j<numHaplotypes;++j)
  {
    jobject haplotypeObject = env->GetObjectArrayElement(haplotypeDataArray, j);
    jbyteArray haplotypeBases = (jbyteArray)env->GetObjectField(haplotypeObject, g_load_time_initializer.m_haplotypeBasesFID);
#ifdef ENABLE_ASSERTIONS
    assert(haplotypeBases && ("haplotypeBases is NULL at index : "+to_string(j)+"\n").c_str());
#endif
    jbyte* haplotypeBasesArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(haplotypeBases, &is_copy);
#ifdef ENABLE_ASSERTIONS
    assert(haplotypeBasesArray && "haplotypeBasesArray not initialized in JNI"); 
    assert(env->GetArrayLength(haplotypeBases) < MCOLS);
#endif
#ifdef DEBUG0_1
    cout << "JNI haplotype length "<<env->GetArrayLength(haplotypeBases)<<"\n";
#endif
    haplotypeBasesArrayVector[j] = make_pair(haplotypeBases, haplotypeBasesArray);
#ifdef DEBUG3
    for(unsigned k=0;k<env->GetArrayLength(haplotypeBases);++k)
      debug_dump("haplotype_bases_jni.txt",to_string((int)haplotypeBasesArray[k]),true);
#endif
  }

  unsigned numTestCases = numReads*numHaplotypes;
  vector<testcase> tc_array;
  tc_array.clear();
  tc_array.resize(numTestCases);
  unsigned tc_idx = 0;
  //Store arrays for release later
  vector<vector<pair<jbyteArray,jbyte*> > > readBasesArrayVector;
  readBasesArrayVector.clear();
  readBasesArrayVector.resize(numReads);
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
      debug_dump("reads_jni.txt",to_string((int)readBasesArray[j]),true);
      debug_dump("reads_jni.txt",to_string((int)readQualsArray[j]),true);
      debug_dump("reads_jni.txt",to_string((int)insertionGOPArray[j]),true);
      debug_dump("reads_jni.txt",to_string((int)deletionGOPArray[j]),true);
      debug_dump("reads_jni.txt",to_string((int)overallGCPArray[j]),true);
    }
#endif

    for(unsigned j=0;j<numHaplotypes;++j)
    {
      jsize haplotypeLength = env->GetArrayLength(haplotypeBasesArrayVector[j].first);
      jbyte* haplotypeBasesArray = haplotypeBasesArrayVector[j].second;
      tc_array[tc_idx].rslen = (int)readLength;
      tc_array[tc_idx].haplen = (int)haplotypeLength;
      tc_array[tc_idx].rs = (char*)readBasesArray;
      tc_array[tc_idx].hap = (char*)haplotypeBasesArray;
      //Can be avoided 
      for(unsigned k=0;k<readLength;++k)
      {
	tc_array[tc_idx].q[k] = (int)readQualsArray[k];
	tc_array[tc_idx].i[k] = (int)insertionGOPArray[k];
	tc_array[tc_idx].d[k] = (int)deletionGOPArray[k];
	tc_array[tc_idx].c[k] = (int)overallGCPArray[k];
      }

      ++tc_idx;  
    }
    RELEASE_BYTE_ARRAY_ELEMENTS(overallGCP, overallGCPArray, JNI_RELEASE_MODE);	//order of GET-RELEASE is important
    RELEASE_BYTE_ARRAY_ELEMENTS(deletionGOP, deletionGOPArray, JNI_RELEASE_MODE);
    RELEASE_BYTE_ARRAY_ELEMENTS(insertionGOP, insertionGOPArray, JNI_RELEASE_MODE);
    RELEASE_BYTE_ARRAY_ELEMENTS(readQuals, readQualsArray, JNI_RELEASE_MODE);

    //Release readBases at end because it is used by compute_full_prob
    readBasesArrayVector[i].clear();
    readBasesArrayVector[i].resize(1);
    readBasesArrayVector[i][0] = make_pair(readBases, readBasesArray);
  }

  jdouble* likelihoodDoubleArray = (jdouble*)GET_DOUBLE_ARRAY_ELEMENTS(likelihoodArray, &is_copy);
#ifdef ENABLE_ASSERTIONS
  assert(likelihoodDoubleArray && "likelihoodArray is NULL");
  assert(env->GetArrayLength(likelihoodArray) == numTestCases);
#endif
#pragma omp parallel for schedule (dynamic,10) private(tc_idx) num_threads(maxNumThreadsToUse) 
  for(tc_idx=0;tc_idx<numTestCases;++tc_idx)
  {
    double result_avxd = GEN_INTRINSIC(compute_full_prob_avx, d)<double>(&(tc_array[tc_idx]));
    double result = log10(result_avxd) - log10(ldexp(1.0, 1020));
    likelihoodDoubleArray[tc_idx] = result;
  }
#ifdef DEBUG
  for(tc_idx=0;tc_idx<numTestCases;++tc_idx)
  {
    debug_dump("return_values_jni.txt",to_string(likelihoodDoubleArray[tc_idx]),true);
  }
#endif
  RELEASE_DOUBLE_ARRAY_ELEMENTS(likelihoodArray, likelihoodDoubleArray, 0);

  //Release read arrays first
  for(int i=readBasesArrayVector.size()-1;i>=0;--i)//note the order - reverse of GET
  {
    for(int j=readBasesArrayVector.size()-1;j>=0;--j)
      RELEASE_BYTE_ARRAY_ELEMENTS(readBasesArrayVector[i][j].first, readBasesArrayVector[i][j].second, JNI_RELEASE_MODE);
    readBasesArrayVector[i].clear();
  }
  readBasesArrayVector.clear();

  //Now release haplotype arrays
  for(int j=numHaplotypes-1;j>=0;--j)	//note the order - reverse of GET
    RELEASE_BYTE_ARRAY_ELEMENTS(haplotypeBasesArrayVector[j].first, haplotypeBasesArrayVector[j].second, JNI_RELEASE_MODE);
  haplotypeBasesArrayVector.clear();
  tc_array.clear();
}

