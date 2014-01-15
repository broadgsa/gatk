#include <jni.h>
#include <assert.h>
#include <stdio.h>
#include "org_broadinstitute_sting_utils_pairhmm_JNILoglessPairHMM.h"

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

#include "template.h"

#include "define-double.h"
#include "shift_template.c"
#include "pairhmm-template-kernel.cc"

using namespace std;


#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5

//#define DEBUG3 1
#define DEBUG 1

template<class T>
string to_string(T obj)
{
  stringstream ss;
  string ret_string;
  ss.clear();
  ss << std::scientific << obj;
  ss >> ret_string;
  ss.clear();
  return ret_string;
}

void debug_dump(string filename, string s, bool to_append, bool add_newline=true)
{
  ofstream fptr;
  fptr.open(filename.c_str(), to_append ? ofstream::app : ofstream::out);
  fptr << s;
  if(add_newline)
    fptr << "\n";
  fptr.close();
}

#define INT_STORE_ARRAY_SIZE 2048
int insertionGOPIntArray[INT_STORE_ARRAY_SIZE];
int deletionGOPIntArray[INT_STORE_ARRAY_SIZE];
int overallGCPIntArray[INT_STORE_ARRAY_SIZE];
int readQualsIntArray[INT_STORE_ARRAY_SIZE];

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
  jbyte* readBasesArray = (env)->GetByteArrayElements(readBases, &is_copy);
  jbyte* haplotypeBasesArray = (env)->GetByteArrayElements(haplotypeBases, &is_copy);
  jbyte* readQualsArray = (env)->GetByteArrayElements(readQuals, &is_copy);
  jbyte* insertionGOPArray = (env)->GetByteArrayElements(insertionGOP, &is_copy);
  jbyte* deletionGOPArray = (env)->GetByteArrayElements(deletionGOP, &is_copy);
  jbyte* overallGCPArray = (env)->GetByteArrayElements(overallGCP, &is_copy);
#ifdef DEBUG
  assert(readBasesArray && "readBasesArray not initialized in JNI"); 
  assert(haplotypeBasesArray && "haplotypeBasesArray not initialized in JNI"); 
  assert(readQualsArray && "readQualsArray not initialized in JNI"); 
  assert(insertionGOPArray && "insertionGOP array not initialized in JNI");
  assert(deletionGOPArray && "deletionGOP array not initialized in JNI");
  assert(overallGCPArray && "OverallGCP array not initialized in JNI");
  assert(readLength < INT_STORE_ARRAY_SIZE);
#endif
  for(unsigned i=0;i<readLength;++i)
  {
    insertionGOPIntArray[i] = (int)insertionGOPArray[i];
    deletionGOPIntArray[i] = (int)deletionGOPArray[i];
    overallGCPIntArray[i] = (int)overallGCPArray[i];
    readQualsIntArray[i] = (int)readQualsArray[i];
  }

  testcase tc;
  tc.rslen = readLength;
  tc.haplen = haplotypeLength;

  tc.rs = (char*)readBasesArray;
  tc.hap = (char*)haplotypeBasesArray;
  tc.q = (int*)readQualsIntArray;
  tc.i = (int*)insertionGOPIntArray;
  tc.d = (int*)deletionGOPIntArray;
  tc.c = (int*)overallGCPIntArray;

  double result_avxd = GEN_INTRINSIC(compute_full_prob_avx, d)<double>(&tc);
  double result = log10(result_avxd) - log10(ldexp(1.0, 1020));
#ifdef DEBUG
  debug_dump("return_values_jni.txt",to_string(result),true);
#endif


  env->ReleaseByteArrayElements(overallGCP, overallGCPArray, JNI_ABORT);
  env->ReleaseByteArrayElements(deletionGOP, deletionGOPArray, JNI_ABORT);
  env->ReleaseByteArrayElements(insertionGOP, insertionGOPArray, JNI_ABORT);
  env->ReleaseByteArrayElements(readQuals, readQualsArray, JNI_ABORT);
  env->ReleaseByteArrayElements(haplotypeBases, haplotypeBasesArray, JNI_ABORT);
  env->ReleaseByteArrayElements(readBases, readBasesArray, JNI_ABORT);

  return 0.0;
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
