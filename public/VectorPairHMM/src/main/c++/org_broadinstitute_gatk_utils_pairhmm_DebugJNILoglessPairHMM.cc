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
#include "org_broadinstitute_gatk_utils_pairhmm_DebugJNILoglessPairHMM.h"
#include "utils.h"
#include "LoadTimeInitializer.h"
#include "jnidebug.h"
DataHolder<double> g_double_dataholder;

using namespace std;

JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_DebugJNILoglessPairHMM_jniInitialize
(JNIEnv* env, jobject thisObject,
 jint readMaxLength, jint haplotypeMaxLength)
{
  static int g_num_init_calls = 0;
#ifdef DEBUG3
  cout << "Entered alloc initialized .. readMaxLength "<<readMaxLength<<" haplotypeMaxLength "<<haplotypeMaxLength<<"\n";
#endif
  g_double_dataholder.initialize(readMaxLength, haplotypeMaxLength);
#ifdef DEBUG3
  debug_dump("lengths_jni.txt", to_string(readMaxLength)+" "+to_string(haplotypeMaxLength),true);
#endif
  ++g_num_init_calls;
}

JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_DebugJNILoglessPairHMM_jniInitializeProbabilities
(JNIEnv* env, jclass thisObject,
 jobjectArray transition, jbyteArray insertionGOP, jbyteArray deletionGOP, jbyteArray overallGCP
 )
{
  jboolean is_copy = JNI_FALSE;
  jsize length = (env)->GetArrayLength(insertionGOP);
#ifdef DEBUG3
  cout << "Entered initializeProbabilities .. length "<<length<<"\n";
#endif
  jbyte* insertionGOPArray = (env)->GetByteArrayElements(insertionGOP, &is_copy);
  jbyte* deletionGOPArray = (env)->GetByteArrayElements(deletionGOP, &is_copy);
  jbyte* overallGCPArray = (env)->GetByteArrayElements(overallGCP, &is_copy);
#ifdef DEBUG
  if(insertionGOPArray == 0)
    cerr << "insertionGOP array not initialized in JNI\n";
  ////assert(insertionGOPArray && "insertionGOP array not initialized in JNI");
  if(deletionGOPArray == 0)
    cerr << "deletionGOP array not initialized in JNI\n";
  ////assert(deletionGOPArray && "deletionGOP array not initialized in JNI");
  assert(overallGCPArray && "OverallGCP array not initialized in JNI");
#endif

  g_double_dataholder.initializeProbabilities(length, insertionGOPArray, deletionGOPArray, overallGCPArray);

  env->ReleaseByteArrayElements(overallGCP, overallGCPArray, JNI_ABORT);
  env->ReleaseByteArrayElements(deletionGOP, deletionGOPArray, JNI_ABORT);
  env->ReleaseByteArrayElements(insertionGOP, insertionGOPArray, JNI_ABORT);
}

JNIEXPORT jdouble JNICALL 
Java_org_broadinstitute_gatk_utils_pairhmm_DebugJNILoglessPairHMM_jniInitializePriorsAndUpdateCells( 
    JNIEnv* env, jobject thisObject,
    jboolean doInitialization, jint paddedReadLength, jint paddedHaplotypeLength,
    jbyteArray readBases, jbyteArray haplotypeBases, jbyteArray readQuals,
    jint hapStartIndex
    )
{
#ifdef DEBUG3
  cout << "Entered mainCompute .. doInitialization "<<(doInitialization == JNI_TRUE)<<" hapStartIndex "<<hapStartIndex<<"\n";
  cout << "mainCompute padded lengths "<< paddedReadLength << " " << paddedHaplotypeLength <<"\n";
#endif
  jboolean is_copy = JNI_FALSE;
  jbyte* readBasesArray = (env)->GetByteArrayElements(readBases, &is_copy);
  jbyte* haplotypeBasesArray = (env)->GetByteArrayElements(haplotypeBases, &is_copy);
  jbyte* readQualsArray = (env)->GetByteArrayElements(readQuals, &is_copy);
#ifdef DEBUG
  assert(readBasesArray && "readBasesArray not initialized in JNI"); 
  assert(haplotypeBasesArray && "haplotypeBasesArray not initialized in JNI"); 
  assert(readQualsArray && "readQualsArray not initialized in JNI"); 
#endif
  testcase tc;

  tc.rslen = paddedReadLength-1;
  tc.haplen = paddedHaplotypeLength-1;

  tc.rs = (char*)readBasesArray;
  tc.hap = (char*)haplotypeBasesArray;
  tc.q = (char*)readQualsArray; //TOASK - q is now char*

  compute_full_prob<double>(&tc, g_double_dataholder.m_matchMatrix, g_double_dataholder.m_insertionMatrix,
      g_double_dataholder.m_deletionMatrix, g_double_dataholder.m_transition, 
    doInitialization == JNI_TRUE, hapStartIndex, NULL);

  env->ReleaseByteArrayElements(readBases, readBasesArray, JNI_ABORT);
  env->ReleaseByteArrayElements(haplotypeBases, haplotypeBasesArray, JNI_ABORT);
  env->ReleaseByteArrayElements(readQuals, readQualsArray, JNI_ABORT);
  return 0.0;
}

JNIEXPORT jdouble JNICALL 
Java_org_broadinstitute_gatk_utils_pairhmm_DebugJNILoglessPairHMM_jniSubComputeReadLikelihoodGivenHaplotypeLog10( 
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
  //assert(readLength < MROWS);
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

  double result_avxd = g_compute_full_prob_double(&tc, 0);
  double result = log10(result_avxd) - log10(ldexp(1.0, 1020));
#ifdef DEBUG
  g_load_time_initializer.debug_dump("return_values_jni.txt",to_string(result),true);
#endif


  RELEASE_BYTE_ARRAY_ELEMENTS(overallGCP, overallGCPArray, JNI_RO_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(deletionGOP, deletionGOPArray, JNI_RO_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(insertionGOP, insertionGOPArray, JNI_RO_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(readQuals, readQualsArray, JNI_RO_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(haplotypeBases, haplotypeBasesArray, JNI_RO_RELEASE_MODE);
  RELEASE_BYTE_ARRAY_ELEMENTS(readBases, readBasesArray, JNI_RO_RELEASE_MODE);

  return 0.0;
}

