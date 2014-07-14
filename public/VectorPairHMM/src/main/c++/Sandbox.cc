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


#include "Sandbox.h"
#include "org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM.h"
#include "utils.h"
#include "jni_common.h"
/*
 * Class:     Sandbox
 * Method:    jniGetMachineType
 * Signature: ()J
 */
JNIEXPORT jlong JNICALL Java_Sandbox_jniGetMachineType
  (JNIEnv * env, jobject thisObj)
{
  return 0;
}

/*
 * Class:     Sandbox
 * Method:    jniInitializeClassFieldsAndMachineMask
 * Signature: (Ljava/lang/Class;Ljava/lang/Class;J)V
 */
JNIEXPORT void JNICALL Java_Sandbox_jniInitializeClassFieldsAndMachineMask
  (JNIEnv* env, jobject thisObject, jclass readDataHolderClass, jclass haplotypeDataHolderClass, jlong mask)
{
  Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniInitializeClassFieldsAndMachineMask(env, thisObject, readDataHolderClass,
      haplotypeDataHolderClass, mask);
}

/*
 * Class:     Sandbox
 * Method:    jniInitializeHaplotypes
 * Signature: (I[LSandbox/JNIHaplotypeDataHolderClass;)V
 */
JNIEXPORT void JNICALL Java_Sandbox_jniInitializeHaplotypes
  (JNIEnv * env, jobject thisObject, jint numHaplotypes, jobjectArray haplotypeDataArray)
{
  Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniInitializeHaplotypes(env, thisObject, numHaplotypes, haplotypeDataArray);
}

/*
 * Class:     Sandbox
 * Method:    jniFinalizeRegion
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_Sandbox_jniFinalizeRegion
  (JNIEnv * env, jobject thisObject)
{
  Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniFinalizeRegion(env, thisObject);
}


/*
 * Class:     Sandbox
 * Method:    jniComputeLikelihoods
 * Signature: (II[LSandbox/JNIReadDataHolderClass;[LSandbox/JNIHaplotypeDataHolderClass;[DI)V
 */
JNIEXPORT void JNICALL Java_Sandbox_jniComputeLikelihoods
  (JNIEnv* env, jobject thisObject, jint numReads, jint numHaplotypes, 
   jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray, jint maxNumThreadsToUse)
{
  Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniComputeLikelihoods(env, thisObject,
      numReads, numHaplotypes, readDataArray, haplotypeDataArray, likelihoodArray, maxNumThreadsToUse);
}
/*
 * Class:     Sandbox
 * Method:    jniClose
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_Sandbox_jniClose
  (JNIEnv* env, jobject thisObject)
{ Java_org_broadinstitute_gatk_utils_pairhmm_VectorLoglessPairHMM_jniClose(env, thisObject); } 

JNIEXPORT void JNICALL Java_Sandbox_doEverythingNative
  (JNIEnv* env, jobject thisObject, jstring fileNameString)
{
  const char* fileName = env->GetStringUTFChars(fileNameString, 0);
  char local_array[800];
  strncpy(local_array, fileName, 200);
  env->ReleaseStringUTFChars(fileNameString, fileName);
  do_compute(local_array, true, 10000, false);
}

