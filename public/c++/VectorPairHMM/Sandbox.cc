#include "Sandbox.h"
#include "org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM.h"
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
  Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniInitializeClassFieldsAndMachineMask(env, thisObject, readDataHolderClass,
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
  Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniInitializeHaplotypes(env, thisObject, numHaplotypes, haplotypeDataArray);
}

/*
 * Class:     Sandbox
 * Method:    jniFinalizeRegion
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_Sandbox_jniFinalizeRegion
  (JNIEnv * env, jobject thisObject)
{
  Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniFinalizeRegion(env, thisObject);
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
  Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniComputeLikelihoods(env, thisObject,
      numReads, numHaplotypes, readDataArray, haplotypeDataArray, likelihoodArray, maxNumThreadsToUse);
}
/*
 * Class:     Sandbox
 * Method:    jniClose
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_Sandbox_jniClose
  (JNIEnv* env, jobject thisObject)
{ Java_org_broadinstitute_sting_utils_pairhmm_VectorLoglessPairHMM_jniClose(env, thisObject); } 

JNIEXPORT void JNICALL Java_Sandbox_doEverythingNative
  (JNIEnv* env, jobject thisObject, jstring fileNameString)
{
  const char* fileName = env->GetStringUTFChars(fileNameString, 0);
  do_compute((char*)fileName);
  env->ReleaseStringUTFChars(fileNameString, fileName);
}

