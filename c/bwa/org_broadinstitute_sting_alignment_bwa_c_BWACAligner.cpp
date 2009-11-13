#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "bwa_gateway.h"
#include "org_broadinstitute_sting_alignment_bwa_c_BWACAligner.h"

static jclass java_alignment_array_class = NULL;
static jclass java_alignment_class = NULL;
static jmethodID java_alignment_constructor = NULL;

typedef void (BWA::*int_setter)(int value);
typedef void (BWA::*float_setter)(float value);

static jstring get_configuration_string(JNIEnv* env, jobject configuration, const char* field_name);

static void set_int_configuration_param(JNIEnv* env, jobject configuration, const char* field_name, BWA* bwa, int_setter setter);
static void set_float_configuration_param(JNIEnv* env, jobject configuration, const char* field_name, BWA* bwa, float_setter setter);

JNIEXPORT jlong JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_create(JNIEnv* env, jobject instance, jobject configuration)
{
  jstring java_ann = get_configuration_string(env,configuration,"annFileName");
  jstring java_amb = get_configuration_string(env,configuration,"ambFileName");
  jstring java_pac = get_configuration_string(env,configuration,"pacFileName");
  jstring java_forward_bwt = get_configuration_string(env,configuration,"forwardBWTFileName");
  jstring java_forward_sa = get_configuration_string(env,configuration,"forwardSAFileName");
  jstring java_reverse_bwt = get_configuration_string(env,configuration,"reverseBWTFileName");
  jstring java_reverse_sa = get_configuration_string(env,configuration,"reverseSAFileName");

  const char* ann_filename = env->GetStringUTFChars(java_ann,JNI_FALSE);
  const char* amb_filename = env->GetStringUTFChars(java_amb,JNI_FALSE);
  const char* pac_filename = env->GetStringUTFChars(java_pac,JNI_FALSE);
  const char* forward_bwt_filename = env->GetStringUTFChars(java_forward_bwt,JNI_FALSE);
  const char* forward_sa_filename = env->GetStringUTFChars(java_forward_sa,JNI_FALSE);
  const char* reverse_bwt_filename = env->GetStringUTFChars(java_reverse_bwt,JNI_FALSE);
  const char* reverse_sa_filename = env->GetStringUTFChars(java_reverse_sa,JNI_FALSE);

  BWA* bwa = new BWA(ann_filename,
		     amb_filename,
		     pac_filename,
		     forward_bwt_filename,
		     forward_sa_filename,
		     reverse_bwt_filename,
		     reverse_sa_filename);

  set_float_configuration_param(env, configuration, "maximumEditDistance", bwa, &BWA::set_max_edit_distance);
  set_int_configuration_param(env, configuration, "maximumGapOpens", bwa, &BWA::set_max_gap_opens);
  set_int_configuration_param(env, configuration, "maximumGapExtensions", bwa, &BWA::set_max_gap_extensions);
  set_int_configuration_param(env, configuration, "disallowIndelWithinRange", bwa, &BWA::set_disallow_indel_within_range);
  set_int_configuration_param(env, configuration, "mismatchPenalty", bwa, &BWA::set_mismatch_penalty);
  set_int_configuration_param(env, configuration, "gapOpenPenalty", bwa, &BWA::set_gap_open_penalty);
  set_int_configuration_param(env, configuration, "gapExtensionPenalty", bwa, &BWA::set_gap_extension_penalty);

  env->ReleaseStringUTFChars(java_ann,ann_filename);
  env->ReleaseStringUTFChars(java_amb,amb_filename);
  env->ReleaseStringUTFChars(java_pac,pac_filename);
  env->ReleaseStringUTFChars(java_forward_bwt,forward_bwt_filename);
  env->ReleaseStringUTFChars(java_forward_sa,forward_sa_filename);
  env->ReleaseStringUTFChars(java_reverse_bwt,reverse_bwt_filename);
  env->ReleaseStringUTFChars(java_reverse_sa,reverse_sa_filename);

  // Cache the class object for an array of alignments.
  java_alignment_array_class = env->FindClass("[Lorg/broadinstitute/sting/alignment/Alignment;");
  java_alignment_class = env->FindClass("org/broadinstitute/sting/alignment/Alignment");
  java_alignment_constructor = env->GetMethodID(java_alignment_class, "<init>", "(IIZI[C[I)V");

  return (jlong)bwa;
}

JNIEXPORT void JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_destroy(JNIEnv* env, jobject instance, jlong java_bwa) 
{
  BWA* bwa = (BWA*)java_bwa;
  delete bwa;
}

JNIEXPORT jobjectArray JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_getAlignments(JNIEnv* env, jobject object, jlong java_bwa, jbyteArray java_bases) {
  BWA* bwa = (BWA*)java_bwa;

  const jsize read_length = env->GetArrayLength(java_bases);
  jbyte *read_bases = env->GetByteArrayElements(java_bases,JNI_FALSE); 

  Alignment* alignments = NULL;
  unsigned num_alignments = 0;
  bwa->align((const char*)read_bases,read_length,alignments,num_alignments);

  jobjectArray java_alignments = env->NewObjectArray(num_alignments, env->FindClass("org/broadinstitute/sting/alignment/Alignment"), NULL);  

  for(unsigned alignment_idx = 0; alignment_idx < (unsigned)num_alignments; alignment_idx++) {
    Alignment& alignment = *(alignments + alignment_idx);

    unsigned cigar_length;
    if(alignment.type == BWA_TYPE_NO_MATCH) cigar_length = 0;
    else if(!alignment.cigar) cigar_length = 1;
    else cigar_length = alignment.n_cigar;

    jcharArray java_cigar_operators = env->NewCharArray(cigar_length);
    jintArray java_cigar_lengths = env->NewIntArray(cigar_length);

    if(alignment.cigar) {
      for(unsigned cigar_idx = 0; cigar_idx < (unsigned)alignment.n_cigar; ++cigar_idx) {
	jchar cigar_operator = "MIDS"[alignment.cigar[cigar_idx]>>14];
	jint  cigar_length = alignment.cigar[cigar_idx]&0x3fff;
	env->SetCharArrayRegion(java_cigar_operators,cigar_idx,1,&cigar_operator);
	env->SetIntArrayRegion(java_cigar_lengths,cigar_idx,1,&cigar_length);
      }
    }
    else {
      if(alignment.type != BWA_TYPE_NO_MATCH) {
	jchar cigar_operator = 'M';
	env->SetCharArrayRegion(java_cigar_operators,0,1,&cigar_operator);
	env->SetIntArrayRegion(java_cigar_lengths,0,1,&read_length);
      }
    }

    jclass java_alignment_class = env->FindClass("org/broadinstitute/sting/alignment/Alignment");
    jmethodID java_alignment_constructor = env->GetMethodID(java_alignment_class, "<init>", "(IIZI[C[I)V");
    jobject java_alignment = env->NewObject(java_alignment_class,
					    java_alignment_constructor,
					    alignment.contig,
					    alignment.pos,
					    alignment.negative_strand,
					    alignment.mapQ,
					    java_cigar_operators,
					    java_cigar_lengths);
    env->SetObjectArrayElement(java_alignments,alignment_idx,java_alignment);

    delete[] alignment.cigar;
  }

  delete[] alignments;

  env->ReleaseByteArrayElements(java_bases,read_bases,0);

  return java_alignments;
}

static jstring get_configuration_string(JNIEnv* env, jobject configuration, const char* field_name) {
  jclass configuration_class = env->GetObjectClass(configuration);
  jfieldID configuration_field = env->GetFieldID(configuration_class, field_name, "Ljava/lang/String;");
  return (jstring)env->GetObjectField(configuration,configuration_field);
}

static void set_int_configuration_param(JNIEnv* env, jobject configuration, const char* field_name, BWA* bwa, int_setter setter) {
  jclass configuration_class = env->GetObjectClass(configuration);
  jfieldID configuration_field = env->GetFieldID(configuration_class, field_name, "Ljava/lang/Integer;");
  jobject boxed_value = env->GetObjectField(configuration,configuration_field);
  if(boxed_value != NULL) {
    jclass int_box_class = env->FindClass("java/lang/Integer");
    jmethodID int_extractor = env->GetMethodID(int_box_class,"intValue", "()I");
    jint value = env->CallIntMethod(boxed_value,int_extractor);
    (bwa->*setter)(value);
  }
}

static void set_float_configuration_param(JNIEnv* env, jobject configuration, const char* field_name, BWA* bwa, float_setter setter) {
  jclass configuration_class = env->GetObjectClass(configuration);
  jfieldID configuration_field = env->GetFieldID(configuration_class, field_name, "Ljava/lang/Float;");
  jobject boxed_value = env->GetObjectField(configuration,configuration_field);
  if(boxed_value != NULL) {
    jclass float_box_class = env->FindClass("java/lang/Float");
    jmethodID float_extractor = env->GetMethodID(float_box_class,"floatValue", "()F");
    jfloat value = env->CallFloatMethod(boxed_value,float_extractor);
    (bwa->*setter)(value);
  }
}
