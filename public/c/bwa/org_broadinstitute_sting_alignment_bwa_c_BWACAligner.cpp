#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "bwa_gateway.h"
#include "org_broadinstitute_sting_alignment_bwa_c_BWACAligner.h"

typedef void (BWA::*int_setter)(int value);
typedef void (BWA::*float_setter)(float value);

static jobject convert_to_java_alignment(JNIEnv* env, const jbyte* read_bases, const jsize read_length, const Alignment& alignment);
static jstring get_configuration_file(JNIEnv* env, jobject configuration, const char* field_name);
static void set_int_configuration_param(JNIEnv* env, jobject configuration, const char* field_name, BWA* bwa, int_setter setter);
static void set_float_configuration_param(JNIEnv* env, jobject configuration, const char* field_name, BWA* bwa, float_setter setter);
static void throw_config_value_exception(JNIEnv* env, const char* field_name, const char* message);

JNIEXPORT jlong JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_create(JNIEnv* env, jobject instance, jobject bwtFiles, jobject configuration)
{
  jstring java_ann = get_configuration_file(env,bwtFiles,"annFile");
  if(java_ann == NULL) return 0L;
  jstring java_amb = get_configuration_file(env,bwtFiles,"ambFile");
  if(java_amb == NULL) return 0L;
  jstring java_pac = get_configuration_file(env,bwtFiles,"pacFile");
  if(java_pac == NULL) return 0L;
  jstring java_forward_bwt = get_configuration_file(env,bwtFiles,"forwardBWTFile");
  if(java_forward_bwt == NULL) return 0L;
  jstring java_forward_sa = get_configuration_file(env,bwtFiles,"forwardSAFile");
  if(java_forward_sa == NULL) return 0L;
  jstring java_reverse_bwt = get_configuration_file(env,bwtFiles,"reverseBWTFile");
  if(java_reverse_bwt == NULL) return 0L;
  jstring java_reverse_sa = get_configuration_file(env,bwtFiles,"reverseSAFile");
  if(java_reverse_sa == NULL) return 0L;

  const char* ann_filename = env->GetStringUTFChars(java_ann,JNI_FALSE);
  if(env->ExceptionCheck()) return 0L;
  const char* amb_filename = env->GetStringUTFChars(java_amb,JNI_FALSE);
  if(env->ExceptionCheck()) return 0L;
  const char* pac_filename = env->GetStringUTFChars(java_pac,JNI_FALSE);
  if(env->ExceptionCheck()) return 0L;
  const char* forward_bwt_filename = env->GetStringUTFChars(java_forward_bwt,JNI_FALSE);
  if(env->ExceptionCheck()) return 0L;
  const char* forward_sa_filename = env->GetStringUTFChars(java_forward_sa,JNI_FALSE);
  if(env->ExceptionCheck()) return 0L; 
  const char* reverse_bwt_filename = env->GetStringUTFChars(java_reverse_bwt,JNI_FALSE);
  if(env->ExceptionCheck()) return 0L; 
  const char* reverse_sa_filename = env->GetStringUTFChars(java_reverse_sa,JNI_FALSE);
  if(env->ExceptionCheck()) return 0L; 

  BWA* bwa = new BWA(ann_filename,
		     amb_filename,
		     pac_filename,
		     forward_bwt_filename,
		     forward_sa_filename,
		     reverse_bwt_filename,
		     reverse_sa_filename);

  Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_updateConfiguration(env,instance,(jlong)bwa,configuration); 
  if(env->ExceptionCheck()) return 0L;

  env->ReleaseStringUTFChars(java_ann,ann_filename);
  if(env->ExceptionCheck()) return 0L; 
  env->ReleaseStringUTFChars(java_amb,amb_filename);
  if(env->ExceptionCheck()) return 0L; 
  env->ReleaseStringUTFChars(java_pac,pac_filename);
  if(env->ExceptionCheck()) return 0L; 
  env->ReleaseStringUTFChars(java_forward_bwt,forward_bwt_filename);
  if(env->ExceptionCheck()) return 0L; 
  env->ReleaseStringUTFChars(java_forward_sa,forward_sa_filename);
  if(env->ExceptionCheck()) return 0L; 
  env->ReleaseStringUTFChars(java_reverse_bwt,reverse_bwt_filename);
  if(env->ExceptionCheck()) return 0L; 
  env->ReleaseStringUTFChars(java_reverse_sa,reverse_sa_filename);
  if(env->ExceptionCheck()) return 0L; 

  return (jlong)bwa;
}

JNIEXPORT void JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_destroy(JNIEnv* env, jobject instance, jlong java_bwa) 
{
  BWA* bwa = (BWA*)java_bwa;
  delete bwa;
}

JNIEXPORT void JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_updateConfiguration(JNIEnv *env, jobject instance, jlong java_bwa, jobject configuration) {
  BWA* bwa = (BWA*)java_bwa;
  set_float_configuration_param(env, configuration, "maximumEditDistance", bwa, &BWA::set_max_edit_distance);
  if(env->ExceptionCheck()) return; 
  set_int_configuration_param(env, configuration, "maximumGapOpens", bwa, &BWA::set_max_gap_opens);
  if(env->ExceptionCheck()) return; 
  set_int_configuration_param(env, configuration, "maximumGapExtensions", bwa, &BWA::set_max_gap_extensions);
  if(env->ExceptionCheck()) return; 
  set_int_configuration_param(env, configuration, "disallowIndelWithinRange", bwa, &BWA::set_disallow_indel_within_range);
  if(env->ExceptionCheck()) return; 
  set_int_configuration_param(env, configuration, "mismatchPenalty", bwa, &BWA::set_mismatch_penalty);
  if(env->ExceptionCheck()) return; 
  set_int_configuration_param(env, configuration, "gapOpenPenalty", bwa, &BWA::set_gap_open_penalty);
  if(env->ExceptionCheck()) return; 
  set_int_configuration_param(env, configuration, "gapExtensionPenalty", bwa, &BWA::set_gap_extension_penalty);
  if(env->ExceptionCheck()) return;
}

JNIEXPORT jobjectArray JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_getPaths(JNIEnv *env, jobject instance, jlong java_bwa, jbyteArray java_bases) 
{
 BWA* bwa = (BWA*)java_bwa;

  const jsize read_length = env->GetArrayLength(java_bases);
  if(env->ExceptionCheck()) return NULL;

  jbyte *read_bases = env->GetByteArrayElements(java_bases,JNI_FALSE); 
  if(read_bases == NULL) return NULL;

  bwt_aln1_t* paths = NULL;
  unsigned num_paths = 0;

  unsigned best_path_count, second_best_path_count;
  bwa->find_paths((const char*)read_bases,read_length,paths,num_paths,best_path_count,second_best_path_count);

  jobjectArray java_paths = env->NewObjectArray(num_paths, env->FindClass("org/broadinstitute/sting/alignment/bwa/c/BWAPath"), NULL);  
  if(java_paths == NULL) return NULL;

  for(unsigned path_idx = 0; path_idx < (unsigned)num_paths; path_idx++) {
    bwt_aln1_t& path = *(paths + path_idx);

    jclass java_path_class = env->FindClass("org/broadinstitute/sting/alignment/bwa/c/BWAPath");
    if(java_path_class == NULL) return NULL;

    jmethodID java_path_constructor = env->GetMethodID(java_path_class, "<init>", "(IIIZJJIII)V");
    if(java_path_constructor == NULL) return NULL;

    // Note that k/l are being cast to long.  Bad things will happen if JNI assumes that they're ints.
    jobject java_path = env->NewObject(java_path_class,
                                       java_path_constructor,
                                       path.n_mm,
                                       path.n_gapo,
                                       path.n_gape,
                                       path.a,
                                       (jlong)path.k,
                                       (jlong)path.l,
                                       path.score,
                                       best_path_count,
                                       second_best_path_count);
    if(java_path == NULL) return NULL;      

    env->SetObjectArrayElement(java_paths,path_idx,java_path);
    if(env->ExceptionCheck()) return NULL;

    env->DeleteLocalRef(java_path_class);
    if(env->ExceptionCheck()) return NULL;
  }

  delete[] paths;

  env->ReleaseByteArrayElements(java_bases,read_bases,JNI_FALSE); 

  return env->ExceptionCheck() ? NULL : java_paths;
}

JNIEXPORT jobjectArray JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_convertPathsToAlignments(JNIEnv *env, jobject instance, jlong java_bwa, jbyteArray java_bases, jobjectArray java_paths) 
{
  BWA* bwa = (BWA*)java_bwa;

  const jsize read_length = env->GetArrayLength(java_bases);
  if(env->ExceptionCheck()) return NULL;

  jbyte *read_bases = env->GetByteArrayElements(java_bases,JNI_FALSE); 
  if(read_bases == NULL) return NULL;

  const jsize num_paths = env->GetArrayLength(java_paths);
  bwt_aln1_t* paths = new bwt_aln1_t[num_paths];
  unsigned best_count = 0, second_best_count = 0;

  for(unsigned path_idx = 0; path_idx < (unsigned)num_paths; path_idx++) {
    jobject java_path = env->GetObjectArrayElement(java_paths,path_idx);
    jclass java_path_class = env->GetObjectClass(java_path);
    if(java_path_class == NULL) return NULL;

    bwt_aln1_t& path = *(paths + path_idx);

    jfieldID mismatches_field = env->GetFieldID(java_path_class, "numMismatches", "I");
    if(mismatches_field == NULL) return NULL;
    path.n_mm = env->GetIntField(java_path,mismatches_field);
    if(env->ExceptionCheck()) return NULL;

    jfieldID gap_opens_field = env->GetFieldID(java_path_class, "numGapOpens", "I");
    if(gap_opens_field == NULL) return NULL;
    path.n_gapo = env->GetIntField(java_path,gap_opens_field);
    if(env->ExceptionCheck()) return NULL;

    jfieldID gap_extensions_field = env->GetFieldID(java_path_class, "numGapExtensions", "I");
    if(gap_extensions_field == NULL) return NULL;
    path.n_gape = env->GetIntField(java_path,gap_extensions_field);
    if(env->ExceptionCheck()) return NULL;

    jfieldID negative_strand_field = env->GetFieldID(java_path_class, "negativeStrand", "Z");
    if(negative_strand_field == NULL) return NULL;
    path.a = env->GetBooleanField(java_path,negative_strand_field);
    if(env->ExceptionCheck()) return NULL;

    jfieldID k_field = env->GetFieldID(java_path_class, "k", "J");
    if(k_field == NULL) return NULL;
    path.k = env->GetLongField(java_path,k_field);
    if(env->ExceptionCheck()) return NULL;

    jfieldID l_field = env->GetFieldID(java_path_class, "l", "J");
    if(l_field == NULL) return NULL;
    path.l = env->GetLongField(java_path,l_field);
    if(env->ExceptionCheck()) return NULL;

    jfieldID score_field = env->GetFieldID(java_path_class, "score", "I");
    if(score_field == NULL) return NULL;
    path.score = env->GetIntField(java_path,score_field);
    if(env->ExceptionCheck()) return NULL;

    jfieldID best_count_field = env->GetFieldID(java_path_class, "bestCount", "I");
    if(best_count_field == NULL) return NULL;
    best_count = env->GetIntField(java_path,best_count_field);
    if(env->ExceptionCheck()) return NULL;

    jfieldID second_best_count_field = env->GetFieldID(java_path_class, "secondBestCount", "I");
    if(second_best_count_field == NULL) return NULL;
    second_best_count = env->GetIntField(java_path,second_best_count_field);
    if(env->ExceptionCheck()) return NULL;
  }

  Alignment* alignments = NULL;
  unsigned num_alignments = 0;
  bwa->generate_alignments_from_paths((const char*)read_bases,read_length,paths,num_paths,best_count,second_best_count,alignments,num_alignments);

  jobjectArray java_alignments = env->NewObjectArray(num_alignments, env->FindClass("org/broadinstitute/sting/alignment/Alignment"), NULL);  
  if(java_alignments == NULL) return NULL;

  for(unsigned alignment_idx = 0; alignment_idx < (unsigned)num_alignments; alignment_idx++) {
    Alignment& alignment = *(alignments + alignment_idx);
    jobject java_alignment = convert_to_java_alignment(env,read_bases,read_length,alignment);
    if(java_alignment == NULL) return NULL;
    env->SetObjectArrayElement(java_alignments,alignment_idx,java_alignment);
    if(env->ExceptionCheck()) return NULL;
  }

  delete[] alignments;
  delete[] paths;

  env->ReleaseByteArrayElements(java_bases,read_bases,JNI_FALSE); 

  return env->ExceptionCheck() ? NULL : java_alignments;
}

JNIEXPORT jobject JNICALL Java_org_broadinstitute_sting_alignment_bwa_c_BWACAligner_getBestAlignment(JNIEnv *env, jobject instance, jlong java_bwa, jbyteArray java_bases) {
  BWA* bwa = (BWA*)java_bwa;

  const jsize read_length = env->GetArrayLength(java_bases);
  if(env->ExceptionCheck()) return NULL;

  jbyte *read_bases = env->GetByteArrayElements(java_bases,JNI_FALSE); 
  if(read_bases == NULL) return NULL;

  Alignment* best_alignment = bwa->generate_single_alignment((const char*)read_bases,read_length);
  jobject java_best_alignment = (best_alignment != NULL) ? convert_to_java_alignment(env,read_bases,read_length,*best_alignment) : NULL;
  delete best_alignment;

  env->ReleaseByteArrayElements(java_bases,read_bases,JNI_FALSE); 

  return java_best_alignment;
}

static jobject convert_to_java_alignment(JNIEnv *env, const jbyte* read_bases, const jsize read_length, const Alignment& alignment) {
  unsigned cigar_length;
  if(alignment.type == BWA_TYPE_NO_MATCH) cigar_length = 0;
  else if(!alignment.cigar) cigar_length = 1;
  else cigar_length = alignment.n_cigar;
  
  jcharArray java_cigar_operators = env->NewCharArray(cigar_length);
  if(java_cigar_operators == NULL) return NULL;
  jintArray java_cigar_lengths = env->NewIntArray(cigar_length);
  if(java_cigar_lengths == NULL) return NULL;
  
  if(alignment.cigar) {
    for(unsigned cigar_idx = 0; cigar_idx < (unsigned)alignment.n_cigar; ++cigar_idx) {
      jchar cigar_operator = "MIDS"[alignment.cigar[cigar_idx]>>14];
      jint  cigar_length = alignment.cigar[cigar_idx]&0x3fff;
      
      env->SetCharArrayRegion(java_cigar_operators,cigar_idx,1,&cigar_operator);
      if(env->ExceptionCheck()) return NULL;
      env->SetIntArrayRegion(java_cigar_lengths,cigar_idx,1,&cigar_length);
      if(env->ExceptionCheck()) return NULL;
    }
  }
  else {
    if(alignment.type != BWA_TYPE_NO_MATCH) {
      jchar cigar_operator = 'M';
      env->SetCharArrayRegion(java_cigar_operators,0,1,&cigar_operator);
      if(env->ExceptionCheck()) return NULL;
      env->SetIntArrayRegion(java_cigar_lengths,0,1,&read_length);
      if(env->ExceptionCheck()) return NULL;
    }
  }
  delete[] alignment.cigar;
    
  jclass java_alignment_class = env->FindClass("org/broadinstitute/sting/alignment/Alignment");
  if(java_alignment_class == NULL) return NULL;
  
  jmethodID java_alignment_constructor = env->GetMethodID(java_alignment_class, "<init>", "(IIZI[C[IILjava/lang/String;IIIII)V");
  if(java_alignment_constructor == NULL) return NULL;

  jstring java_md = env->NewStringUTF(alignment.md);
  if(java_md == NULL) return NULL;
  delete[] alignment.md;
  
  jobject java_alignment = env->NewObject(java_alignment_class,
                                          java_alignment_constructor,
                                          alignment.contig,
                                          alignment.pos,
                                          alignment.negative_strand,
                                          alignment.mapping_quality,
                                          java_cigar_operators,
                                          java_cigar_lengths,
                                          alignment.edit_distance,
                                          java_md,
                                          alignment.num_mismatches,
                                          alignment.num_gap_opens,
                                          alignment.num_gap_extensions,
                                          alignment.num_best,
                                          alignment.num_second_best);
  if(java_alignment == NULL) return NULL;      
  
  env->DeleteLocalRef(java_alignment_class);
  if(env->ExceptionCheck()) return NULL;

  return java_alignment;
}

static jstring get_configuration_file(JNIEnv* env, jobject configuration, const char* field_name) {
  jclass configuration_class = env->GetObjectClass(configuration);
  if(configuration_class == NULL) return NULL;

  jfieldID configuration_field = env->GetFieldID(configuration_class, field_name, "Ljava/io/File;");
  if(configuration_field == NULL) return NULL;

  jobject configuration_file = (jobject)env->GetObjectField(configuration,configuration_field); 

  jclass file_class = env->FindClass("java/io/File");
  if(file_class == NULL) return NULL;

  jmethodID path_extractor = env->GetMethodID(file_class,"getAbsolutePath", "()Ljava/lang/String;");
  if(path_extractor == NULL) return NULL;
  
  jstring path = (jstring)env->CallObjectMethod(configuration_file,path_extractor);
  if(path == NULL) return NULL;

  env->DeleteLocalRef(configuration_class);
  env->DeleteLocalRef(file_class);
  env->DeleteLocalRef(configuration_file);

  return path;
}

static void set_int_configuration_param(JNIEnv* env, jobject configuration, const char* field_name, BWA* bwa, int_setter setter) {
  jclass configuration_class = env->GetObjectClass(configuration);
  if(configuration_class == NULL) return;

  jfieldID configuration_field = env->GetFieldID(configuration_class, field_name, "Ljava/lang/Integer;");
  if(configuration_field == NULL) return;

  jobject boxed_value = env->GetObjectField(configuration,configuration_field);
  if(env->ExceptionCheck()) return;

  if(boxed_value != NULL) {
    jclass int_box_class = env->FindClass("java/lang/Integer");
    if(int_box_class == NULL) return;

    jmethodID int_extractor = env->GetMethodID(int_box_class,"intValue", "()I");
    if(int_extractor == NULL) return;

    jint value = env->CallIntMethod(boxed_value,int_extractor);
    if(env->ExceptionCheck()) return;

    if(value < 0) 
    {
      throw_config_value_exception(env,field_name,"cannot be set to a negative value");
      return;
    }

    (bwa->*setter)(value);

    env->DeleteLocalRef(int_box_class);
  }

  env->DeleteLocalRef(boxed_value);
  env->DeleteLocalRef(configuration_class);
}

static void set_float_configuration_param(JNIEnv* env, jobject configuration, const char* field_name, BWA* bwa, float_setter setter) 
{
  jclass configuration_class = env->GetObjectClass(configuration);
  if(configuration_class == NULL) return;

  jfieldID configuration_field = env->GetFieldID(configuration_class, field_name, "Ljava/lang/Float;");
  if(configuration_field == NULL) return;

  jobject boxed_value = env->GetObjectField(configuration,configuration_field);
  if(boxed_value != NULL) {
    jclass float_box_class = env->FindClass("java/lang/Float");
    if(float_box_class == NULL) return;

    jmethodID float_extractor = env->GetMethodID(float_box_class,"floatValue", "()F");
    if(float_extractor == NULL) return;

    jfloat value = env->CallFloatMethod(boxed_value,float_extractor);
    if(env->ExceptionCheck()) return;

    if(value < 0) 
    {
      throw_config_value_exception(env,field_name,"cannot be set to a negative value");
      return;
    }

    (bwa->*setter)(value);

    env->DeleteLocalRef(float_box_class);
  }

  env->DeleteLocalRef(boxed_value);
  env->DeleteLocalRef(configuration_class);
}

static void throw_config_value_exception(JNIEnv* env, const char* field_name, const char* message) 
{
  char* buffer = new char[strlen(field_name)+1+strlen(message)+1];
  sprintf(buffer,"%s %s",field_name,message);
  jclass sting_exception_class = env->FindClass("org/broadinstitute/sting/utils/StingException");
  if(sting_exception_class == NULL) return;
  env->ThrowNew(sting_exception_class, buffer);
  delete[] buffer;
}
