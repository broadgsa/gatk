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


#ifndef JNI_COMMON_H
#define JNI_COMMON_H

/*#define SINGLE_THREADED_ONLY 1*/
#include <jni.h>
/*#define ENABLE_ASSERTIONS 1*/
#ifdef SINGLE_THREADED_ONLY
#define DO_PROFILING 1
#endif
/*#define DEBUG0_1 1*/
/*#define DEBUG3 1*/
/*#define DUMP_TO_SANDBOX 1*/


/*#define DIRECT_ACCESS_TO_JAVA_HEAP_MEMORY 1*/

#ifdef DIRECT_ACCESS_TO_JAVA_HEAP_MEMORY
//Gets direct access to Java arrays
#define GET_BYTE_ARRAY_ELEMENTS env->GetPrimitiveArrayCritical
#define RELEASE_BYTE_ARRAY_ELEMENTS env->ReleasePrimitiveArrayCritical
#define JNI_RO_RELEASE_MODE JNI_ABORT
#define GET_DOUBLE_ARRAY_ELEMENTS env->GetPrimitiveArrayCritical
#define RELEASE_DOUBLE_ARRAY_ELEMENTS env->ReleasePrimitiveArrayCritical

#else
//Likely makes copy of Java arrays to JNI C++ space
#define GET_BYTE_ARRAY_ELEMENTS env->GetByteArrayElements
#define RELEASE_BYTE_ARRAY_ELEMENTS env->ReleaseByteArrayElements
#define JNI_RO_RELEASE_MODE JNI_ABORT
#define GET_DOUBLE_ARRAY_ELEMENTS env->GetDoubleArrayElements
#define RELEASE_DOUBLE_ARRAY_ELEMENTS env->ReleaseDoubleArrayElements

#endif		//ifdef DIRECT_ACCESS_TO_JAVA_HEAP_MEMORY

#endif  //ifndef JNI_COMMON_H
