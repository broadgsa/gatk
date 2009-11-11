#!/bin/sh
export BWA_HOME="/Users/mhanna/src/bwa"
export JAVA_INCLUDE="/System/Library/Frameworks/JavaVM.framework/Headers"
export TARGET_LIB="libbwa.dylib"
export EXTRA_LIBS="-lc -lz -lsupc++"
export LIBTOOL_COMMAND="libtool -dynamic"
make 
