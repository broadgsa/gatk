#!/bin/sh
export BWA_HOME="${PWD}/bwasvn47"
export JAVA_INCLUDE="${JAVA_HOME}/include"
export JAVA_PLATFORM_INCLUDE="${JAVA_HOME}/include/darwin"
export TARGET_LIB="libbwa.dylib"
export EXTRA_LIBS="-lc -lz -lstdc++"
export LIBTOOL_COMMAND="libtool -dynamic"
make 
