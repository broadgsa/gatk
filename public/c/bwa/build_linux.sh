#!/bin/sh
export BWA_HOME="${PWD}/bwasvn47"
export JAVA_INCLUDE="${JAVA_HOME}/include"
export JAVA_PLATFORM_INCLUDE="${JAVA_HOME}/include/linux"
export TARGET_LIB="libbwa.so"
export EXTRA_LIBS="-lc -lz -lstdc++ -lpthread"
export LIBTOOL_COMMAND="g++ -shared -Wl,-soname,libbwa.so"
make 
