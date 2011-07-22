#!/bin/sh
export BWA_HOME="/humgen/gsa-scr1/hanna/src/bwa-trunk/bwa"
export JAVA_INCLUDE="/broad/tools/Linux/x86_64/pkgs/jdk_1.6.0_12/include -I/broad/tools/Linux/x86_64/pkgs/jdk_1.6.0_12/include/linux"
export TARGET_LIB="libbwa.so"
export EXTRA_LIBS="-lc -lz -lstdc++ -lpthread"
export LIBTOOL_COMMAND="g++ -shared -Wl,-soname,libbwa.so"
make 
