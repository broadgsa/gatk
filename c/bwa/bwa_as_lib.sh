CFLAGS="-g -Wall -O2 -m64"
BWA_INCLUDE="/Users/mhanna/src/bwa-0.5.3"
BWA_SRC="/Users/mhanna/src/bwa-0.5.3"

g++ $CFLAGS -c -I$BWA_INCLUDE -I/System/Library/Frameworks/JavaVM.framework/Headers org_broadinstitute_sting_alignment_bwa_BWACAligner.cpp
g++ $CFLAGS -c -I$BWA_INCLUDE bwa_gateway.cpp -o bwa_gateway.o
gcc $CFLAGS -c  $BWA_SRC/bntseq.c -o bntseq.o
gcc $CFLAGS -c $BWA_SRC/bwase.c -o bwase.o
gcc $CFLAGS -c $BWA_SRC/bwt.c -o bwt.o
gcc $CFLAGS -c $BWA_SRC/bwtaln.c -o bwtaln.o
gcc $CFLAGS -c $BWA_SRC/bwtgap.c -o bwtgap.o
gcc $CFLAGS -c $BWA_SRC/bwtio.c -o bwtio.o
gcc $CFLAGS -c $BWA_SRC/bwaseqio.c -o bwaseqio.o
gcc $CFLAGS -c $BWA_SRC/cs2nt.c -o cs2nt.o
gcc $CFLAGS -c $BWA_SRC/kstring.c -o kstring.o
gcc $CFLAGS -c $BWA_SRC/stdaln.c -o stdaln.o
gcc $CFLAGS -c $BWA_SRC/utils.c -o utils.o
g++ -dynamiclib -o libbwa.dylib org_broadinstitute_sting_alignment_bwa_BWACAligner.o bwa_gateway.o bntseq.o bwase.o bwt.o bwtaln.o bwtgap.o bwtio.o bwaseqio.o cs2nt.o kstring.o stdaln.o utils.o -framework JavaVM -lz