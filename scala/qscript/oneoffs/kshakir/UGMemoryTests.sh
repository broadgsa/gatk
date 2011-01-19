#!/bin/sh

STING_HOME=/humgen/gsa-hpprojects/dev/kshakir/src/Sting_patches
TMP_DIR=/broad/shptmp/kshakir
JOB_QUEUE=gsa

if [ "$1" == "debug" ]; then
    JAVA_DEBUG="-Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=n,address=8555"
    shift
fi

java $JAVA_DEBUG -Djava.io.tmpdir="$TMP_DIR" -jar "$STING_HOME"/dist/Queue.jar -jobPrefix QTest -S "$STING_HOME"/scala/qscript/kshakir/UGMemoryTests.scala -Y UGMemoryTests.yaml -gatk "$STING_HOME"/dist/GenomeAnalysisTK.jar -jobQueue $JOB_QUEUE $@
