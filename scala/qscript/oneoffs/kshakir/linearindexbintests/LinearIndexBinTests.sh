#!/bin/sh

STING_HOME=/humgen/gsa-hpprojects/dev/kshakir/src/Sting_CountRod
TMP_DIR=/broad/shptmp/kshakir/CountRod
JOB_QUEUE=
STATUS_TO=
#JOB_QUEUE="-jobQueue week"
#STATUS_TO="-statusTo kshakir"

if [ "$1" == "debug" ]; then
    JAVA_DEBUG="-Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=n,address=8555"
    shift
fi

java $JAVA_DEBUG -Djava.io.tmpdir=${TMP_DIR} -jar ${STING_HOME}/dist/Queue.jar -jobPrefix QCountRodTest -S LinearIndexBinTests.scala -gatk ${STING_HOME}/dist/GenomeAnalysisTK.jar -jobQueue ${JOB_QUEUE} ${STATUS_TO} -bsub $@
