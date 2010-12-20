#!/bin/sh

DIRECTORY=.

find $DIRECTORY -name CountRod.out -exec grep -H "Max Memory" {} \; | sort | sed -e 's!'$DIRECTORY'/\(.*\)_\([.0-9]*\)mfpb_\(.*\)g_run\(.*\)/CountRod.out: * Max Memory *: *\([.0-9]*\) MB!\1\t\2\t\3\t\4\t\5!' > $DIRECTORY/memory.txt

find $DIRECTORY -name CountRod.out -exec grep -H "CPU" {} \; | sort | sed -e 's!'$DIRECTORY'/\(.*\)_\([.0-9]*\)mfpb_\(.*\)g_run\(.*\)/CountRod.out: * CPU time *: *\([.0-9]*\) sec.!\1\t\2\t\3\t\4\t\5!' > $DIRECTORY/cpu.txt

find $DIRECTORY -name \*.done -o -name \*.fail | sort | sed -e 's!'$DIRECTORY'/\(.*\)_\([.0-9]*\)mfpb_\(.*\)g_run\(.*\)/\.CountRod.txt\.*\(.*\)!\1\t\2\t\3\t\4\t\5!' > $DIRECTORY/success.txt

TAB="	"
echo "set${TAB}max_features_per_bin${TAB}memory_limit_gb${TAB}run_number${TAB}max_memory_mb${TAB}cpu_s${TAB}job_success" > $DIRECTORY/mfpb.txt
paste $DIRECTORY/memory.txt $DIRECTORY/cpu.txt $DIRECTORY/success.txt | cut -f 1-5,10,15 >> $DIRECTORY/mfpb.txt
