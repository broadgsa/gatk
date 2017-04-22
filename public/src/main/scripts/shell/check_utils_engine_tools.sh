#!/bin/sh

# Exit with an error if:
# - utils contains a reference to engine or tools
# - engine contains a reference to tools

sh -c \
    "grep -Rn \
      -e 'org.broadinstitute.gatk.tools' \
      -e 'org.broadinstitute.gatk.engine' \
      */*/src/*/*/org/broadinstitute/gatk/utils | \
      grep -v dependencyanalyzer && \
    grep -Rn \
      -e 'org.broadinstitute.gatk.tools' \
      */*/src/*/*/org/broadinstitute/gatk/engine" | \
  sed -e 's/:/:'$'\x1B\x5B\x35\x6d''/2' -e 's/$/'$'\x1B\x5B\x6d''/' | \
  grep gatk

RESULT=$?
if [[ ${RESULT} -eq 0 ]]; then
    echo "Fix the above errors. Do not import tools nor engine into the utils, and do not import tools into the engine." >&2
    exit 1
else
    exit 0
fi
