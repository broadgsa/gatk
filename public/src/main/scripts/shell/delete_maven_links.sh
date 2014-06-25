#!/bin/bash
#
# Delete symlinks to testdata/qscripts created by maven. This is necessary
# in cases where "mvn clean" fails to delete these links itself and exits
# with an error.
#
# Should be run from the root directory of your git clone.
#

find -L . -type l -name testdata -print0 | xargs -0 rm
find -L . -type l -name qscript -print0 | xargs -0 rm

exit 0
