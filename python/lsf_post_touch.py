import sys
import os
status = sys.argv[1]
directories = list()

for j in range(2,len(sys.argv)) :
    directories.append(sys.argv[j])

if ( status == "0" ):
    os.system("touch "+" ".join(directories))
