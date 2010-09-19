import sys
import os
directories = list()
status = os.getenv("LSB_JOBEXIT_STAT")

for j in range(1,len(sys.argv)) :
    directories.append(sys.argv[j])

if ( status == "0" or status == 0):
    os.system("touch "+" ".join(directories))
