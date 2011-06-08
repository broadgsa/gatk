import sys
import os

in_list = sys.argv[1]
out_list = sys.argv[2]

inp = open(in_list)
out = open(out_list,'w')

for line in inp.readlines():
    path = line.strip()
    #print(path)
    tries = 0
    index = path.strip(".bam")+".bai"
    finished = path.rsplit("/",1)[0]+"/finished.txt"
    while ( not (os.path.exists(path) and os.path.exists(finished) and os.path.exists(index) ) and tries < 15 ):
        base = path.rsplit("/",1)[0]
        vers = base.rsplit("/",1)[1]
        base = base.rsplit("/",1)[0]
        bam = path.rsplit("/",1)[1]
        vers = "v%d" % (int(vers.lstrip("v"))+1)
        path = "%s/%s/%s" % (base,vers,bam)
        finished = "%s/%s/%s" % (base,vers,"finished.txt")
        index = path.strip(".bam")+".bai"
        #print(path)
        tries += 1
    if ( os.path.exists(path) and os.path.exists(finished) and os.path.exists(index) ):
        out.write(path+"\n")
    else:
        print("No filepath on the file system contains bam, index, and finished.txt for entry "+path)
