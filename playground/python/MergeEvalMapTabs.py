#!/usr/bin/env python

# merges .tab files produced by EvalMapping.py into a tabular format

#from __future__ import print_function # python 2.6 / 3.0 print as a function
import os,sys,re
from rpy2 import *

def ordered_eval_files(aligner, mutation_type, file_ext):
    tab_files = [f for f in os.listdir(".") if f.endswith("."+aligner+"."+file_ext) and (mutation_type in f or "None_den." in f)]
    
    tabs = []
    for f in tab_files:
        num = int(re.search(r'den\.(\d\d?)_', f).groups()[0])
        if "None_den" in f:
            num = 0
        tabs.append( (num,f) )
    tabs.sort()
    return tabs

def anchor_hist(filename):
    fin = open(filename)
    offsets = []
    read_poses = []
    
    for line in fin:
        id = None
        fields = line.split()
        if len(fields) > 0 and ".read." in fields[0]:
            id = line.split()[0]
            read_poses.append(id.split(".")[-1])
            offsets.append(id.split(":")[1].split(".")[0])

        #if "read name :" in line:
        #    id = line.split()[3]
            
    import rpy2.robjects as robj
    if len(read_poses):
        print len(read_poses)
        return robj.r.hist(robj.FloatVector(read_poses), breaks=77, plot=False)
    else:
        print 0
        return robj.IntVector([])

def estimateCPUTime(filename):
    times = []
    if os.path.exists(filename):
        regex = re.compile('CPU time\s*:\s+([\d\.]+) sec')
        #print ('Filename', filename)
        def extractTime1(line):
            sobj = regex.search(line)
            if sobj <> None:
                time = float(sobj.group(1))
                #print (line, sobj, sobj.groups(), time)
                return time
            else:
                return None
        times = [time for time in map(extractTime1, open(filename)) if time <> None]
        #print('times', times)
        
    return max(times)

def print_eval_tables(aligner, mut_type, num_var, filename, min_qual_cutoff=0):
    filebase = os.path.splitext(filename)[0]
    cpuTime = estimateCPUTime(filebase + '.stdout')

    lines = open(filename).readlines()
    data = None
    if len(lines) > 1:
        da_line = lines[0]
        for line in lines:
            line = line.rstrip()
            fields = line.split("\t")
            try:
                qual_cutoff = int(fields[0])
            except:
                qual_cutoff = 1000
            if qual_cutoff <= min_qual_cutoff:
                break
            else:
                da_line = line
            
        #print ("%2d\t" % num_var, file=sys.stderr),
        data = [aligner, mut_type, min_qual_cutoff, num_var] + map(int, da_line.split()) + [cpuTime]
        print "%s\t%s\t%2d\t%s\t%.2f" % (aligner, mut_type, num_var, da_line, cpuTime)
    return data

def plot_anchor_hists():
    import rpy2.robjects as robj
    #robj.r.par( mfrow = robj.RVector([2,2]) )
    #robj.r('X11(width=24, height=15)')
    robj.r('png("lotto_histo.png", width=1500, height=900)')
    robj.r('par(mfrow=c(4,6 ), cin=c(3,3))')
    
    for (aligner, cutoff) in (("maq", 30), ("ilt",-1), ("merlin",-1), ("swmerlin",-1)):
        print aligner, ">"+str(cutoff)
        
        num_file = ordered_eval_files(aligner, mut_type, file_ext="eval")
        #print (*num_file, sep="\n")
        for num, file in num_file:
            
            h = anchor_hist( file )
            title = "{aligner} {num} bp insert".format(aligner=aligner, num=num)
            robj.r.plot(h, xlab="Anchor position", main=title, col="blue", ylim=robj.FloatVector([0,40]), xlim=robj.FloatVector([0,90]))

    robj.r("dev.off()")
    #raw_input("Press any key to exit...")


def analyzeData( data, mut_type ):
    import rpy2.robjects as robj
    # print name / data name
    aligners = [('MAQ', 'maq>30'), ('Merlin', 'swmerlin>-1'), ('ILT', 'ilt>-1'), ('BWA', 'bwa32>30')] 

    def selectData( key ):
        def keepDataP( data ):
            return data[0] == key and data[1] == mut_type
        return filter( keyDataP, data )

    def placeRatePlot():
        robj.r('png("x.png", width=1500, height=900)')
        x = robj.FloatVector(mutDensity)
        y = robj.FloatVector(placeRate)
        robj.r.plot(x, y, xlab="Mutation density / size", col="blue")
        robj.r("dev.off()")

    return 0

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "merge_tabs.py MUTATION_TYPE (e.g. INSERTION, DELETION, SNP)"
        sys.exit(1)
    else:
        mut_type = sys.argv[1]

    if False:
        plot_anchor_hists()
    else:
        data = []
        for (aligner, cutoff) in (("maq", 30), ("maq", 0), ("maq",-1), ("ilt",-1), ("merlin",-1), ("swmerlin",-1), ("bwa", 30), ("bwa", 0), ("bwa", -1), ("bwa32", 30), ("bwa32", 0), ("bwa32", -1)):
            name = aligner + ">" + str(cutoff)
            #print (aligner, ">"+str(cutoff))
            
            num_file = ordered_eval_files(aligner, mut_type, file_ext="tab")
            #num_file = ordered_eval_files(aligner, mut_type, file_ext="eval")
            for num, file in num_file:
                datum = print_eval_tables( name, mut_type, num, file, cutoff )
                if datum <> None:
                    data.append(datum)
                #print
    
        analyzeData( data, mut_type )

