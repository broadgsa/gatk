import sys
import os
import subprocess
import shlex

from optparse import OptionParser

def parseInput(fList,ignoreExt = None):
    inNames = []
    for ele in fList:
        if ( ignoreExt != None and ele.endswith(ignoreExt) ):
            inFileNames.append(ele)
        if ( os.path.exists(ele) ):
            for line in open(ele).readlines():
                inNames.append(line.strip())
    return inNames

def bamsWithSamples(bamList):
    cmdbase = "samtools view -H %s | grep SM | tr '\\t' '\\n' | grep SM | sed 's/SM://g' | uniq"
    sam2bam = dict()
    for bf in bamList:
        if ( not os.path.exists(bf) ):
            raise IOError("Bam file "+bf+" does not exist")
        cmd = cmdbase % bf
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        proc.wait()
        stdout = proc.stdout.readlines()
        if ( len(stdout) > 1 ):
            raise RuntimeError("The bam file "+bf+" contains multiple different sample entries")
        sm = stdout[0].strip()
        sam2bam[sm]=bf
    return sam2bam

def runMain(opt,arg):
    bamFiles = bamsWithSamples(parseInput(opt.bam_files,"bam"))
    caseNames = set(parseInput(opt.cases))
    controlNames = set(parseInput(opt.controls))
    output = open(opt.output,'w')
    output.write("samples:")
    sample_base = "- id: %s\n  properties:\n    cohort: %s\n    bam: %s"
    for s in bamFiles.keys():
        cc = "Unknown"
        if ( s in caseNames ):
            cc = "case"
        if ( s in controlNames ):
            cc = "control"
        output.write("\n" + sample_base % (s,cc,bamFiles[s]))
    output.close()

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-I","--bams",dest="bam_files",help="the bam files, as multiple arguments or a simple newline-delimited file",action="append")
    parser.add_option("-A","--case",dest="cases",action="append",help="A list of the case samples, as multiple arguments or a simple newline-delimited file")
    parser.add_option("-O","--control",dest="controls",action="append",help="A list of the control samples, multiple arguments or a newline-delimited file")
    parser.add_option("-o","--out",dest="output",action="store",help="Name of the output metadata file to write to")
    (options,args) = parser.parse_args()
    runMain(options,args)

if __name__ == "__main__":
    main()
