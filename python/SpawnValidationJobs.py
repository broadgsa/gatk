#!/usr/bin/env python

import getopt, sys, os, string
from farm_commands import *

def picardCMD(name, **keywords ):
    cmd = name
    for key, value in keywords.iteritems():
        cmd += ' ' + key + "=" + str(value)
    return cmd

def spawnValidationJob( input_file, output_head, farm, maxErrors):
    validate_exe = "ValidateSAM"
    output_file = output_head + '.stdout'

    if regenExistingFiles or not os.path.exists(output_file):
        cmd_str = picardCMD( validate_exe, I=input_file, M=maxErrors )
        if farm == "":
             cmd_str += " > " + output_file
        cmd(cmd_str, farm, output_head, just_print_commands=justPrintCommands)
    
def usage():
    print "Required arguments:"
    print "  -d         Directory to grab all sam/bam files from"
    print
    print "Optional arguments:"
    print "  -f QUEUE   Farm jobs to QUEUE on LSF"
    print
    print "  -m MAXERRORS  Maximum number of errors to detect before aborting"
    print


def get_all_sam_files(dir):
    files = []
    
    for dirpath, dirnames, filenames in os.walk(dir):
        for filename in filenames:
            base, ext = os.path.splitext(filename)
            if ext.lower() in ['.sam', '.bam']:
                files.append( os.path.join( dirpath, filename ) )
            #print filename, base, ext

    return files

def output_filename( input_file ):
    parts = filter(lambda x: x.strip() <> '', input_file.split("/"))
    print parts
    return ".".join(parts) + ".validation" 

justPrintCommands = False
regenExistingFiles = False

if __name__ == "__main__":
    opts = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:f:m:r", ["dir","farm","maxErrors", "regenExistingFiles"])
    except getopt.GetoptError:
        print sys.argv
        usage()
        sys.exit(2)

    dir = ""
    mapper_str = "all"
    farm_sub = False
    maxErrors = 1000

    for opt, arg in opts:
        print opt, arg
        if opt in ("-d", "--dir"):
            dir = arg
        if opt in ("-f", "--farm"):
            farm_sub = arg
        if opt in ("-m", "--maxErrors"):
            maxErrors = arg
        if opt in ("-r", "--regenExistingFiles"):
            regenExistingFiles = True            
    if dir == "":
        usage()
        sys.exit(2)

    input_files = get_all_sam_files(dir)
    print 'Processing files: N=', len(input_files)
    for input_file in input_files:
        print '  ->', input_file
        
    for input_file in input_files:
        output_file = output_filename( input_file )
        print input_file, "=>", output_file
        spawnValidationJob( input_file, output_file, farm_sub, maxErrors )
        
        
        