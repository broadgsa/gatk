import farm_commands
import os.path
import sys
from optparse import OptionParser
import string

defaultCommands = ['CountReads', 'Pileup']

def indexBAM(reads, queue=None):
    cmd = "samtools index " + reads
    farm_commands.cmd(cmd, queue, None)    

def readIntervalFile(file):
    def notHeaderP(line):
        return line[0][0] <> '@'
    lines = filter(notHeaderP, map( string.split, [line for line in open(file)]))
    
    byContig = {}
    for entry in lines:
        elt = [int(entry[1]), int(entry[2])]
        if entry[0] in byContig:
            byContig[entry[0]].append(elt)
        else:
            byContig[entry[0]] = [elt]
    
    #print byContig
    return byContig

def isIntervalFile(filename):
    return os.path.splitext(filename)[1] == '.interval_list'

def filterByInterval(intervalMap, unfilteredPileupFile, filteredPileupFile):
    print 'Intervals', intervalMap

    def maybeWrite1(line, offsetMap, out):
        debug = False
        parts = line.split()
        contig = parts[0]
        pos = int(parts[1])
        if debug: print 'pileup', contig, pos, line,

        if contig in intervalMap:
            intervals = intervalMap[contig]
            
            while (offsetMap[contig] < len(intervals)):
                   offset = offsetMap[contig]
                   #print intervals[offset]
                   curStart, curStop = intervals[offset][0:2]
                   if debug: print contig, pos, curStart, curStop
                   if pos >= curStart:
                       if pos <= curStop:
                           #print line,
                           print >> out, line,
                           break
                       else:
                           if debug: print 'Progressing', contig, pos, curStart, curStop
                           offsetMap[contig] += 1
                   else:
                       break
        return offsetMap
        
    def allDone( offsetMap, intervalMap ):
        def oneDone( contig ):
            return offsetMap[contig] >= len(intervalMap[contig])
        return all( map(oneDone, offsetMap.iterkeys()) )
            
    i = 0
    offsetMap = dict([(key, 0) for key in intervalMap.keys()])
    out = open(filteredPileupFile, 'w')
    for line in open(unfilteredPileupFile):
        offsetMap = maybeWrite1(line, offsetMap, out)
        i += 1
        if i % 10000 == 0 and allDone(offsetMap, intervalMap):
            break
    out.close()
                    
    

def main():    
    global OPTIONS, ROOT

    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file", dest="validateInfoFile",
                        type="string", default=None,
                        help="File name of the validation set")
    parser.add_option("-d", "--dataDir", dest="dataDir",
                        type="string", default=None,
                        help="Directory to the data files")
    parser.add_option("-o", "--outputDir", dest="outputDir",
                        type="string", default=None,
                        help="Directory to put the output")
    parser.add_option("-q", "--farmQueue", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to submit jobs to.  Leave blank for local processing")
    parser.add_option("-e", "--ignoreExistingFiles", dest="ignoreExistingFiles",
                        action='store_true', default=False,
                        help="Ignores the existing files")
    parser.add_option("-a", "--rebuildAllFiles", dest="rebuildAllFiles",
                        action='store_true', default=False,
                        help="If provided, all intermediate files (BAM and pileups) will be regenerated")
    parser.add_option("-n", "--dontValidate", dest="dontValidate",
                        action='store_true', default=False,
                        help="If provided, won't try to validate, just make sure the validation data files are ready")
    parser.add_option("-g", "--gatk", dest="gatkPath",
                        type="string", default="~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar",
                        help="Path to the gatk")
    parser.add_option("-x", "--args", dest="extraArgs",
                        type="string", default="",
                        help="Extra args to pass to the GATK")
    parser.add_option("-p", "--justPrint", dest="justPrint",
                        action='store_true', default=False,
                        help="Don't actually run GATK, just setup data files")



    (OPTIONS, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")

    if OPTIONS.validateInfoFile == None:
        parser.error("f option is required")
    if OPTIONS.dataDir == None:
        parser.error("d option is required")
    if OPTIONS.outputDir == None:
        parser.error("o option is required")

    if not os.path.exists(OPTIONS.outputDir):
        os.mkdir(OPTIONS.outputDir)
    if not os.path.exists(OPTIONS.dataDir):
        os.mkdir(OPTIONS.dataDir)
 
    for line in open(OPTIONS.validateInfoFile):
        if line.strip() == "" or line[0] == '#':
            # comment line
            continue
    
        # line is of the format:
        # <reads.bam> <ref.fasta> <processingRegion>
        originalReads, ref, region = line.split()
        
        originalReads = os.path.expanduser(originalReads)
        ref   = os.path.expanduser(ref)

        if not os.path.exists(originalReads) or not os.path.exists(ref):
            print 'Input files do not exist!', originalReads, ref
            sys.exit(1)        

        head, read_filename = os.path.split(originalReads)
        filebase = os.path.splitext(read_filename)[0]

        reads = os.path.join(OPTIONS.dataDir, read_filename)
        readsIndex = reads + '.bai'
        subBAM = os.path.join(OPTIONS.dataDir, filebase + '.selected.bam')
        pileup = os.path.join(OPTIONS.dataDir, filebase + '.selected.pileup')
        validationOutput = os.path.join(OPTIONS.outputDir, filebase + '.validate.output')        

        #print 'reads', reads
        if not os.path.exists(reads) or OPTIONS.rebuildAllFiles:
            farm_commands.cmd("ln -s " + os.path.abspath(originalReads) + " " + reads)

        if not os.path.exists(readsIndex) or OPTIONS.rebuildAllFiles:
            indexBAM(reads)
            
        if not os.path.exists(subBAM) or OPTIONS.rebuildAllFiles:
            if region == '*' or isIntervalFile(region):
                farm_commands.cmd("ln -s " + os.path.abspath(reads) + " " + subBAM)
            else:
                cmd = "samtools view -b " + reads + " " + region + " > " + subBAM
                farm_commands.cmd(cmd, None, None)

        subBAMIndex =  subBAM+'.bai'
        if not os.path.exists(subBAMIndex) or OPTIONS.rebuildAllFiles:
            if region == '*' or isIntervalFile(region):
                farm_commands.cmd("ln -s " + os.path.abspath(readsIndex) + " " + subBAMIndex)
            else:
                indexBAM(subBAM)

        if not os.path.exists(pileup) or OPTIONS.rebuildAllFiles:
            if isIntervalFile(region):
                filteredPileup = pileup + ".unfiltered"
                if not os.path.exists(filteredPileup) or OPTIONS.rebuildAllFiles:
                    cmd = "samtools pileup -cf " + ref + " " + subBAM + " > " + filteredPileup
                    #farm_commands.cmd(cmd, None, None)
                filterByInterval(readIntervalFile(region), filteredPileup, pileup)
            else:
                cmd = "samtools pileup -cf " + ref + " " + subBAM + " > " + pileup
                farm_commands.cmd(cmd, None, None)
 
        if not OPTIONS.dontValidate and (not os.path.exists(validationOutput) or OPTIONS.ignoreExistingFiles):
            print validationOutput, 'does not exist'
            analysis = "ValidatingPileup"
            cmd = "java -ea -Xmx1024m -jar " + OPTIONS.gatkPath + " -T " + analysis + " -I " + subBAM + " -R " + ref + " -l INFO -S SILENT -U " + OPTIONS.extraArgs

            if isIntervalFile(region):
                cmd += " --intervals " + region

            cmd += " -B pileup,SAMPileup," + pileup
            print cmd
            farm_commands.cmd(cmd, OPTIONS.farmQueue, outputFile=validationOutput, just_print_commands=OPTIONS.justPrint)

if __name__ == "__main__":
    main()
