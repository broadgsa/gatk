import farm_commands
import os.path
import sys
from optparse import OptionParser

defaultCommands = ['CountReads', 'Pileup']

def indexBAM(reads):
    cmd = "samtools index " + reads
    farm_commands.cmd(cmd, None, None)    

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
                        action='store_true', default=True,
                        help="Ignores the existing files")
    parser.add_option("-a", "--rebuildAllFiles", dest="rebuildAllFiles",
                        action='store_true', default=False,
                        help="If provided, all intermediate files (BAM and pileups) will be regenerated")

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
        readsIndex = os.path.join(OPTIONS.dataDir, filebase + '.selected.bam.bai')
        subBAM = os.path.join(OPTIONS.dataDir, filebase + '.selected.bam')
        pileup = os.path.join(OPTIONS.dataDir, filebase + '.selected.pileup')
        validationOutput = os.path.join(OPTIONS.outputDir, filebase + '.validate.output')        

        #print 'reads', reads
        if not os.path.exists(reads) or OPTIONS.rebuildAllFiles:
            farm_commands.cmd("ln -s " + os.path.abspath(originalReads) + " " + reads)

        if not os.path.exists(readsIndex) or OPTIONS.rebuildAllFiles:
            indexBAM(reads)
            
        if not os.path.exists(subBAM) or OPTIONS.rebuildAllFiles:
            if region == '*':
                farm_commands.cmd("ln -s " + os.path.abspath(reads) + " " + subBAM)
                farm_commands.cmd("ln -s " + os.path.abspath(readsIndex) + " " + subBAM+'.bai')
            else:
                cmd = "samtools view -b " + reads + " " + region + " > " + subBAM
                farm_commands.cmd(cmd, None, None)
 
        if not os.path.exists(pileup) or OPTIONS.rebuildAllFiles:
            cmd = "samtools pileup -cf " + ref + " " + subBAM + " > " + pileup
            farm_commands.cmd(cmd, None, None)
 
        if not os.path.exists(validationOutput) or OPTIONS.ignoreExistingFiles:
            analysis = "ValidatingPileup"
            cmd = "java -jar ~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -T " + analysis + " -I " + subBAM + " -R " + ref + " -l INFO -S SILENT -U -B pileup SAMPileup " + pileup
            print cmd
            farm_commands.cmd(cmd, OPTIONS.farmQueue, outputFile=validationOutput)

if __name__ == "__main__":
    main()
