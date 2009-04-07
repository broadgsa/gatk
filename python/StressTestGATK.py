import farm_commands
import os.path
import sys
import getopt

defaultCommands = ['CountReads', 'Pileup']

def usage():
    print "Optional arguments:"
    print "  -f QUEUE   Farm jobs to QUEUE on LSF"
    print "  -c cmd1,cmd2   Walkers to execute, otherwise", ' '.join(defaultCommands)
    print "  -e             Ignore existing files", ' '.join(defaultCommands)
    
if __name__ == "__main__":
    opts = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:c:a:e", ["farm", "commands", "args", "ignoreExistingFiles"])
    except getopt.GetoptError:
        print sys.argv
        usage()
        sys.exit(2)

    farm_sub = False
    commandsList = defaultCommands
    ignoreExistingFiles = False
    extraArgs = ''

    for opt, arg in opts:
        if opt in ("-f", "--farm"):
            farm_sub = arg
        if opt in ("-c", "--commands"):
            commandsList = arg.split(',')
        if opt in ("-e", "--ignoreExistingFiles"):
            ignoreExistingFiles = True
        if opt in ("-a", "--args"):
            extraArgs = arg

    directory = args[1]
    for line in open(args[0]):
        lineParts = line.split()
        lane = lineParts[0].strip()
        
        ref = '/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta'
        if ( len(lineParts) > 1 ):
            ref = lineParts[1].strip()
            
        if not os.path.exists(lane):
            print 'Input SAM/BAM file: "', lane, '" does not exist, skipping...'
            continue        

        head, lane_filename = os.path.split(lane)
        filebase = os.path.splitext(lane_filename)[0]
    
        # convert the fasta
        for analysis in commandsList:
            output = os.path.join(directory, filebase + '.' + analysis + '.output')
            if ignoreExistingFiles or not os.path.exists(output):
                cmd = "java -jar ~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -T " + analysis + " -I " + lane + " -R " + ref + " -o " + output + " -l INFO " + extraArgs
                print cmd
                farm_commands.cmd(cmd, farm_sub, output)


