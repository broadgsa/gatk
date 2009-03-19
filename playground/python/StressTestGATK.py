import farm_commands
import os.path
import sys
import getopt

def usage():
    print "Optional arguments:"
    print "  -f QUEUE   Farm jobs to QUEUE on LSF"

if __name__ == "__main__":
    opts = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:", ["farm"])
    except getopt.GetoptError:
        print sys.argv
        usage()
        sys.exit(2)

    farm_sub = False

    for opt, arg in opts:
        if opt in ("-f", "--farm"):
            farm_sub = arg

    for line in open(sys.argv[1]):
        lane = line.strip()
        head, lane_filename = os.path.split(lane)
        filebase = os.path.splitext(lane_filename)[0]
    
        # convert the fasta
        for analysis in ['CountLoci', 'Pileup']:
            output = filebase + '.' + analysis + '.output'
            if not os.path.exists(output):
                cmd = "java -jar ~/dev/GenomeAnalysisTK/trunk/playground/java/dist/GenomeAnalysisTK.jar T=" + analysis + " I= " + lane + " R= /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
                print cmd
                farm_commands.cmd(cmd, farm_sub, output, just_print_commands=True)


