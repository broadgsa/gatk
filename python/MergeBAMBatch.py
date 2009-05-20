import farm_commands
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import ValidateGATK

MERGE_BIN = '/seq/software/picard/current/bin/MergeSamFiles.jar'
bam_ext = '.bam'

if __name__ == "__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
    parser.add_option("-d", "--dir", dest="output_dir",
                        type="string", default="./",
                        help="Output directory")
    parser.add_option("-i", "--ignoreExistingFiles", dest="ignoreExistingFiles",
                        action='store_true', default=False,
                        help="Ignores already written files, if present")

    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    directory = OPTIONS.output_dir
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    
    today = date.today()
    time_stamp = today.isoformat()
    
    for line in open(args[0]):
        s = line.split()
        if ( s <> [] and s[0] <> '#' ):
            merged_filename = s[0]
            output = os.path.join(directory, merged_filename + '.stdout')
            output_filename = os.path.join(directory, merged_filename + bam_ext)
            output_index = output_filename + ".bai"
            sources = reduce( operator.__add__, map( glob.glob, s[1:] ), [] )
            
            if OPTIONS.ignoreExistingFiles or not os.path.exists(output_filename):
                cmd = 'java -Xmx4096m -jar ' + MERGE_BIN + ' AS=true SO=coordinate O=' + output_filename + ' VALIDATION_STRINGENCY=SILENT ' + ' I=' + (' I='.join(sources))
                print cmd
                farm_commands.cmd(cmd, OPTIONS.farmQueue, output)

            if OPTIONS.ignoreExistingFiles or not os.path.exists(output_index):
                ValidateGATK.indexBAM(output_filename, OPTIONS.farmQueue)

