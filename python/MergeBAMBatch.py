import farm_commands
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import ValidateGATK
import picard_utils

MERGE_BIN = '/seq/software/picard/current/bin/MergeSamFiles.jar'
bam_ext = '.bam'

def readNAIdMap(NAIdFile):
    m = dict()
    for data in [line.split() for line in open(NAIdFile)]:
        naid, pop = data[0:2]
        print naid, ' => ', pop
        assert naid not in m
        m[naid] = pop
    print 'Read NAID->population map'
    print 'Contains', len(m), 'id -> population mappings'
    print 'Distinct populations:', picard_utils.unique(m.values())
    return m

class MergeFilesSpec:
    def __init__(self, sources, pop, merged_filename_base ):
        self.sourceFiles = sources
        self.pop = pop
        self.merged_filename_base = merged_filename_base

    def sources(self):
        return self.sourceFiles

    def filename(self):
        return self.merged_filename_base + '.' + self.pop

def splitSourcesByPopulation(allSources, merged_filename_base, NAID2Pop):
    if NAID2Pop == None:
        return [MergeFilesSpec(allSources, '', merged_filename_base)]
    else:
        specs = dict()
        for source in allSources:
            spec = None
            for naid, pop in NAID2Pop.iteritems():
                if source.find(naid) <> -1:
                    if pop in specs:
                        spec = specs[pop]
                    else:
                        spec = MergeFilesSpec([], pop, merged_filename_base)
                        specs[pop] = spec
                    #print 'Mapping', source, naid, pop
                    spec.sourceFiles.append(source)
            if spec == None:
                sys.exit('File contains an unknown NAID: ' + source)
        return specs.values()

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
    parser.add_option("-m", "--mergeBin", dest="mergeBin",
                        type="string", default=MERGE_BIN,
                        help="Path to merge binary")
    parser.add_option("-n", "--naIDPops", dest="NAIDS2POP",
                        type="string", default=None,
                        help="Path to file contains NAID POP names.  If provided, input files will be merged by population")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    directory = OPTIONS.output_dir
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    
    NAID2Pop = None
    if OPTIONS.NAIDS2POP <> None:
        NAID2Pop = readNAIdMap(OPTIONS.NAIDS2POP)
    
    today = date.today()
    time_stamp = today.isoformat()
    
    for line in open(args[0]):
        s = line.split()
        if ( s <> [] and s[0] <> '#' ):
            merged_filename_base = s[0]
            allSources = reduce( operator.__add__, map( glob.glob, s[1:] ), [] )
            print 'Merging info:'
            for mergeFilesSpec in splitSourcesByPopulation(allSources, merged_filename_base, NAID2Pop):
                print '-----'
                print ' Population', mergeFilesSpec.pop
                print ' Filename', mergeFilesSpec.filename()
                print ' N sources', len(mergeFilesSpec.sources())
                print ' sources', mergeFilesSpec.sources()

                output = os.path.join(directory, mergeFilesSpec.filename() + '.stdout')
                output_filename = os.path.join(directory, mergeFilesSpec.filename() + bam_ext)
                output_index = output_filename + ".bai"
                
                jobid = None
                if OPTIONS.ignoreExistingFiles or not os.path.exists(output_filename):
                    #cmd = 'java -Xmx4096m -jar ' + OPTIONS.mergeBin + ' MSD=true AS=true SO=coordinate O=' + output_filename + ' VALIDATION_STRINGENCY=SILENT ' + ' I=' + (' I='.join(sources))
                    cmd = picard_utils.mergeBAMCmd(output_filename, mergeFilesSpec.sources(), OPTIONS.mergeBin)
                    print cmd
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, output)
    
                if OPTIONS.ignoreExistingFiles or not os.path.exists(output_index):
                    cmd = "samtools index " + output_filename
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, output, waitID = jobid)

