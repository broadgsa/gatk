import farm_commands
import os.path
import sys
from optparse import OptionParser
import string
import re
import glob

#lanes = ["30JW3AAXX.6", "30KRNAAXX.1", "30KRNAAXX.6", "30PYMAAXX.5"]
#idsList = ['NA12843', 'NA19065', 'NA19064', 'NA18637']

lanes = ["30JW3AAXX.6", "30PYMAAXX.5", "30PNUAAXX.8", "30PPJAAXX.5"]
idsList = ['NA12843', 'NA18637', "NA19058", "NA12842"]
ids = dict(zip(lanes, idsList))
gatkPath = "~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar"
ref = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
analysis = "CombineDuplicates"

def getPicardPath(lane, picardRoot = '/seq/picard/'):
    flowcell, laneNo = lane.split('.')
    filePat = os.path.join(picardRoot, flowcell, '*', laneNo, '*')
    dirs = glob.glob(filePat)
    print dirs
    if len(dirs) > 1:
        system.exit("Bad lane -- too many directories matching pattern " + filePat)
    return dirs[0]

def getReferenceGenotypeFileFromConcordanceFile(concordFile):
    # REFERENCE_GENOTYPES=/seq/references/reference_genotypes/hapmap/Homo_sapiens_assembly18/NA19058.geli
    p = re.compile('REFERENCE_GENOTYPES=([/.\w]+)')
    for line in open(concordFile):
        match = p.search(line)
        print 'Match is', line, match
        if match <> None:
            return match.group(1)
    return None

def concord(options, geli, output, genotypeFile):
    return ("java -jar /seq/software/picard/current/bin/CollectGenotypeConcordanceStatistics.jar OPTIONS_FILE=%s INPUT=%s OUTPUT=%s REFERENCE_GENOTYPES=%s MINIMUM_LOD=5.0" % ( options, geli, output, genotypeFile ) )

def readPicardConcordance(file):
    p = re.compile('HOMOZYGOUS_REFERENCE|HETEROZYGOUS|HOMOZYGOUS_NON_REFERENCE')
# CATEGORY    OBSERVATIONS        AGREE   DISAGREE        PCT_CONCORDANCE
# HOMOZYGOUS_REFERENCE    853     853     0       1
# HETEROZYGOUS    416     413     3  0.992788
# HOMOZYGOUS_NON_REFERENCE        235     231     4       0.982979
    types = [str, int, int, int, float]
    def parse1(line):
        return [f(x) for f, x in zip(types, line.split())]
    data = [parse1(line) for line in open(file) if p.match(line) <> None]
    return data
