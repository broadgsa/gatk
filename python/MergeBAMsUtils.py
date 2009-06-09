import farm_commands
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import ValidateGATK
import picard_utils
import operator

MERGE_BIN = '/seq/software/picard/current/bin/MergeSamFiles.jar'
bam_ext = 'bam'
bam_index_ext = 'bai'

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

_abbrevs = [
    (1<<50L, 'P'),
    (1<<40L, 'T'), 
    (1<<30L, 'G'), 
    (1<<20L, 'M'), 
    (1<<10L, 'k'),
    (1, '')
    ]

def greek(size):
    """Return a string representing the greek/metric suffix of a size"""
    for factor, suffix in _abbrevs:
        if size > factor:
            break
    return '%.1f%s' % (float(size)/factor, suffix)

class MergeFilesSpec:
    def __init__(self, sources, group, merged_filename_base, path = '' ):
        self.sourceFiles = sources
        self.groupName = group
        self.merged_filename_base = merged_filename_base
        self.path = ''
    
    def __str__(self):
        return 'MergeFilesSpec: ' + self.group()
    
    def group(self):
        return self.groupName

    def sources(self):
        return self.sourceFiles

    def filename(self):
        if self.merged_filename_base <> None:
            return self.merged_filename_base + '.' + self.group()
        else:
            return self.group()

    def pprint(self):
        print '--------------------------------------------------------------------------------'
        print ' Population:         ', self.group()
        print ' Merged filename:    ', self.getMergedBAM()
        print ' N sources:          ', len(self.sources())
        print ' Sources:            ', self.sources()
        print ' Sizes:              ', self.sourceSizes(humanReadable=True)
        print ' Est. merged size:   ', greek(reduce(operator.__add__, self.sourceSizes(), 0))
        
    def setPath(self, path):
        self.path = path
        
    def getMergedBase(self):
        return os.path.join(self.path, self.filename())

    def getMergedBAM(self):
        return self.getMergedBase() + '.' + bam_ext

    def getMergedBAMIndex(self):
        return self.getMergedBase() + '.' + bam_ext + '.' + bam_index_ext
                
    def sourceSizes(self, humanReadable=False):
        sizes = map( os.path.getsize, self.sources() )
        if humanReadable:
            sizes = map(greek, sizes)
        return sizes
                
    def mergeCmd(self, mergeBin = None, MSD = True, useSamtools = False):
        if mergeBin == None:
            mergeBin = MERGE_BIN
            
        return picard_utils.mergeBAMCmd(self.getMergedBAM(), self.sources(), mergeBin, MSD = MSD, useSamtools = useSamtools)
        
    def getIndexCmd(self):
        return "samtools index " + self.getMergedBAM()
       
# Very general-purpose operation that returns merge specs given two lists of pairs:
# The first list, sources, contains superkey / sourceFile pairs.  
# The second list, groups, contains key / group pairs.
#
# This function walks over all source pairs, and for each superkey, it 
# looks for any key within groups contained within superkey.  If it finds it,
# it associates sourceFile with the merge group in groups.  
#
# The system requires that superkey match one (and only) one key in groups.  It also
# requires that each group string be unique.  The system can handle groups provided
# as a list of pairs of a dictionary.
#
# The function returns a list of MergeFileSpecs, one for each group with 
# at least one associated sourceFile.
#
def groupSources(sources, groups, merged_filename_base):
    if groups == None:
        return [MergeFilesSpec(map( lambda x: x[1], sources), '', merged_filename_base)]
    else:
        specs = dict()
        if type(groups) == list: groups = dict(groups)
        
        for superkey, sourceFile in sources:
            spec = None
            for key, group in groups.iteritems():
                #print 'Examining', superkey, key, group
                if superkey.find(key) <> -1:
                    if group in specs:
                        spec = specs[group]
                    else:
                        spec = MergeFilesSpec([], group, merged_filename_base)
                        specs[group] = spec
                    print 'Mapping', group, key, superkey, sourceFile
                    spec.sourceFiles.append(sourceFile)
            if spec == None:
                sys.exit('File contains an unknown superkey: ' + superkey)
        v = specs.values()
        v.sort(key = MergeFilesSpec.group)
        return v

import unittest
class TestMergeBAMsUtils(unittest.TestCase):
    def setUp(self):
        import cStringIO
        groupsString = """NA10846 ceu CEPH1
NA10847 ceu CEPH1
NA12144 ceu CEPH1
NA12145 ceu CEPH1
NA12146 yri CEPH1
NA12239 yri CEPH1
NA07029 ceu CEPH1
NA07019 ceu CEPH1
NA06994 ceu CEPH1
NA07000 ceu CEPH1
NA07022 ceu CEPH1
NA07056 ceu CEPH1
NA07048 ceu CEPH1
NA06991 ceu CEPH1
NA07034 ceu CEPH1
"""
        lanesString = """NA10846 30GA9AAXX 1 Paired CEPH 30GA9AAXX.1.observed_genotypes.geli
NA10846 30GA9AAXX 6 Paired CEPH 30GA9AAXX.6.observed_genotypes.geli
NA10847 30GA9AAXX 7 Paired CEPH 30GA9AAXX.7.observed_genotypes.geli
NA12146 30JLTAAXX 2 Paired CEPH 30JLTAAXX.2.observed_genotypes.geli
NA12239 30PNVAAXX 1 Paired CEPH 30PNVAAXX.1.observed_genotypes.geli
NA12144 30PYMAAXX 1 Paired CEPH 30PYMAAXX.1.observed_genotypes.geli
NA12146 30PYMAAXX 7 Paired CEPH 30PYMAAXX.7.observed_genotypes.geli
"""

        self.laneIDCounts = dict([["NA10846", 2], ["NA10847", 1], ["NA12146", 2], ["NA12239", 1], ["NA12144", 1]])

        pops = [line.strip() for line in cStringIO.StringIO(groupsString)]
        lanes = [line.strip() for line in cStringIO.StringIO(lanesString)]
        
        print pops
        print lanes
        
        self.ids2pop = [line.split()[0:2] for line in pops]
        self.ids2gelis = [[line.split()[0], line.split()[5]] for line in lanes]
        self.ids2ids = dict([[line.split()[0]] * 2 for line in lanes])

    def testPopGroups(self):
        specs = groupSources(self.ids2gelis, self.ids2pop, 'foo')
        print 'Specs', specs
        self.assertEqual(len(specs), 2)
        self.assertEqual(specs[0].group(), 'ceu')
        self.assertEqual(specs[1].group(), 'yri')
        
        ceu = specs[0]
        yri = specs[1]
        #print ceu.sources()
        self.assertEqual(len(ceu.sources()), 4)
        self.assertEqual(len(yri.sources()), 3)
        self.assertEqual(ceu.getMergedBAM(), 'foo.ceu.bam')
        self.assertEqual(ceu.getMergedBAMIndex(), 'foo.ceu.bam.bai')
        self.assertEqual(yri.getMergedBAM(), 'foo.yri.bam')
        self.assertEqual(yri.getMergedBAMIndex(), 'foo.yri.bam.bai')

    def testIDGroups(self):
        specs = groupSources(self.ids2gelis, self.ids2ids, 'foo')
        self.assertEqual(len(specs), 5)
        for spec in specs:
            print 'Spec', spec
            self.assertEqual(len(spec.sources()), self.laneIDCounts[spec.group()])
            self.assertEqual(spec.getMergedBAM(), 'foo.' + spec.group() + '.bam')

if __name__ == '__main__':
    unittest.main()