# import farm_commands2
import os.path
import sys
from optparse import OptionParser
import glob
import operator
import itertools
import re
import vcfReader
import string

def average(l):
    sum = reduce(operator.add, l, 0)
    return sum / len(l)

def printHeaderSep():
    print
    print ''.join(['-'] * 80)

class Sample:
    def __init__(self, name):
        self.name = name
        self.rawBases = 0
        self.mappedBases = 0
        self.nSNPs = 0
        self.nIndels = 0
        
    def getName(self): return self.name
    def getNSNPs(self): return self.nSNPs
    def getNIndels(self): return self.nIndels
        
    def __str__(self):
        return '[%s rawBases=%d mappedBases=%d percentMapped=%.2f nSNPs=%d nIndels=%d]' % (self.name, self.rawBases, self.mappedBases, (self.mappedBases * 100.0) / max(self.rawBases,1), self.nSNPs, self.nIndels)
    __repr__ = __str__
    
def flatFileIterator(file, fields = None, skip = 0):
    count = 0
    for line in open(file):
        count += 1
        if count > skip:
            s = map(string.strip, line.split('\t'))
            if ( fields != None ):
                s = map(lambda field: s[field], fields)
    
            if len(s) == 1: s = s[0]
            yield s

# 1.  FASTQ_FILE, path to fastq file on ftp site  
# 2.  MD5, md5sum of file
# 3.  RUN_ID, SRA/ERA run accession
# 4.  STUDY_ID, SRA/ERA study accession   
# 5.  STUDY_NAME, Name of stury   
# 6.  CENTER_NAME, Submission centre name 
# 7.  SUBMISSION_ID, SRA/ERA submission accession 
# 8.  SUBMISSION_DATE, Date sequence submitted, YYYY-MM-DAY
# 9.  SAMPLE_ID, SRA/ERA sample accession 
# 10. SAMPLE_NAME, Sample name
# 11. POPULATION, Sample population
# 12. EXPERIMENT_ID, Experiment accession 
# 13. INSTRUMENT_PLATFORM, Type of sequencing machine
# 14. INSTRUMENT_MODEL, Model of sequencing machine 
# 15. LIBRARY_NAME, Library name  
# 16. RUN_NAME, Name of machine run
# 17. RUN_BLOCK_NAME, Name of machine run sector  
# 18. INSERT_SIZE, Submitter specifed insert size 
# 19. LIBRARY_LAYOUT, Library layout, this can be either PAIRED or SINGLE 
# 20. PAIRED_FASTQ, Name of mate pair file if exists (Runs with failed mates will have 
#     a library layout of PAIRED but no paired fastq file)
# 21. WITHDRAWN, 0/1 to indicate if the file has been withdrawn, only present if a file has been withdrawn
# 22. WITHDRAWN_DATE, date of withdrawal, this should only be defined if a file is 
#     withdrawn
# 23. COMMENT, comment about reason for withdrawal
# 24. READ_COUNT, read  count for the file
# 25. BASE_COUNT, basepair count for the file
def countBases(samples, seqIndex):
    total = 0
    
    for project, sampleID, withdrawnP, bases in flatFileIterator(seqIndex, [3,9,20,24]):
        if ( withdrawnP == "0" and useProject(project) and sampleID in samples ):
            if OPTIONS.verbose: print project, sampleID, withdrawnP, bases
            sample = samples[sampleID]
            sample.rawBases += int(bases)
            total += int(bases)
    
    printStatus(samples)
    print 'Total raw bases', total
    return total

def printStatus(samples):
    for sample in samples.itervalues():
        print sample

def findVariantEvalResults(key, file, type=str):
    def capture1(line):
        if key in line:
            s = line.split()
            return type(s[len(s)-1])
        else:
            return None

    return [val for val in map(capture1, open(file)) if val != None]


def getDBSNPRate(file):
    if file != None:
        key = "[evaluation_name=eval].[comparison_name=dbsnp].[jexl_expression=none].[filter_name=called].[novelty_name=all].[analysis=Comp Overlap].[data_point=% evals at comp]"
        return findVariantEvalResults(key, file, float)[0]
    else:
        return -1

def useProject(project):
    return OPTIONS.project == None or project == OPTIONS.project

def countMappedBases(samples, alignmentIndex):
    if ( OPTIONS.coverageFile != None ): 
        # read from summary file, looking for the line:
        # Total   340710  1187.14 N/A     N/A     N/A
        for parts in map( string.split, open(OPTIONS.coverageFile) ):
            if parts[0] == "Total":
                return -1, int(parts[1])
    else:
        return readMappedBasesFromBAS(samples, alignmentIndex)

def readMappedBasesFromBAS(samples, alignmentIndex):
    totalBases = 0 
    totalMapped = 0
    
    for project, sampleID, basFile in flatFileIterator(alignmentIndex, [2,3,6]):
        #print project, sampleID, basFile
        if ( useProject(project) and sampleID in samples ):
            if OPTIONS.verbose: print project, sampleID, basFile
            sample = samples[sampleID]
    
            for rawBases, mappedBases in flatFileIterator(os.path.join(OPTIONS.root, basFile), [7, 8], skip=1):
                #print '  ->', rawBases, mappedBases
                if OPTIONS.rawBasesFromBas:
                    sample.rawBases += int(rawBases)
                    totalBases += int(rawBases)
                sample.mappedBases += int(mappedBases)
                totalMapped += int(mappedBases)
                #print '  totals', totalBases, totalMapped
    
    printStatus(samples)
    print 'Total raw    bases', totalBases
    print 'Total mapped bases', totalMapped
    
    return totalBases, totalMapped

def countSNPs(samples, snpsVCF, useIndels = False):
    total = 0
    
    header, columnNames, remainingLines = vcfReader.readVCFHeader(open(snpsVCF))
    sampleIDs = columnNames[9:]

    for header, vcf, counter in vcfReader.lines2VCF(remainingLines, extendedOutput = True, decodeAll = False):
        if ( counter > OPTIONS.maxRecords and OPTIONS.maxRecords != -1 ):
            break
        if vcf.passesFilters():
            if ( vcf.isVariant() ): 
                total += 1
                
                if ( OPTIONS.verbose and total % 10000 == 0 ):
                    print '  progress', vcf.getChrom(), vcf.getPos()
                
                genotypes = vcf.rest[1:]
                for sampleID, genotypeField in itertools.izip(sampleIDs, genotypes):
                    #print sampleID, samples
                    if sampleID in samples:
                        genotype = genotypeField.split(':')[0]
                        variant = genotype != "0/0" and genotype != "0|0" and genotype != "0\0" and genotype != "./."
                        #print '  => ', vcf, sampleID, genotype, variant
                        if variant:
                            if ( useIndels ):
                                samples[sampleID].nIndels += 1
                            else:
                                samples[sampleID].nSNPs += 1

    printStatus(samples)
    return total

def countIndels(samples, indelsVCF):
    total = 0
    
    if ( indelsVCF != None ):
        return countSNPs(samples, indelsVCF, True)
    
    return total

def readSamples(vcf):
    print 'Reading samples for', OPTIONS.population
    header, columnNames, remainingLines = vcfReader.readVCFHeader(open(vcf))
    samples = map(Sample, columnNames[9:])
    if ( OPTIONS.onlySample != None ):
        samples = filter( lambda x: x.getName() == OPTIONS.onlySample, samples )
    
    print 'No. samples: ', len(samples)
    print 'Samples:     ', map(Sample.getName, samples)

    return dict(map( lambda x: (x.getName(), x), samples))
    
if __name__ == "__main__":
    usage = "usage: %prog"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--alignmentIndex", dest="alignmentIndex",type='string', default=None, help="1KG formated alignment index file")
    parser.add_option("-s", "--sequenceIndex", dest="sequenceIndex", type='string', default=None, help="1KG formated sequence index file")
    parser.add_option("", "--onlySample", dest="onlySample", type='string', default=None, help="If provide, only this sample will be processed")
    parser.add_option("", "--snps", dest="snps", type='string', default=None, help="SNPs VCF")
    parser.add_option("", "--snpsEval", dest="snpsVE", type='string', default=None, help="SNPs VCF VariantEval")
    parser.add_option("", "--indels", dest="indels", type='string', default=None, help="Indels VCF")
    parser.add_option("", "--indelsEval", dest="indelsVE", type='string', default=None, help="Indels VCF VariantEval")
    parser.add_option("", "--totalGenome", dest="totalGenome", type='float', default=2.96e9, help="Size, in bp, of the callable genome")
    parser.add_option("", "--calledGenome", dest="calledGenome", type='float', default=None, help="Size, in bp, of the callable genome")
    parser.add_option("-p", "--pop", dest="population", type='string', default="Anonymous", help="Population")
    parser.add_option("", "--project", dest="project", type='string', default=None, help="If provided, will only include fastq/BAM files that match this project in the stats calculations")
    parser.add_option("-r", "--root", dest="root",type='string', default=".", help="Path to the 1KG data")
    parser.add_option("-M", "--maxRecords", dest="maxRecords", type='int', default=-1, help="If provided, will only include fastq/BAM files that match this regex in the stats calculations")
    parser.add_option("-v", "--verbose", dest="verbose", action='store_true', default=False, help="If provided, will be verbose during output")
    parser.add_option("", "--rawBasesFromBas", dest="rawBasesFromBas", action='store_true', default=False, help="If provided, we'll take our raw base counts from the BAS file")
    parser.add_option("-o", "--output", dest="output",type='string', default=None, help="Path to the 1KG data")
    parser.add_option("-c", "--coverageFile", dest="coverageFile",type='string', default=None, help="Path to GATK DoC .sample_summary file")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")

    samples = readSamples(OPTIONS.snps)
    nSamples = len(samples)    

    ignore, totalMappedBases = countMappedBases(samples, OPTIONS.alignmentIndex)    
    totalBases = countBases(samples, OPTIONS.sequenceIndex)
    meanMappedDepth = totalMappedBases / OPTIONS.totalGenome / nSamples
    totalSNPs = countSNPs(samples, OPTIONS.snps)
    totalIndels = countIndels(samples, OPTIONS.indels)
    
    snpNoveltyRate = 100 - getDBSNPRate(OPTIONS.snpsVE)
    indelNoveltyRate = 100 - getDBSNPRate(OPTIONS.indelsVE)

    out = sys.stdout
    if ( OPTIONS.output != None ): out = open(OPTIONS.output, 'w')
    print >> out, 'number of samples', nSamples
    print >> out, 'total raw bases', totalBases
    print >> out, 'total mapped bases', totalMappedBases

    # mean mapped depth is total bases mapped divided by acgt reference base count divided by number of individuals, after rmdup: for exons this is calculated on the target region only
    
    print >> out, 'mean mapped depth', meanMappedDepth
    
    print >> out, 'bases called (fraction ref genome) %f (%.2f%%)' % (OPTIONS.calledGenome, 100.0 * OPTIONS.calledGenome / OPTIONS.totalGenome)
    print >> out, 'number of SNP sites (%% novel) %d (%.2f%%)' % (totalSNPs, snpNoveltyRate)
    print >> out, 'average # SNP sites per individual', average(map(Sample.getNSNPs, samples.itervalues()))
    print >> out, 'number of indel sites (%% novel) %d (%.2f%%)' % (totalIndels, indelNoveltyRate)
    print >> out, 'average # indel sites per individual', average(map(Sample.getNIndels, samples.itervalues()))
    out.close()
    
