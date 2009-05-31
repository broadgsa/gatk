import farm_commands
import os.path
import sys
from optparse import OptionParser
import string
import re
import glob
import picard_utils
import itertools

gatkPath = "~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar"
ref = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
analysis = "CombineDuplicates"
   
def main():    
    global OPTIONS, ROOT

    usage = "usage: %prog lanes.list nIndividuals [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-q", "--farmQueue", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to submit jobs to.  Leave blank for local processing")
    parser.add_option("-l", "--lod", dest="lod",
                        type="float", default=5,
                        help="minimum lod for calling a variant")
    parser.add_option("-k", "--column", dest="column",
                        type="int", default=1,
                        help="Column in the file with the geli file path")
    parser.add_option("-o", "--output", dest="output",
                        type="string", default='/dev/stdout',
                        help="x")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    lines = [line.split() for line in open(args[0])]
    nIndividuals = int(args[1])
    gelis = map( lambda x: x[OPTIONS.column-1], lines )
    variantsOut = map( lambda geli: os.path.split(geli)[1] + '.calls', gelis)

    print gelis
    print variantsOut

    nTotalSnps = 0
    nNovelSnps = 0
    for geli in gelis:
        root, flowcellDotlane, ext = picard_utils.splitPath(geli)
        dbsnp_matches = os.path.join(root, flowcellDotlane) + '.dbsnp_matches'
        TOTAL_SNPS, NOVEL_SNPS, PCT_DBSNP, NUM_IN_DB_SNP = picard_utils.read_dbsnp(dbsnp_matches)
        nTotalSnps += int(TOTAL_SNPS)
        nNovelSnps += int(NOVEL_SNPS)
        print 'DATA:    ', flowcellDotlane, TOTAL_SNPS, NOVEL_SNPS, PCT_DBSNP, NUM_IN_DB_SNP, dbsnp_matches
    print 'DATA:    TOTAL SNP CALLS SUMMED ACROSS LANES, NOT ACCOUNT FOR IDENTITY', nTotalSnps
    print 'DATA:    NOVEL SNP CALLS SUMMED ACROSS LANES, NOT ACCOUNT FOR IDENTITY ', nNovelSnps
    print 'DATA:    AVERAGE DBSNP RATE ACROSS LANES ', float(nTotalSnps - nNovelSnps) / nTotalSnps

    jobid = None
    for geli, variantOut in zip(gelis, variantsOut):
        if not os.path.exists(variantOut):
            cmd = ("GeliToText.jar I=%s | awk '$7 > %f' > %s" % ( geli, OPTIONS.lod, variantsOut) )
            #jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, just_print_commands=False)

    cmd = ("cat %s | awk '$1 !~ \"@\" && $1 !~ \"#Sequence\" && $0 !~ \"GeliToText\"' | sort -k 1 -k 2 -n > tmp.calls" % ( ' '.join(variantsOut) ) )
    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, just_print_commands=False, waitID = jobid)

    sortedCallFile = 'all.sorted.calls'
    cmd = ("~/dev/GenomeAnalysisTK/trunk/perl/sortByRef.pl -k 1 tmp.calls ~/work/humanref/Homo_sapiens_assembly18.fasta.fai > %s"  % ( sortedCallFile ) )    
    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, just_print_commands=False, waitID = jobid)
     
    sortedCalls = [line.split() for line in open(sortedCallFile)]
    aggregratedCalls = picard_utils.aggregateGeliCalls(sortedCalls)
    
    outputFile = open(OPTIONS.output, 'w')
    print >> outputFile, 'loc ref alt EM_alt_freq discovery_likelihood discovery_null discovery_prior discovery_lod EM_N n_ref n_het n_hom'

    for snp in map( lambda x: picard_utils.aggregatedGeliCalls2SNP(x, nIndividuals), aggregratedCalls ):
        if snp == None: continue    # ignore bad calls
        #print snp
        #sharedCalls = list(sharedCallsGroup)
        #genotype = list(sharedCalls[0][5])
        print >> outputFile, '%s   %s %s %.6f -420.0 -420.0 0.000000 100.0 %d %d %d %d' % (snp.loc, snp.ref, snp.alt(), snp.q(), nIndividuals, snp.nRefGenotypes(), snp.nHetGenotypes(), snp.nHomVarGenotypes())
    outputFile.close()
    
                
if __name__ == "__main__":
    main()
