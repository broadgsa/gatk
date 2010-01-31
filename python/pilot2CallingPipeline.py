import farm_commands
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import faiReader
import math
import shutil

GATK_STABLE = 'java -ea -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar -l INFO -R /humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta '
GATK_DEV = 'java -ea -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -l INFO -R /humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta '
GATK = GATK_STABLE

class CallTarget:
    def __init__(self, name, hetero, minQ = 50, depth = 120, truthGFFName = None, variantEvalArgs = '', otherVCF = None):
        self.name = name
        self.hetero = hetero
        self.minQ = minQ
        self.depth = depth
        if truthGFFName <> None:
            self.truthGFFName = truthGFFName
        else:
            self.truthGFFName = name
        self.variantEvalArgs = variantEvalArgs
        self.otherVCF = otherVCF
        self.unionVCF = None
        if otherVCF != None:
            self.unionVCF = self.name + '.gatk.glftrio.union.filtered.vcf'
        self.listFile = None
        self.releaseVCF = None

CEU_HET = 0.79e-3
YRI_HET = 1.0e-3

targets = [
    CallTarget('NA12878', CEU_HET),
    CallTarget('NA12891', CEU_HET),
    CallTarget('NA12892', CEU_HET),
    CallTarget('NA19238', YRI_HET),
    CallTarget('NA19239', YRI_HET),
    CallTarget('NA19240', YRI_HET),

    CallTarget('ceu.trio', CEU_HET, 50, 360, 'NA12878', '--sampleName NA12878', '/humgen/gsa-hpprojects/1kg/1kg_pilot2/currentBestProjectCalls/CEU_1kg_pilot2.vcf'),
    CallTarget('yri.trio', YRI_HET, 50, 360, 'NA19240', '--sampleName NA19240', '/humgen/gsa-hpprojects/1kg/1kg_pilot2/currentBestProjectCalls/YRI_1kg_pilot2.vcf'),

#    CallTarget('NA19240.alltechs.solid_original', YRI_HET, 50, 120, 'NA19240'),
#    CallTarget('NA12878.alltechs.solid_original', CEU_HET, 50, 120, 'NA12878'),

    CallTarget('NA19240.alltechs.solid_recal', YRI_HET, 50, 120, 'NA19240'),
    CallTarget('NA12878.alltechs.solid_recal', CEU_HET, 50, 120, 'NA12878'),
]

sets = ['Intersection', 'filteredInBoth', 'gatk', 'gatk-filteredInOther', 'glftrio', 'glftrio-filteredInOther', 'Intersection', ['gatk-unique', 'gatk.*'], ['glftrio-unique', 'glftrio.*']]

def main():
    global OPTIONS
    usage = "usage: %prog stage [options]"
    parser = OptionParser(usage=usage)
#     parser.add_option("-q", "--farm", dest="farmQueue",
#                         type="string", default=None,
#                         help="Farm queue to send processing jobs to")
    parser.add_option("", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, nothing actually gets run, just a dry run")
    parser.add_option("-s", "--sample", dest="useSample",
                        type='string', default=None,
                        help="If provided, only run pipeline for this sample")
    parser.add_option("-c", "--minQ", dest="minQ",
                        type='float', default=None,
                        help="If provided, will actually use this Q threshold for calls")
    parser.add_option("-d", "--dir", dest="dir",
                        type='string', default="",
                        help="If provided, this is the root where files are read and written")
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
                       
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    stages = args[0].split(",")

    allJobs = []
    for callTarget in targets:
        if callTarget.name == OPTIONS.useSample or OPTIONS.useSample == None:
            lastJobs = None
            
            if OPTIONS.minQ != None:
                callTarget.minQ = OPTIONS.minQ
            
            for stage in stages:
                print 'STAGE', stage
                target = callTarget.name
                callTarget.listFile = os.path.join("lists", target + ".list")
                unfilteredVCFBaseName = target + '.gatk.ug.vcf'
                unfilteredVCF = os.path.join(OPTIONS.dir, unfilteredVCFBaseName)
                filteredVCF = os.path.join(OPTIONS.dir, target + '.gatk.ug.filtered.vcf')
    
                tmpdir = os.path.join(OPTIONS.dir, "intermediates", unfilteredVCFBaseName + ".scatter")
                if not os.path.exists(tmpdir):
                    os.makedirs(tmpdir)
    
                if callTarget.unionVCF != None:
                    callTarget.releaseVCF = os.path.join(OPTIONS.dir, callTarget.name + ".gatk_glftrio.intersection.annotated.filtered.vcf")
    
                print 'Heading into stage', stage, lastJobs
    
                newJobs = []
                if stage == 'CALL':
                    newJobs = callSNPs(target, callTarget, callTarget.listFile, tmpdir)
        
                if stage == 'MERGE':
                    newJobs = mergeSNPs(target, lastJobs, unfilteredVCF, tmpdir)
        
                if stage == 'FILTER':
                    newJobs = filterSNPs(callTarget, lastJobs, unfilteredVCF, filteredVCF)
        
                if stage == 'UNION':
                    newJobs = unionSNPs(callTarget, lastJobs, filteredVCF )
    
                if stage == 'SELECT_INTERSECT':
                    if ( callTarget.unionVCF != None ):
                        # we are one fo the release targets
                        newJobs = finalize1KGRelease(callTarget, lastJobs, os.path.join(OPTIONS.dir, callTarget.unionVCF ), callTarget.releaseVCF)
    
                if stage == 'EVAL':
                    newJobs = evalSNPs(callTarget, lastJobs, unfilteredVCF, filteredVCF)
    
                if stage == 'RELEASE':
                    dir = '/humgen/gsa-scr1/pub/1000GenomesPilot2'
                    subdir = '1000GenomesPilot2SNPs_GATK_glftrio_' + date.today().strftime("%m%d%y")
                    releaseSNPs(os.path.join(dir, subdir), callTarget, filteredVCF )
        
                if stage == 'CLEAN':
                    shutil.rmtree("intermediates", True)
                    for file in [filteredVCF, unfilteredVCF]:
                        if os.path.exists(file): os.remove(file)
        
                print 'New jobs'
                for job in newJobs:
                    print '    ', job
                allJobs.append(newJobs)
                if newJobs != []: 
                    lastJobs = newJobs

    print 'EXECUTING JOBS'
    farm_commands.executeJobs(allJobs, farm_queue = OPTIONS.farmQueue, just_print_commands = OPTIONS.dry) 

hg18 = ['chr' + str(i) for i in range(1,23)] + ['chrX']
b36 = [str(i) for i in range(1,23)] + ['X']

def autosomePlusX(name):
    return name in b36 or name in hg18

def partition(faiRecords, maybeChrom, n):
# 1       247249719       3       60      61
# 2       242951149       251370554       60      61
    chromAndSize = [[x[0], int(x[1])] for x in faiRecords if autosomePlusX(x[0])]
    #print chromAndSize
    if maybeChrom <> None:
        chromAndSize = filter(lambda x: x[0] == maybeChrom, chromAndSize)
    
    def r(chrom, chrStart, chrEnd):
        outputVCF = 'calls.chr%s.good.%s.vcf' % (chrom, chrStart)
        L = '-L %s:%d-%d' % (chrom, chrStart, chrEnd)
        return outputVCF, chrom, chrStart, chrEnd, L
    
    for chrom, size in chromAndSize:
        #print 'SIZE', chrom, size
        if n == 1:
            yield r(chrom, 1, size)        
        else:
            idealBp = float(size-1) / n
            chrEnd = None
            chrStart = 1
            for i in range(1, n):
                chrEnd = 1 + int(round(idealBp * (i + 1)))
                #print chrom, chrStart, chrEnd, idealBp
                result = r(chrom, chrStart, chrEnd)
                chrStart = chrEnd + 1
                yield result
            if chrEnd <> size:
                raise Exception('X')

def callSNPs(target, callTarget, listFile, tmpdir):
    extras = "-mmq 10 -mbq 10 -pl SOLID --heterozygosity %e" % (callTarget.hetero)
    fai = "/humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta.fai"
    NWaysParallel = 5

    bins = partition(faiReader.readFAI(fai), None, NWaysParallel)
    jobs = list()
    for outputVCF, chrom, chrStart, chrEnd, L in bins:
        outputVCF = os.path.join(tmpdir, outputVCF)
        #print outputVCF, L, os.path.exists(outputVCF)
        #if not os.path.exists(outputVCF):
        print 'Enqueuing job for', outputVCF, L
        cmd = 'java -Xmx2048m -jar /home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar ' + \
            '-T UnifiedGenotyper -R /humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -mrl 500000 ' + \
            '-I %s -confidence %d %s -varout %s -vf VCF -L %s -gm JOINT_ESTIMATE' % (listFile, callTarget.minQ, extras, outputVCF, L)
        #jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry)
        jobs.append(farm_commands.FarmJob(cmd, jobName = "CALL_%s_c%s_s%s" % (target, str(chrom), str(chrStart))))

    return jobs

def mergeSNPs(target, lastJobs, snpFile, tmpdir):
    #print 'LastJobs = ', lastJobs
    cmd = "python ~/dev/GenomeAnalysisTK/trunk/python/mergeVCFs.py -a -f /humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta.fai %s/*.vcf > %s" % (tmpdir, snpFile)
    return [farm_commands.FarmJob(cmd, jobName = "MERGE_%s" % (target), dependencies = lastJobs)]

def filterSNPs(callTarget, lastJobs, unfilteredVCF, filteredVCF):
    target = callTarget.name
    #expression = "AB > 0.75 || DP > %s" % depth
    expression1 = ['GATK_STANDARD', "AB > 0.75 || DP > %s || MQ0 > 40 || SB > -0.10" % callTarget.depth]
    expression2 = ['HARD_TO_VALIDATE', "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"]
    cmd = GATK + '-T VariantFiltration -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -B variant,VCF,%s --clusterWindowSize 10 -o %s ' % (unfilteredVCF, filteredVCF)
    
    for name, exp in [expression1, expression2]:
        cmd += '--filterName %s --filterExpression "%s" ' % ( name, exp )
    
    #jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry)
    return [farm_commands.FarmJob(cmd, jobName = "FILTER_%s" % (target), dependencies = lastJobs)]

def evalSNPs(callTarget, lastJobs, unfilteredVCF, filteredVCF):
    evalRoot = os.path.join(OPTIONS.dir, "eval")
    if not os.path.exists(evalRoot):
        os.makedirs(evalRoot)

    def eval1(vcfFullPath, namePostfix = "", args = ""):
        vcf = os.path.split(vcfFullPath)[1]
        out = os.path.join(OPTIONS.dir, "eval", vcf + namePostfix + ".eval")
        cmd = GATK + "-T VariantEval -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -B 1kg_ceu,VCF,/humgen/gsa-hpprojects/1kg/1kg_pilot1/SNPCalls/Joint/RC1/CEU.2and3_way.annotated.vcf -B 1kg_yri,VCF,/humgen/gsa-hpprojects/1kg/1kg_pilot1/SNPCalls/Joint/RC1/YRI.2and3_way.annotated.vcf -B eval,VCF,%s -G -A -o %s -L %s" % ( vcfFullPath, out, '\;'.join(map(str, b36)) )
        hapmap3 = os.path.join("../../hapmap3Genotypes", callTarget.truthGFFName + ".b36.gff") 
        if os.path.exists(hapmap3):
            cmd += " -B hapmap-chip,GFF,%s %s %s" % (hapmap3, callTarget.variantEvalArgs, args)
        #jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry)
        return farm_commands.FarmJob(cmd, jobName = "EVAL_%s_%s" % (callTarget.name, namePostfix), dependencies = lastJobs)
    
    jobs = []
    jobs.append(eval1(filteredVCF))
    jobs.append(eval1(unfilteredVCF))
    if callTarget.unionVCF != None:
        jobs.append(eval1(os.path.join(OPTIONS.dir, callTarget.unionVCF), ".union"))
        for set in sets:
            if type(set) == list:
                name, selector = set
            else:
                name, selector = set, set
            jobs.append(eval1(os.path.join(OPTIONS.dir, callTarget.unionVCF), "." + name, " -vcfInfoSelector set=\"" + selector + "\""))

    if callTarget.releaseVCF != None:
        jobs.append(eval1(callTarget.releaseVCF, ".release"))


    return jobs

def unionSNPs(callTarget, lastJobs, filteredVCF ):
    if callTarget.otherVCF != None:
        cmd = GATK + "-T VCFCombine -B GATK,VCF,%s -B glfTrio,VCF,%s -O %s -type UNION -priority GATK,glfTrio -A" % ( filteredVCF, callTarget.otherVCF, os.path.join(OPTIONS.dir, callTarget.unionVCF) ) 
        return [farm_commands.FarmJob(cmd, jobName = "UNION_%s" % (callTarget.name), dependencies = lastJobs)]
    else:
        return []
#        jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry)

def finalize1KGRelease(callTarget, lastJobs, inputVCF, outputVCF):
    commands = []
    lastJob = lastJobs
    
#     subdir = os.path.join(OPTIONS.dir, "1kg_release")
#     if not os.path.exists(subdir):
#         os.makedirs(subdir)

    # first, filter out the intersection
    cmd = GATK + '-T VCFSelect -B variant,VCF,%s -o %s -match "(set eq \'Intersection\' || set eq \'filteredInBoth\')" ' % (inputVCF, outputVCF)
    lastJob = farm_commands.FarmJob(cmd, jobName = "SELECT_%s" % (callTarget.name), dependencies = lastJob)
    commands.append(lastJob)
    
    # call variant annotator to annotate all of the calls
    # annotatedVCF = os.path.join(subdir, callTarget.name + ".gatk_glftrio.intersection.annotated.vcf")
#     cmd = GATK + '-T VariantAnnotator -I %s -standard -B variant,VCF,%s -vcf %s -L 1:1-1,000,000 ' % (callTarget.listFile, intersectVCF, annotatedVCF)
#     lastJob = farm_commands.FarmJob(cmd, jobName = "ANNOTATE_%s" % (callTarget.name), dependencies = lastJob)
#     commands.append(lastJob)
#     
#     # filter the calls
#     filteredVCF = os.path.join(subdir, callTarget.name + ".gatk_glftrio.intersection.annotated.filtered.vcf")
#     commands = commands + filterSNPs(callTarget, lastJob, annotatedVCF, filteredVCF)
#     
    return commands

def releaseSNPs(dir, callTarget, filteredVCF ):
    if not os.path.exists(dir):
        os.makedirs(dir)

    print 'Copying files into ', dir, filteredVCF, callTarget.unionVCF, callTarget.releaseVCF
    if not OPTIONS.dry: shutil.copy(filteredVCF, dir)
    if callTarget.unionVCF != None:
        if not OPTIONS.dry: shutil.copy(os.path.join(OPTIONS.dir, callTarget.unionVCF), dir)
    if callTarget.releaseVCF != None:
        if not OPTIONS.dry: shutil.copy(callTarget.releaseVCF, dir)

if __name__ == "__main__":
    main()


# java -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -T VCFCombine -R /humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta -B GATK,VCF,ceu.trio.gatk.ug.filtered.vcf -B glfTrio,VCF,/humgen/gsa-hpprojects/1kg/1kg_pilot2/currentBestProjectCalls/CEU_1kg_pilot2.vcf -O test.vcf -type UNION -priority GATK,glfTrio -l INFO  -A
# java -ea -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -l INFO -R /humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta -T VariantEval -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -B eval,VCF,test.vcf -B hapmap-chip,GFF,../../hapmap3Genotypes/NA12878.b36.gff --sampleName NA12878 -vcfInfoSelector set=gatk-filtered


# if ( $1 == 3.1 ) then
# cat ceu.trio.calls.allTechs.mmq10_mbq10_q200.filtered.vcf | cut -f 1-10 > ceu.trio.calls.allTechs.mmq10_mbq10_q200.NA12878only.filtered.vcf
# endif
# 
# if ( $1 == 4 ) then
# foreach callset ( NA12878.allTechs.mmq10_mbq10_q200 ceu.trio.calls.allTechs.mmq10_mbq10_q200.NA12878only )
# java -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar -T CallsetConcordance -R /humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta -B GATK,VCF,$callset.filtered.vcf -B glfTrio,VCF,CEU_1kg_pilot2.na12878.vcf -CT SimpleVenn -CO ${callset}_v_CEU_1kg_pilot2.filtered.vcf -l INFO
# cat ${callset}_v_CEU_1kg_pilot2.filtered.vcf | awk '$1 ~ "#" || $8 ~ "callset2_only"' > ${callset}_v_CEU_1kg_pilot2.filtered.CEU_1kg_pilot2Unique.vcf
# cat ${callset}_v_CEU_1kg_pilot2.filtered.vcf | awk '$1 ~ "#" || $8 ~ "callset1_only"' > ${callset}_v_CEU_1kg_pilot2.filtered.${callset}Unique.vcf
# cat ${callset}_v_CEU_1kg_pilot2.filtered.vcf | awk '$1 ~ "#" || $8 ~ "concordant"' > ${callset}_v_CEU_1kg_pilot2.filtered.concordant.vcf
# end
# endif
# 
# if ( $1 == 5 ) then
# mkdir /humgen/gsa-scr1/pub/1000Pilot2_010710
# 
# foreach file ( NA12878.allTechs.mmq10_mbq10_q200.filtered.vcf NA12878.SLX.mmq10_mbq10_q50.filtered.vcf NA12891.calls.mmq10_mbq10_q50.filtered.vcf NA12892.calls.mmq10_mbq10_q50.filtered.vcf ceu.trio.calls.allTechs.mmq10_mbq10_q200.filtered.vcf )
# echo $file
# cp $file /humgen/gsa-scr1/pub/1000Pilot2_010710
# end
# 
# endif
