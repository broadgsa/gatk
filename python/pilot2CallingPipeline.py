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

CEU_HET = 0.79e-3
YRI_HET = 1.0e-3

targets = [
    CallTarget('NA12878', CEU_HET),
    CallTarget('NA12891', CEU_HET),
    CallTarget('NA12892', CEU_HET),
    CallTarget('NA19238', YRI_HET),
    CallTarget('NA19239', YRI_HET),
    CallTarget('NA19240', YRI_HET),

    CallTarget('NA19240.alltechs', YRI_HET, 50, 120, 'NA19240'),
    CallTarget('NA12878.alltechs', CEU_HET, 50, 120, 'NA12878'),

    CallTarget('ceu.trio', CEU_HET, 50, 360, 'NA12878', '--sampleName NA12878', '/humgen/gsa-hpprojects/1kg/1kg_pilot2/currentBestProjectCalls/CEU_1kg_pilot2.vcf'),
    CallTarget('yri.trio', YRI_HET, 50, 360, 'NA19240', '--sampleName NA19240', '/humgen/gsa-hpprojects/1kg/1kg_pilot2/currentBestProjectCalls/YRI_1kg_pilot2.vcf'),
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
    parser.add_option("-d", "--dir", dest="dir",
                        type='string', default="",
                        help="If provided, this is the root where files are read and written")
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
                       
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    stage = args[0]

    for callTarget in targets:
        target = callTarget.name
        listFile = os.path.join("lists", target + ".list")
        unfilteredVCF = os.path.join(OPTIONS.dir, target + '.gatk.ug.vcf')
        filteredVCF = os.path.join(OPTIONS.dir, target + '.gatk.ug.filtered.vcf')
        tmpdir = os.path.join(OPTIONS.dir, "intermediates", unfilteredVCF + ".scatter")
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)

        if stage == 'CALL':
            callSNPs(target, callTarget, listFile, tmpdir)

        if stage == 'MERGE':
            mergeSNPs(target, unfilteredVCF, tmpdir)

        if stage == 'FILTER':
            filterSNPs(target, callTarget.depth, unfilteredVCF, filteredVCF)

        if stage == 'UNION':
            unionSNPs(callTarget, filteredVCF )

        if stage == 'RELEASE':
            dir = '/humgen/gsa-scr1/pub/1000GenomesPilot2'
            subdir = '1000GenomesPilot2SNPs_GATK_glftrio_' + date.today().strftime("%m%d%y")
            releaseSNPs(os.path.join(dir, subdir), callTarget, filteredVCF )

        if stage == 'CLEAN':
            shutil.rmtree("intermediates", True)
            for file in [filteredVCF, unfilteredVCF]:
                if os.path.exists(file): os.remove(file)

        if stage == 'EVAL':
            evalSNPs(callTarget, unfilteredVCF, filteredVCF)

GATK_STABLE = 'java -ea -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar -l INFO -R /broad/1KG/reference/human_b36_both.fasta '
GATK_DEV = 'java -ea -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -l INFO -R /broad/1KG/reference/human_b36_both.fasta '
GATK = GATK_DEV

def callSNPs(target, callTarget, listFile, tmpdir):
    cmd = "python ./runmeCalls.py -q %s -N 5 -d %s -I %s /broad/1KG/reference/human_b36_both.fasta.fai -e \"-mmq 10 -mbq 10 -pl SOLID --heterozygosity %e\" -Q %s" % (OPTIONS.farmQueue, tmpdir, listFile, callTarget.hetero, callTarget.minQ)
    print 'Enqueuing job for', target, callTarget, listFile, tmpdir
    jobid = farm_commands.cmd(cmd, None, None, just_print_commands = OPTIONS.dry)

def mergeSNPs(target, snpFile, tmpdir):
    cmd = "python ~/dev/GenomeAnalysisTK/trunk/python/mergeVCFs.py -a -f /broad/1KG/reference/human_b36_both.fasta.fai %s/*.vcf > %s" % (tmpdir, snpFile)
    jobid = farm_commands.cmd(cmd, None, None, just_print_commands = OPTIONS.dry)

def filterSNPs(target, depth, unfilteredVCF, filteredVCF):
    #expression = "AB > 0.75 || DP > %s" % depth
    expression = "AB > 0.75 || DP > %s || MQ0 > 40 || SB > -0.10" % depth
    cmd = GATK + '-T VariantFiltration -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -B variant,VCF,%s --clusterWindowSize 10 -o %s --filterExpression "%s"' % (unfilteredVCF, filteredVCF, expression)
    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry)

def evalSNPs(callTarget, unfilteredVCF, filteredVCF):
    def eval1(vcf, namePostfix = "", args = ""):
        out = os.path.join(OPTIONS.dir, "eval", vcf + namePostfix + ".eval")
        cmd = GATK + "-T VariantEval -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -B 1kg_ceu,VCF,/humgen/gsa-hpprojects/1kg/1kg_pilot1/SNPCalls/Joint/RC1/CEU.2and3_way.annotated.vcf -B 1kg_yri,VCF,/humgen/gsa-hpprojects/1kg/1kg_pilot1/SNPCalls/Joint/RC1/YRI.2and3_way.annotated.vcf -B eval,VCF,%s -G -A -o %s -L %s" % ( vcf, out, '\;'.join(map(str, range(1,23))) )
        hapmap3 = os.path.join("../../hapmap3Genotypes", callTarget.truthGFFName + ".b36.gff") 
        if os.path.exists(hapmap3):
            cmd += " -B hapmap-chip,GFF,%s %s %s" % (hapmap3, callTarget.variantEvalArgs, args)
        jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry)
    eval1(filteredVCF)
    eval1(unfilteredVCF)
    if callTarget.unionVCF != None:
        eval1(callTarget.unionVCF, ".union")    
        for set in sets:
            if type(set) == list:
                name, selector = set
            else:
                name, selector = set, set
            eval1(callTarget.unionVCF, "." + name, " -vcfInfoSelector set=\"" + selector + "\"")    

def unionSNPs(target, filteredVCF ):
    if target.otherVCF != None:
        cmd = GATK + "-T VCFCombine -B GATK,VCF,%s -B glfTrio,VCF,%s -O %s -type UNION -priority GATK,glfTrio -A" % ( filteredVCF, target.otherVCF, target.unionVCF ) 
        jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry)

def releaseSNPs(dir, callTarget, filteredVCF ):
    if not os.path.exists(dir):
        os.makedirs(dir)

    print dir, filteredVCF
    #shutil.copy(filteredVCF, dir)
    if callTarget.unionVCF != None:
        shutil.copy(callTarget.unionVCF, dir)

if __name__ == "__main__":
    main()


# java -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -T VCFCombine -R /broad/1KG/reference/human_b36_both.fasta -B GATK,VCF,ceu.trio.gatk.ug.filtered.vcf -B glfTrio,VCF,/humgen/gsa-hpprojects/1kg/1kg_pilot2/currentBestProjectCalls/CEU_1kg_pilot2.vcf -O test.vcf -type UNION -priority GATK,glfTrio -l INFO  -A
# java -ea -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -l INFO -R /broad/1KG/reference/human_b36_both.fasta -T VariantEval -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -B eval,VCF,test.vcf -B hapmap-chip,GFF,../../hapmap3Genotypes/NA12878.b36.gff --sampleName NA12878 -vcfInfoSelector set=gatk-filtered


# if ( $1 == 3.1 ) then
# cat ceu.trio.calls.allTechs.mmq10_mbq10_q200.filtered.vcf | cut -f 1-10 > ceu.trio.calls.allTechs.mmq10_mbq10_q200.NA12878only.filtered.vcf
# endif
# 
# if ( $1 == 4 ) then
# foreach callset ( NA12878.allTechs.mmq10_mbq10_q200 ceu.trio.calls.allTechs.mmq10_mbq10_q200.NA12878only )
# java -Xmx4096m -jar /home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar -T CallsetConcordance -R /broad/1KG/reference/human_b36_both.fasta -B GATK,VCF,$callset.filtered.vcf -B glfTrio,VCF,CEU_1kg_pilot2.na12878.vcf -CT SimpleVenn -CO ${callset}_v_CEU_1kg_pilot2.filtered.vcf -l INFO
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
