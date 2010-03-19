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

IS_OFFSET = True
CHR_OFFSET = 5
START_OFFSET = 6
END_OFFSET = 7
TYPE_OFFSET = 11

class RepeatInfo:
    def __init__(self, name):
        self.name = name
        self.count = 0
        self.coveredBases = 0

def badChr(excludes, chr):
    if excludes == None:
        return False

    for exclude in excludes:
        if chr.find(exclude) != -1:
            return True
    return False

def main():
    global OPTIONS
    usage = "usage: %prog stage [options]"
    parser = OptionParser(usage=usage)
#     parser.add_option("-q", "--farm", dest="farmQueue",
#                         type="string", default=None,
#                         help="Farm queue to send processing jobs to")
    parser.add_option("", "--header", dest="header",
                        type='string', default=None,
                        help="interval_list file")
    parser.add_option("-o", "--output", dest="output",
                        type='string', default=None,
                        help="output interval_list filename")
    parser.add_option("-r", "--ref", dest="ref",
                        type='string', default=None,
                        help="reference name -- either hg18 or b36")
    parser.add_option("-x", "--exclude", dest="excludes",
                        action="append", type='string',
                        help="If provided, only run pipeline for this sample")
    parser.add_option("-m", "--maxRecords", dest="maxRecords",
                        type='int', default=None,
                        help="If provided, max. number of records to process")
    parser.add_option("-s", "--skip", dest="skip",
                        type='int', default=None,
                        help="If provided, only process every skip records")
    parser.add_option("", "--excludeChr", dest="excludeChrs",
                        action="append", type='string',
                        help="If provided, don't include chr matching this string")
                       
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    repeatFile = args[0]

    # open the output file
    out = open(OPTIONS.output, 'w')

    # write out the header
    for line in open(OPTIONS.header):
        out.write(line)

    info = dict()
    nRecords = 1
    i = 1
    for line in open(repeatFile):
        if OPTIONS.maxRecords != None and nRecords > OPTIONS.maxRecords:
            break
    
        if len(line) == 0 or line[0] == '#':
            continue
    
        parts = line.split()
        type = parts[TYPE_OFFSET]
        start = int(parts[START_OFFSET]) + 1
        end = int(parts[END_OFFSET]) + 1
        chr = parts[CHR_OFFSET]

        if OPTIONS.ref == 'b36':
            chr = chr.replace('chrM', 'chrMT')
            chr = chr.replace('chr', '')

        if (OPTIONS.excludes == None or type not in OPTIONS.excludes) and not badChr(OPTIONS.excludeChrs, chr):
            name = 'repeat_' + str(i) + '_' + type
            strand = '+'

            if OPTIONS.skip == None or i % OPTIONS.skip == 0:
                print >> out, '\t'.join(map(str, [chr, start, end, strand, name]))
                nRecords += 1

            i += 1

            if type not in info:
                info[type] = RepeatInfo(type)
    
            typeInfo = info[type]
            typeInfo.count += 1
            typeInfo.coveredBases += end - start

    out.close()

    for typeInfo in info.values():
        print "%20s\t%20d\t%20d" % ( typeInfo.name, typeInfo.count, typeInfo.coveredBases )
        
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
