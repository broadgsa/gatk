#!/usr/bin/python

import sys, os
from farm_commands import cmd

#if __name__ == "__main__":
from optparse import OptionParser
class OptionParser (OptionParser):

    def check_required (self, opts):
        for opt in opts:
            option = self.get_option(opt)

            # Assumes the option's 'default' is set to None!
            if getattr(self.values, option.dest) is None:
                self.error("%s option not supplied" % option)

def exists_and_non_zero(file):
    return os.path.exists(file) and os.path.getsize(file) != 0

# Global files
gatk_exe = "java -Xmx8192m -jar ~/dev/Sting/trunk/dist/GenomeAnalysisTK.jar "
ref = " -R /broad/1KG/reference/human_b36_both.fasta"

parser = OptionParser()
parser.add_option("-b", "--bam-file", help="Input BAM filename", dest="bam_file", default=None)

(options, args) = parser.parse_args()
parser.check_required(("-b",))    
output_base = os.path.splitext(options.bam_file)[0]
print "Input BAM file:",options.bam_file

# Assure that hapmap chip file is here
hapmapchip_file = output_base+".hapmapchip.gff"
if not exists_and_non_zero(hapmapchip_file):
    sys.exit("Expected to find "+hapmapchip_file+" but it's not around")

# Run SSG calls
ssg_outfile = output_base+".SSG.calls"
if not exists_and_non_zero(ssg_outfile):
    print "Running SingleSampleGenotyper"
    Cmd = gatk_exe + ref + " -T SingleSampleGenotyper"+\
        " -I "+options.bam_file+\
        " -varout "+ssg_outfile+\
        " -L 1"+\
        " -fmq0 -l INFO "
    cmd(Cmd)

# Run VariantEval on SSG calls
varianteval_outdir = "VariantEval"
varianteval_outbase = "VariantEval/"+output_base
varianteval_outfile = "VariantEval/"+output_base+".all.genotype_concordance"

if not exists_and_non_zero(varianteval_outfile):
    print "Running VarianteEval with output to "+varianteval_outbase
    if not os.path.exists(varianteval_outdir):
        os.mkdir(varianteval_outdir)
    Cmd = gatk_exe + ref + "-T VariantEval"+\
        " -B dbsnp,dbsnp,/humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod"+\
        " -B eval,Variants,"+ssg_outfile+\
        " -hc "+hapmapchip_file+\
        " -o "+varianteval_outbase+\
        " -L /humgen/gsa-scr1/andrewk/coverage/b36.hapmap_intervals.chr1"+\
        " -evalContainsGenotypes"+\
        " -minDiscoveryQ 5"+\
        " -fmq0 -l INFO"
    cmd(Cmd)

# Run CoverageHist on BAM file
hist_outfile = output_base+".hapmap.coverage_hist"
if not exists_and_non_zero(hist_outfile):
    print "Running CoverageHist on "+options.bam_file
    Cmd = gatk_exe + ref + "-T CoverageHistogram"+\
        " -I "+options.bam_file+\
        " -o "+hist_outfile+\
        " -L /humgen/gsa-scr1/andrewk/coverage/b36.hapmap_intervals.chr1"+\
        " -fmq0 -l INFO"
        #" -hc ~/coverage/NA12878.b36.gff"+\
        #" -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod"+\ 
        
    cmd(Cmd)#,"gsa", output_base)

# Now show CoverageEval predictions
Cmd = "CoverageEval.py -e"+\
    " -s /humgen/gsa-scr1/projects/1kg_pilot2/coverage_eval/fixed_stats/NA12878.CLEANED.chr1.full.hapmap.stats"+\
    " -g "+hist_outfile
cmd(Cmd)
    

    
    
    
    
    
    
    
    
    