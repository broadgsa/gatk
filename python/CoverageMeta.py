#!/usr/bin/env python

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
    
def process_bam(bam_file):
    gatk_exe = "java -Xmx8192m -jar ~/dev/Sting/trunk/dist/GenomeAnalysisTK.jar "
    output_basepath = os.path.splitext(bam_file)[0]
    output_basefile = os.path.basename(output_basepath)

    print "Input BAM file:",bam_file
    
    using_hapmapchip = None
    if False:
        # Config for chr1 of pilot 2 people
        ref = " -R /broad/1KG/reference/human_b36_both.fasta"
        #interval_string = "1"
    
        # Assure existance of intervals file
        intervals_file = output_basepath+".hapmap.intervals.chr1"
        if not exists_and_non_zero(intervals_file):
            Cmd = "awk -F' ' '{ if ($1 == \""+interval_string+"\") {print $1 \":\" $4} }' "+hapmapchip_file+" > "+intervals_file
            cmd(Cmd, die_on_fail = True)
    
        # Setup hapmapchip GFFs; assure that hapmap chip file is here
        using_hapmapchip = True
        hapmapchip_file = output_basepath+".hapmapchip.gff"
        if not exists_and_non_zero(hapmapchip_file):
            sys.exit("Expected to find "+hapmapchip_file+" but it's not around")
        
    else:
        # Config for pilot 3 
        ref = " -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
        intervals_file = "/seq/references/HybSelOligos/thousand_genomes_alpha_redesign/thousand_genomes_alpha_redesign.targets.interval_list"
    
    # Assure that BAM index exists
    bam_index_file = output_basepath+".bam.bai"
    if not exists_and_non_zero(bam_index_file):
        Cmd = "samtools index "+bam_file
        cmd(Cmd, die_on_fail = True)
        
    # Run SSG calls
    ssg_outfile = output_basepath+".SSG.calls"
    if not exists_and_non_zero(ssg_outfile):
        print "Running SingleSampleGenotyper"
        Cmd = gatk_exe + ref + " -T SingleSampleGenotyper"+\
            " -I "+bam_file+\
            " -varout "+ssg_outfile+\
            " -L "+intervals_file+\
            " -fmq0 -l INFO "
            #" --genotype -lod 5"
        
            #+interval_string+\
        cmd(Cmd, die_on_fail=True)
    
    # Run VariantEval on SSG calls
    varianteval_outdir = output_basepath+".VariantEval"
    varianteval_outbase = output_basepath+".VariantEval/"+output_basepath
    varianteval_outfile = output_basepath+".VariantEval/"+output_basepath #+".all.genotype_concordance"
    
    if not exists_and_non_zero(varianteval_outfile):
        print "Running VariantEval with output to "+varianteval_outbase
        if not os.path.exists(varianteval_outdir):
            os.mkdir(varianteval_outdir)
        Cmd = gatk_exe + ref + " -T VariantEval"+\
            " -B dbsnp,dbsnp,/humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod"+\
            " -B eval,Variants,"+ssg_outfile+\
            " -o "+varianteval_outbase+\
            " -L "+intervals_file+\
            " --evalContainsGenotypes"+\
            " -minConfidenceScore 5"+\
            " -fmq0 -l INFO"
        if using_hapmapchip:
            Cmd = Cmd + " -hc "+hapmapchip_file
    
        cmd(Cmd, die_on_fail=True)
    
    # Run CoverageHist on BAM file
    hist_outfile = output_basepath+".hapmap.coverage_hist"
    if not exists_and_non_zero(hist_outfile):
        print "Running CoverageHist on "+bam_file
        Cmd = gatk_exe + ref + " -T CoverageHistogram"+\
            " -I "+bam_file+\
            " -o "+hist_outfile+\
            " -L "+intervals_file+\
            " -fmq0 -l INFO"
            
        cmd(Cmd, die_on_fail=True)#,"gsa", output_basepath)
    
    # Now do CoverageEval predictions if they haven't been done
    oneline_stats_file = output_basepath+".oneline_stats"
    if not exists_and_non_zero(oneline_stats_file):
        Cmd = "CoverageEval.py -e"+\
            " -s /humgen/gsa-scr1/projects/1kg_pilot2/coverage_eval/fixed_stats/current/NA12878.chr1.full.stats"+\
            " -g "+hist_outfile+\
            " -v "+varianteval_outfile+\
            " -o "+oneline_stats_file
        cmd(Cmd, die_on_fail=True)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--bam-file", help="Input BAM filename", dest="bam_file", default=None)
    parser.add_option("-a", "--all_bams", help="Process all BAMs in current directory", default=False, dest="all_bams", action="store_true")
    
    (options, args) = parser.parse_args()
    
    if not options.all_bams:
        parser.check_required(("-b",))    
        process_bam(options.bam_file)
    else:
        bams = [f for f in os.listdir(".") if f.endswith(".bam")]
        print bams
        for bam in bams:
            process_bam(bam)
            #print

    
    
    
    
    
    
    
    
    