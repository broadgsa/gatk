import farm_commands
import os.path
import sys
from optparse import OptionParser
from gatkConfigParser import *
import glob
import itertools

if __name__ == "__main__":
    usage = """usage: %prog [-c config.cfg]*"""

    parser = OptionParser(usage=usage)
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
    parser.add_option("-c", "--config", dest="configs",
                        action="append", type="string", default=[],
                        help="Configuration file")                        
    parser.add_option("-w", "--wait", dest="initialWaitID",
                        type="string", default=None,
                        help="If providedm the first job dispatched to LSF will use this id as it ended() prerequisite")
    parser.add_option("", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, nothing actually gets run, just a dry run")
    parser.add_option("-i", "--ignoreExistingFiles", dest="ignoreExistingFiles",
                        action='store_true', default=False,
                        help="Ignores already written files, if present")
    parser.add_option("-d", "--dir", dest="outputdir",
                        type="string", default="./",
                        help="Output directory")

    (OPTIONS, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")

    config = gatkConfigParser(OPTIONS.configs)

    if not os.path.exists(OPTIONS.outputdir):
        os.mkdir(OPTIONS.outputdir)

    def outputDir(suffix):
        return os.path.join(OPTIONS.outputdir, suffix)
        
    intervals_dir = outputDir("intervals")
    cleaner_output = outputDir("cleaner")
    injector_output = outputDir("bams")
    snp_output = outputDir("calls/unfiltered_snps")
    filter_output = outputDir("calls/filtered_snps")
    indel_output = outputDir("calls/indels")

    samples = ["NA12878","NA12891","NA12892","NA19238","NA19239","NA19240"]
    DoCs = [76,88,66,56,68,86]
    chrs = range(1, 22) + ["X"]

    for sample, DoC in zip(samples, DoCs):
        #
        # Actually do some work here
        #
        myChrs = chrs
        if sample in ["NA12891", "NA19239"]:
            myChrs = chrs + ["Y"]

        for chr in myChrs:

            def filepath(tech):
                return "/broad/1KG/DCC/ftp/pilot_data/%s/alignment/%s.chrom%s.%s.SRP000032.2009_07.bam" % ( sample, sample, chr, tech )

            techs = ["SLX", "SOLID", "454"]
            MQs = [100,5,5]
            solexaBamFile = None
            solidBamFile = None
            for tech,mappingQuality in zip(techs,MQs):
                bam = filepath(tech)
                def outputFile(root, name):
                    return os.path.join(root, "%s.chr%s.%s.%s" % ( sample, chr, tech, name ))

                def MismatchIntervals(bam, outputFile, intervals): 
                    return config.gatkCmd('MismatchIntervals') + " -o " + outputFile + " -L " + intervals + " -I " + bam
                def IndelIntervals(bam, outputFile, intervals): 
                    return config.gatkCmd('IndelIntervals') + " -o " + outputFile + " -L " + intervals + " -I " + bam
                def MergeIntervals(bam, files, outputFile, intervals): 
                    return config.gatkCmd('MergeIntervals') + " -o " + outputFile + " ".join(map( lambda x: " -intervals " + x, files )) + " -L " + intervals + " -I " + bam
                def CleanIntervals(bam, outputFile, intervals, snpfile): 
                    return config.gatkCmd('IntervalCleaner') + " -O " + outputFile + " -L " + intervals + " -I " + bam
                def Injector(bam, outputFile, intervals, inputfile): 
                    return config.gatkCmd('CleanedReadInjector') + " --output_bam " + outputFile + " -L " + intervals + " -I " + bam + " --cleaned_reads " + inputfile
                def IndelCaller(bam, outputFile, intervals, fraction): 
                    return config.gatkCmd('IndelGenotyper') + " -O " + outputFile + " -L " + intervals + " -I " + bam + " -minFraction " + fraction
                def SnpCaller(bam, outputFile, intervals): 
                    return config.gatkCmd('SingleSampleGenotyper') + " -o " + outputFile + " -L " + intervals + " -I " + bam
                def VarFiltration(bam, outputHead, intervals, snpcalls, badsnps, indelcalls, depth, mq): 
                    return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -L " + intervals + " -I " + bam + " -B variant,Variants," + snpcalls + ",cleaned,CleanedOutSnp," + badsnps + ",indels,SimpleIndel," + indelcalls + " -X DepthOfCoverage:" + depth + " -X MappingQualityZero:" + mq
                def VarFiltration454(bam, outputHead, intervals, snpcalls, depth, mq): 
                    return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -L " + intervals + " -I " + bam + " -B variant,Variants," + snpcalls + " -X DepthOfCoverage:" + depth + " -X MappingQualityZero:" + mq


                jobid = None
                indelsFileLow = None

                if ( tech != "454" ):
                    mismatchIntervalsFile = outputFile(intervals_dir, "mismatches.intervals")
                    cmd = MismatchIntervals(bam, mismatchIntervalsFile, str(chr))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, mismatchIntervalsFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    indelIntervalsFile = outputFile(intervals_dir, "indels.intervals")
                    cmd = IndelIntervals(bam, indelIntervalsFile, str(chr))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, indelIntervalsFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    mergedIntervalsFile = outputFile(intervals_dir, "merged.intervals")
                    cmd = MergeIntervals(bam, [mismatchIntervalsFile, indelIntervalsFile], mergedIntervalsFile, str(chr))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, mergedIntervalsFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    cleanedFile = outputFile(cleaner_output, "bam")
                    badsnpsFile = outputFile(cleaner_output, "realigner.badsnps")
                    cmd = CleanIntervals(bam, cleanedFile, mergedIntervalsFile, badsnpsFile)
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, cleanedFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    injectedFile = outputFile(injector_output, "bam")
                    cmd = Injector(bam, injectedFile, str(chr), cleanedFile)
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, injectedFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    if ( tech == "SLX" ):
                        solexaBamFile = injectedFile
                    if ( tech == "SOLID" ):
                        solidBamFile = injectedFile

                    cmd = "samtools index " + injectedFile
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, injectedFile + ".bai", just_print_commands = OPTIONS.dry, waitID = jobid)

                    indelsFileLow = outputFile(indel_output, "low.calls")
                    cmd = IndelCaller(injectedFile, indelsFileLow, str(chr), "0.1")
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, indelsFileLow, just_print_commands = OPTIONS.dry, waitID = jobid)

                    indelsFileHigh = outputFile(indel_output, "high.calls")
                    cmd = IndelCaller(injectedFile, indelsFileHigh, str(chr), "0.3")
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, indelsFileHigh, just_print_commands = OPTIONS.dry, waitID = jobid)

                    snpsFile = outputFile(snp_output, "calls")
                    cmd = SnpCaller(injectedFile, snpsFile, str(chr))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, snpsFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    varFiltFile = os.path.join(filter_output, "%s.chr%s.%s" % ( sample, chr, tech ))
                    cmd = VarFiltration(injectedFile, varFiltFile, str(chr), snpsFile, badsnpsFile, indelsFileLow, str(DoC), str(mappingQuality))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, varFiltFile, just_print_commands = OPTIONS.dry, waitID = jobid)
                else:
                    snpsFile = outputFile(snp_output, "calls")
                    cmd = SnpCaller(bam, snpsFile, str(chr))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, snpsFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    varFiltFile = os.path.join(filter_output, "%s.chr%s.%s" % ( sample, chr, tech ))
                    cmd = VarFiltration454(bam, varFiltFile, str(chr), snpsFile, str(DoC), str(mappingQuality))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, varFiltFile, just_print_commands = OPTIONS.dry, waitID = jobid)



            def outputFile(root, techs, name):
                return os.path.join(root, "%s.chr%s.%s.%s" % ( sample, chr, techs, name ))
            def SnpCaller(bams, outputFile, intervals): 
                return config.gatkCmd('SingleSampleGenotyper') + " -o " + outputFile + " -L " + intervals + " ".join(map( lambda x: " -I " + x, bams ))
            def VarFiltration(bams, outputHead, intervals, snpcalls, badsnps, indelcalls, depth, mq): 
                return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -L " + intervals + " -B variant,Variants," + snpcalls + ",cleaned,CleanedOutSnp," + badsnps + ",indels,SimpleIndel," + indelcalls + " -X DepthOfCoverage:" + depth + " -X MappingQualityZero:" + mq + " ".join(map( lambda x: " -I " + x, bams ))

            solid454SnpsFile = outputFile(snp_output, "454-SOLID", "calls")
            cmd = SnpCaller([solidBamFile,filepath("454")], solid454SnpsFile, str(chr))
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, solid454SnpsFile, just_print_commands = OPTIONS.dry, waitID = jobid)

            allSnpsFile = outputFile(snp_output, "allTechs", "calls")
            cmd = SnpCaller([solexaBamFile,solidBamFile,filepath("454")], allSnpsFile, str(chr))
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, allSnpsFile, just_print_commands = OPTIONS.dry, waitID = jobid)

