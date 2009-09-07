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
    final_bam_dir = "/humgen/gsa-hphome1/projects/1kg_pilot2/useTheseBamsForAnalyses"

    samples = ["NA12878","NA12891","NA12892","NA19238","NA19239","NA19240"]
    techs = ["SLX"]
    chrs = range(1, 22) + ["X"]
    DoCs = [82,91,70,56,68,86]

# Official genome-wide Depth of Coverage tables for pilot 2, freeze 5:
#        NA12878 NA12891 NA12892 NA19238 NA19239 NA19240
# 454:     36                                      18
# SLX:     82      91      70      56      68      86
# SOLID:   37                                      64
# 454+SLD: 64                                      77
# ALL:     xx                                      xx

    for sample, DoC in zip(samples, DoCs):
        #
        # Actually do some work here
        #
        MQs = [100,5,5]
        def finalBam(tech):
            return os.path.join(final_bam_dir, "%s.%s.bam" % ( sample, tech ))
        def outputFile(root, tech, name):
            return os.path.join(root, "%s.%s.%s" % ( sample, tech, name ))
        def badSnps( tech ):
            return outputFile(cleaner_output, tech, "realigner.badsnps")
        def indelsForFiltering( tech ):
            return outputFile(indel_output, tech, "low.calls")

        myTechs = techs
        if sample in ["NA12878", "NA19240"]:
            myTechs = techs + ["SOLID","454"]

        for tech,mappingQuality in zip(myTechs,MQs):

            myChrs = chrs
            if sample in ["NA12891", "NA19239"]:
                myChrs = chrs + ["Y"]
            def badSnps( tech, chr ):
                return os.path.join(cleaner_output, "%s.chr%s.%s.realigner.badsnps" % ( sample, chr, tech ))

            for chr in myChrs:

                bam = "/broad/1KG/DCC/ftp/pilot_data/%s/alignment/%s.chrom%s.%s.SRP000032.2009_07.bam" % ( sample, sample, chr, tech )

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

                jobid = None
                mergeid = None
                bamToCallFrom = finalBam(tech)

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
                    badsnpsFile = badSnps(tech)
                    cmd = CleanIntervals(bam, cleanedFile, mergedIntervalsFile, badsnpsFile)
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, cleanedFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    injectedFile = outputFile(injector_output, "bam")
                    cmd = Injector(bam, injectedFile, str(chr), cleanedFile)
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, injectedFile, just_print_commands = OPTIONS.dry, waitID = jobid)

# HOW DO I CALL MergeBAMBatch FROM HERE?
#                def MergeBams(): 
#                    return "MergeBAMBatch.py -d " + final_bam_dir + " -q " + OPTIONS.farmQueue + " -s '" + bamToCallFrom + "\t" + os.path.join(injector_output, "%s.chr*.%s.bam" % ( sample, tech )) + "'"
#                    cmd = MergeBams()
#                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, bamToCallFrom, just_print_commands = OPTIONS.dry, waitID = jobid)
                    mergeid = jobid

            cmd = "cat "
            for chr in myChrs:
                cmd = cmd + " " + badSnps(tech, chr)
            cmd = cmd + " > " + badSnps(tech)
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, badSnps(tech), just_print_commands = OPTIONS.dry, waitID = mergeid)

            def IndelCaller(bam, outputFile, fraction): 
                return config.gatkCmd('IndelGenotyper') + " -O " + outputFile + " -I " + bam + " -minFraction " + fraction
            def SnpCaller(bam, outputFile): 
                return config.gatkCmd('SingleSampleGenotyper') + " -o " + outputFile + " -I " + bam
            def VarFiltration(bam, outputHead, snpcalls, badsnps, indelcalls, depth, mq): 
                return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -I " + bam + " -B variant,Variants," + snpcalls + ",cleaned,CleanedOutSnp," + badsnps + ",indels,SimpleIndel," + indelcalls + " -X DepthOfCoverage:" + depth + " -X MappingQualityZero:" + mq
            def VarFiltration454(bam, outputHead, snpcalls, depth, mq): 
                return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -I " + bam + " -B variant,Variants," + snpcalls + " -X DepthOfCoverage:" + depth + " -X MappingQualityZero:" + mq


            indelsFileHigh = outputFile(indel_output, tech, "high.calls")
            cmd = IndelCaller(bamToCallFrom, indelsFileHigh, "0.3")
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, indelsFileHigh, just_print_commands = OPTIONS.dry, waitID = mergeid)

            indelsFileLow = indelsForFiltering(tech)
            cmd = IndelCaller(bamToCallFrom, indelsFileLow, "0.1")
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, indelsFileLow, just_print_commands = OPTIONS.dry, waitID = mergeid)

            snpsFile = outputFile(snp_output, tech, "calls")
            cmd = SnpCaller(bamToCallFrom, snpsFile)
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, snpsFile, just_print_commands = OPTIONS.dry, waitID = jobid)  # wait on the low indel calls

            varFiltFile = os.path.join(filter_output, "%s.%s" % ( sample, tech ))
            if ( tech != "454" ):
                cmd = VarFiltration(bamToCallFrom, varFiltFile, snpsFile, badsnpsFile, indelsFileLow, str(DoC), str(mappingQuality))
            else:
                cmd = VarFiltration454(bamToCallFrom, varFiltFile, snpsFile, str(DoC), str(mappingQuality))
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, varFiltFile, just_print_commands = OPTIONS.dry, waitID = jobid)

        def SnpCaller(bams, outputFile): 
            return config.gatkCmd('SingleSampleGenotyper') + " -o " + outputFile + " ".join(map( lambda x: " -I " + x, bams ))
        def VarFiltration(bams, outputHead, snpcalls, badsnps, indelcalls, depth, mq): 
            return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -B variant,Variants," + snpcalls + ",cleaned,CleanedOutSnp," + badsnps + ",indels,SimpleIndel," + indelcalls + " -X DepthOfCoverage:" + depth + " -X MappingQualityZero:" + mq + " ".join(map( lambda x: " -I " + x, bams ))

#
# HOW DO I MAKE THESE JOBS DEPEND ON THE MERGE IDS OF THE INDIVIDUAL SAMPLES???
# (Or until everything else is done?)
#
        solid454SnpsFile = outputFile(snp_output, "454-SOLID", "calls")
        cmd = SnpCaller([finalBam("SOLID"),finalBam("454")], solid454SnpsFile)
        jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, solid454SnpsFile, just_print_commands = OPTIONS.dry, waitID = allMergeIds)

        solid454VarFiltFile = os.path.join(filter_output, "%s.454-SOLID" % ( sample ))
        cmd = VarFiltration([finalBam("SOLID"),finalBam("454")], solid454VarFiltFile, solid454SnpsFile, badSnps("SOLID"), indelsForFiltering("SOLID"), str(DoC), str(mappingQuality))
        jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, allVarFiltFile, just_print_commands = OPTIONS.dry, waitID = jobid)

        allSnpsFile = outputFile(snp_output, "allTechs", "calls")
        cmd = SnpCaller([finalBam("SLX"),finalBam("SOLID"),finalBam("454")], solid454SnpsFile)
        jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, allSnpsFile, just_print_commands = OPTIONS.dry, waitID = allMergeIds)
        allVarFiltFile = os.path.join(filter_output, "%s.allTechs" % ( sample ))
        cmd = VarFiltration([finalBam("SLX"),finalBam("SOLID"),finalBam("454")], varFiltFile, allSnpsFile, badSnps("SLX"), indelsForFiltering("SLX"), str(DoC), str(mappingQuality))
        jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, allVarFiltFile, just_print_commands = OPTIONS.dry, waitID = jobid)
#
# HOW DO I COMBINE ALL OF THE MQs AND DOCs?
#


