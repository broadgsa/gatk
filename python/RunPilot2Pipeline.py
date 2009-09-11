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
        
# Official genome-wide Depth of Coverage tables for pilot 2, freeze 5:
#        NA12878 NA12891 NA12892 NA19238 NA19239 NA19240
# 454:     36                                      18
# SLX:     82      91      70      56      68      86
# SOLID:   37                                      64
# 454+SLD: 64                                      77
# ALL:     138                                     150
    DoC_454 = {"NA12878":36,"NA19240":18}
    DoC_slx = {"NA12878":82,"NA12891":91,"NA12892":70,"NA19238":56,"NA19239":68,"NA19240":86}
    DoC_solid = {"NA12878":37,"NA19240":64}
    DoC_454solid = {"NA12878":64,"NA19240":77}
    DoC_all = {"NA12878":138,"NA19240":150}
    DoC_hash = {"454":DoC_454,"SLX":DoC_slx,"SOLID":DoC_solid,"454SOLID":DoC_454solid,"ALL":DoC_all}
    MQ_hash = {"SLX":100,"SOLID":5,"454":5,"454SOLID":10,"ALL":110}

    intervals_dir = outputDir("intervals")
    cleaner_output = outputDir("cleaner")
    injector_output = outputDir("bams")
    snp_output = outputDir("calls/unfiltered_snps")
    filter_output = outputDir("calls/filtered_snps")
    indel_output = outputDir("calls/indels")
    final_bam_dir = outputDir("useTheseBamsForAnalyses")
    #final_bam_dir = "/humgen/gsa-hphome1/projects/1kg_pilot2/useTheseBamsForAnalyses"

    samples = ["NA12878","NA12891","NA12892","NA19238","NA19239","NA19240"]
    techs = ["SLX"]
    chrs = range(1, 23) + ["X"]

    for sample in samples:
        #
        # Actually do some work here
        #
        def finalBam(tech):
            return os.path.join(final_bam_dir, "%s.%s.bam" % ( sample, tech ))
        def outputFileTech(root, tech, name):
            return os.path.join(root, "%s.%s.%s" % ( sample, tech, name ))
        def badSnps( tech ):
            return os.path.join(cleaner_output, "%s.%s.realigner.badsnps" % ( sample, tech ))
        def indelsForFiltering( tech ):
            return outputFileTech(indel_output, tech, "low.calls")

        myTechs = techs
        if sample in ["NA12878", "NA19240"]:
            myTechs = techs + ["SOLID","454"]

        for tech in myTechs:

            if ( tech != "454" ):
                myChrs = chrs
                if sample in ["NA12891", "NA19239"]:
                    myChrs = chrs + ["Y"]
                def badSnpsChr( tech, chr ):
                    return os.path.join(cleaner_output, "%s.chr%s.%s.realigner.badsnps" % ( sample, chr, tech ))

                def makeJobName(suffix):
                    return sample + "." + tech + "." + suffix
                def makeJobClass(suffix):
                    return sample + ".*." + suffix

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

                    mismatchIntervalsFile = outputFile(intervals_dir, "mismatches.intervals")
                    cmd = MismatchIntervals(bam, mismatchIntervalsFile, str(chr))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, mismatchIntervalsFile, just_print_commands = OPTIONS.dry, waitID = None, jobName = makeJobName("phase1." + str(chr)))

                    indelIntervalsFile = outputFile(intervals_dir, "indels.intervals")
                    cmd = IndelIntervals(bam, indelIntervalsFile, str(chr))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, indelIntervalsFile, just_print_commands = OPTIONS.dry, waitID = None, jobName = makeJobName("phase1." + str(chr)))

                    mergedIntervalsFile = outputFile(intervals_dir, "merged.intervals")
                    cmd = MergeIntervals(bam, [mismatchIntervalsFile, indelIntervalsFile], mergedIntervalsFile, str(chr))
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, mergedIntervalsFile, just_print_commands = OPTIONS.dry, waitID = makeJobName("phase1." + str(chr)))

                    cleanedFile = outputFile(cleaner_output, "bam")
                    badsnpsFile = badSnpsChr(tech, str(chr))
                    cmd = CleanIntervals(bam, cleanedFile, mergedIntervalsFile, badsnpsFile)
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, cleanedFile, just_print_commands = OPTIONS.dry, waitID = jobid)
                    injectedFile = outputFile(injector_output, "bam")
                    cmd = Injector(bam, injectedFile, str(chr), cleanedFile)
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, injectedFile, just_print_commands = OPTIONS.dry, waitID = jobid)

                    cmd = "samtools index " + injectedFile
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, injectedFile + ".bai", just_print_commands = OPTIONS.dry, waitID = jobid, jobName = makeJobName("phase2"))

                def MergeBams(outputFile): 
                    return "MergeBAMBatch.py -d " + cleaner_output + " -q " + OPTIONS.farmQueue + " -s '" + outputFile + "\t" + os.path.join(cleaner_output, "%s.chr*.%s.bam" % ( sample, tech )) + "'"
 
                cmd = MergeBams(finalBam(tech))
                mergeJobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, finalBam(tech), just_print_commands = OPTIONS.dry, waitID = makeJobName("phase2"), jobName = makeJobName("phase3"))

                cmd = "cat "
                for chr in myChrs:
                    cmd = cmd + " " + badSnpsChr(tech, chr)
                cmd = cmd + " > " + badSnps(tech)
                badsnpsJobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, badSnps(tech), just_print_commands = OPTIONS.dry, waitID = makeJobName("phase2"), jobName = makeJobName("phase4"))

            def IndelCaller(bam, outputFile, fraction): 
                return config.gatkCmd('IndelGenotyper') + " -O " + outputFile + " -I " + bam + " -minFraction " + fraction
            def SnpCaller(bam, outputFile): 
                return config.gatkCmd('SingleSampleGenotyper') + " -varout " + outputFile + " -I " + bam
            def VarFiltration(bam, outputHead, snpcalls, badsnps, indelcalls, depth, mq): 
                return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -I " + bam + " -B variant,Variants," + snpcalls + ",cleaned,CleanedOutSNP," + badsnps + ",indels,SimpleIndel," + indelcalls + " -X DepthOfCoverage:max=" + depth + " -X MappingQualityZero:max=" + mq
            def VarFiltration454(bam, outputHead, snpcalls, depth, mq): 
                return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -I " + bam + " -B variant,Variants," + snpcalls + " -X DepthOfCoverage:max=" + depth + " -X MappingQualityZero:max=" + mq

            waitid = makeJobName("phase3")
            if ( tech == "454" ):
                waitid = None

            bamToCallFrom = finalBam(tech)
            indelsFileHigh = outputFileTech(indel_output, tech, "high.calls")
            cmd = IndelCaller(bamToCallFrom, indelsFileHigh, "0.3")
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, indelsFileHigh, just_print_commands = OPTIONS.dry, waitID = waitid)

            indelsFileLow = indelsForFiltering(tech)
            cmd = IndelCaller(bamToCallFrom, indelsFileLow, "0.1")
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, indelsFileLow, just_print_commands = OPTIONS.dry, waitID = waitid, jobName = makeJobName("phase4"))

            snpsFile = outputFileTech(snp_output, tech, "calls")
            cmd = SnpCaller(bamToCallFrom, snpsFile)
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, snpsFile, just_print_commands = OPTIONS.dry, waitID = waitid, jobName = makeJobName("phase4"))

            varFiltFile = os.path.join(filter_output, "%s.%s" % ( sample, tech ))
            if ( tech != "454" ):
                cmd = VarFiltration(bamToCallFrom, varFiltFile, snpsFile, badSnps(tech), indelsFileLow, str(DoC_hash[tech][sample]), str(MQ_hash[tech]))
            else:
                cmd = VarFiltration454(bamToCallFrom, varFiltFile, snpsFile, str(DoC_hash[tech][sample]), str(MQ_hash[tech]))
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, varFiltFile, just_print_commands = OPTIONS.dry, waitID = makeJobName("phase4"))

        def SnpCaller(bams, outputFile): 
            return config.gatkCmd('SingleSampleGenotyper') + " -varout " + outputFile + " ".join(map( lambda x: " -I " + x, bams ))
        def VarFiltration(bams, outputHead, snpcalls, badsnps, indelcalls, depth, mq): 
            return config.gatkCmd('VariantFiltration') + " -VOH " + outputHead + " -B variant,Variants," + snpcalls + ",cleaned,CleanedOutSNP," + badsnps + ",indels,SimpleIndel," + indelcalls + " -X DepthOfCoverage:max=" + depth + " -X MappingQualityZero:max=" + mq + " ".join(map( lambda x: " -I " + x, bams ))

        if sample in ["NA12878", "NA19240"]:

            solid454SnpsFile = outputFileTech(snp_output, "454-SOLID", "calls")
            cmd = SnpCaller([finalBam("SOLID"),finalBam("454")], solid454SnpsFile)
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, solid454SnpsFile, just_print_commands = OPTIONS.dry, waitID = makeJobClass("phase3"))

            solid454VarFiltFile = os.path.join(filter_output, "%s.454-SOLID" % ( sample ))
            cmd = VarFiltration([finalBam("SOLID"),finalBam("454")], solid454VarFiltFile, solid454SnpsFile, badSnps("SOLID"), indelsForFiltering("SOLID"), str(DoC_hash["454SOLID"][sample]), str(MQ_hash["454SOLID"]))
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, solid454VarFiltFile, just_print_commands = OPTIONS.dry, waitID = jobid)

            allSnpsFile = outputFileTech(snp_output, "allTechs", "calls")
            cmd = SnpCaller([finalBam("SLX"),finalBam("SOLID"),finalBam("454")], allSnpsFile)
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, allSnpsFile, just_print_commands = OPTIONS.dry, waitID = makeJobClass("phase3"))
            allVarFiltFile = os.path.join(filter_output, "%s.allTechs" % ( sample ))
            cmd = VarFiltration([finalBam("SLX"),finalBam("SOLID"),finalBam("454")], allVarFiltFile, allSnpsFile, badSnps("SLX"), indelsForFiltering("SLX"), str(DoC_hash["ALL"][sample]), str(MQ_hash["ALL"]))
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, allVarFiltFile, just_print_commands = OPTIONS.dry, waitID = jobid)

