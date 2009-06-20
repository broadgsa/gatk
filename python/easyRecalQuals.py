#
#
# 
# java -Xmx2048m -jar ~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -R ~/work/humanref/Homo_sapiens_assembly18.fasta --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -l INFO -T CountCovariates -I NA07048.bam --OUTPUT_FILEROOT test/NA07048_test --CREATE_TRAINING_DATA --MIN_MAPPING_QUALITY 1 -L chr1:1-10,000,000 -collapsePos -collapseDinuc && cat test/NA07048_test.raw_data.csv
# 
# java -Xmx2048m -jar ~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -R ~/work/humanref/Homo_sapiens_assembly18.fasta -T TableRecalibration -I NA07048.bam -params test/NA07048_test.raw_data.csv -outputBAM NA07048.test.bam -l INFO -compress 1 -L chr1:1-10,000,000
# samtools index NA07048.test.bam
# 
# java -Xmx2048m -jar ~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -R ~/work/humanref/Homo_sapiens_assembly18.fasta --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -l INFO -T CountCovariates -I NA07048.test.bam --OUTPUT_FILEROOT test/NA07048_test.recal --CREATE_TRAINING_DATA --MIN_MAPPING_QUALITY 1 -L chr1:1-10,000,000 -collapsePos -collapseDinuc && cat test/NA07048_test.recal.raw_data.csv
# 
import farm_commands
import os.path
import sys
from optparse import OptionParser
import picard_utils
from gatkConfigParser import *
import glob

if __name__ == "__main__":
    usage = """usage: %prog config.cfg* input.bam output.bam"""

    parser = OptionParser(usage=usage)
    parser.add_option("-A", "--args", dest="args",
                        type="string", default="",
                        help="arguments to GATK")
    parser.add_option("-C", "--CovariateArgs", dest="CovariateArgs",
                        type="string", default="",
                        help="arguments to GATK")
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
    parser.add_option("-d", "--dir", dest="scratchDir",
                        type="string", default="test",
                        help="Output directory")
    parser.add_option("", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, nothing actually gets run, just a dry run")
    parser.add_option("", "--plot", dest="plot",
                        action='store_true', default=False,
                        help="If provided, will call R to generate convenient plots about the Q scores of the pre and post calibrated files")
    parser.add_option("-i", "--ignoreExistingFiles", dest="ignoreExistingFiles",
                        action='store_true', default=False,
                        help="Ignores already written files, if present")

    (OPTIONS, args) = parser.parse_args()
    #if len(args) != 3:
    #    parser.error("incorrect number of arguments")

    configFiles = args[0:len(args)-2]
    config = gatkConfigParser(configFiles)
    inputBAM = args[len(args)-2]
    outputBAM = args[len(args)-1]
    rootname = os.path.split(os.path.splitext(outputBAM)[0])[1]

    covariateRoot = os.path.join(OPTIONS.scratchDir, rootname)
    covariateInitial = covariateRoot + '.init'
    initDataFile = covariateInitial + '.raw_data.csv'
    covariateRecal = covariateRoot +  '.recal'
    recalDataFile = covariateRecal + '.raw_data.csv'

    if not os.path.exists(OPTIONS.scratchDir):
        os.mkdir(OPTIONS.scratchDir)
    
    def covariateCmd(bam, outputDir, ignoreAdds):
        add = " -I %s --OUTPUT_FILEROOT %s" % (bam, outputDir)
        if not ignoreAdds:
            add += " " + OPTIONS.CovariateArgs
        return config.gatkCmd('CountCovariates') + add 

    def recalibrateCmd(inputBAM, dataFile, outputBAM):
        return config.gatkCmd('TableRecalibration') + " -I %s -params %s -outputBAM %s" % (inputBAM, dataFile, outputBAM)

    def runCovariateCmd(inputBAM, dataFile, dir, jobid, ignoreAdds = False):
        if OPTIONS.ignoreExistingFiles or not os.path.exists(dataFile):
            cmd = covariateCmd(inputBAM, dir, ignoreAdds)
            return farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry, waitID = jobid)

    if OPTIONS.plot:
        def plotCmd(cmd):
            farm_commands.cmd(cmd, None, None, just_print_commands = OPTIONS.dry)
            
        Rscript = config.getOption('R', 'Rscript', 'input_file')
        plotEmpStated = config.getOption('R', 'PlotQEmpStated', 'input_file')
        plotByCycleDinuc = config.getOption('R', 'PlotQByCycleDinuc', 'input_file')
        
        empQualFiles = map( lambda x: glob.glob(''.join([x,'.RG_*.empirical_v_reported_quality.csv'])), [covariateInitial, covariateRecal])
        empQualFiles = empQualFiles[0] + empQualFiles[1]
        readGroupRoots = map(lambda x: x.replace(".empirical_v_reported_quality.csv", ""), empQualFiles)
        print readGroupRoots
        
        for file in empQualFiles:
            plotCmd(' '.join([Rscript, plotEmpStated, file]))
            
        for root in readGroupRoots:
            plotCmd(' '.join([Rscript, plotByCycleDinuc, root]))
            
    else:
        jobid = None
        jobid = runCovariateCmd(inputBAM, initDataFile, covariateInitial, jobid, False)
        
        if OPTIONS.ignoreExistingFiles or not os.path.exists(outputBAM):
            cmd = recalibrateCmd(inputBAM, initDataFile, outputBAM)
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry, waitID = jobid)
            jobid = farm_commands.cmd('samtools index ' + outputBAM, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry, waitID = jobid)
            
        jobid = runCovariateCmd(outputBAM, recalDataFile, covariateRecal, jobid, True)
        