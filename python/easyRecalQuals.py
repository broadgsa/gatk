import farm_commands
import os.path
import sys
from optparse import OptionParser
from gatkConfigParser import *
import glob

if __name__ == "__main__":
    usage = """usage: %prog [-c config.cfg]* input.bam output.bam"""

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
    parser.add_option("-c", "--config", dest="configs",
                        action="append", type="string", default=[],
                        help="Configuration file")                        
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
    if len(args) != 2:
        parser.error("incorrect number of arguments")

    config = gatkConfigParser(OPTIONS.configs)
    inputBAM = args[0]
    outputBAM = args[1]
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

    #
    # Actually do some work here
    #
    jobid = None
    jobid = runCovariateCmd(inputBAM, initDataFile, covariateInitial, jobid, False)
    
    if OPTIONS.ignoreExistingFiles or not os.path.exists(outputBAM):
        cmd = recalibrateCmd(inputBAM, initDataFile, outputBAM)
        jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry, waitID = jobid)
        jobid = farm_commands.cmd('samtools index ' + outputBAM, OPTIONS.farmQueue, None, just_print_commands = OPTIONS.dry, waitID = jobid)
        
    jobid = runCovariateCmd(outputBAM, recalDataFile, covariateRecal, jobid, True)
        