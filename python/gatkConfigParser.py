#
# GATK configuration parser
# 
import ConfigParser
import os.path
import sys
        
defaultRequiredOptions = {}
def addRequiredOption(name, type):
    defaultRequiredOptions[name] = type
    
addRequiredOption('jar', 'input_file')
addRequiredOption('reference', 'input_file')
addRequiredOption('referenceIndex', 'input_file')
addRequiredOption('referenceDict', 'input_file')
addRequiredOption('java', 'file')
addRequiredOption('jvm_args', str)
addRequiredOption('args', str)
addRequiredOption('tmp', 'output_file')

class gatkConfigParser(ConfigParser.SafeConfigParser):
    GATK = 'DEFAULT'

    def __init__(self, configFiles):
        ConfigParser.SafeConfigParser.__init__(self)
        files = filter(None, configFiles)
        print 'Reading configuration file(s):', files
        print self.read(files)
        self.validateRequiredOptions()

    def validateRequiredOptions(self):
        for key, value in defaultRequiredOptions.iteritems():
            self.validateOption(self.GATK, key, value)
    
    def validateOption(self, section, name, type = str):
        v = self.getOption(section, name, type)
        print '  => Validated option', name, v

    def getGATKOption(self, name, type = str):
        return self.getOption(self.GATK, name, type) 

    def getGATKModeOption(self, name, mode, type = str):
        return self.getOption(mode, name, type) 
        
    def getOption(self, section, name, typeF = None):
        if not self.has_option(section, name):
            raise "Option %s not found in section %s" % (name, section)
        else:
            val = self.get(section, name)
            if typeF == 'input_file' or typeF == 'output_file':
                path = os.path.abspath(os.path.expanduser(val))
                if typeF == 'input_file':
                    if not os.path.exists(path):
                        raise "Input file does not exist", path
                    if not os.access(path, os.R_OK):
                        raise "Input file cannot be read", path
                if typeF == 'output_file':
                    if not os.access(path, os.W_OK):
                        raise "Output file cannot be written", path
                return path
            elif type(typeF) == str:
                return str(val)
            elif typeF == None:
                return val
            else:
                return typeF(val)
        
    def java(self): return self.getOption(self.GATK, 'java')
    def jvm_args(self): return self.getOption(self.GATK, 'jvm_args')
    def jar(self): return self.getOption(self.GATK, 'jar')
    def gatk_args(self): return self.getOption(self.GATK, 'args')
    def reference(self): return self.getOption(self.GATK, 'reference')
    
    def gatkCmd(self, mode):
        cmd = ' '.join([self.java(), self.jvm_args(), '-jar', self.jar(), self.gatk_args(), '-R', self.reference()])
        cmd += ' ' + ' '.join(['-T', mode, self.getGATKModeOption('args', mode)])
        return cmd
        
import unittest
class TestMergeBAMsUtils(unittest.TestCase):
    def setUp(self):
        configFile = os.path.join(os.path.split(sys.argv[0])[0] + "/../testdata/defaultGATKConfig.cfg")
        self.config = gatkConfigParser(configFile)

    def testValidate(self):
        self.config.validateRequiredOptions()

    def testCmd(self):
        s = "java -ea -Xmx2048m -jar ~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -l INFO -L 1:1-10,000,000 -R /home/radon01/depristo/work/humanref/Homo_sapiens_assembly18.fasta -T CountCovariates --CREATE_TRAINING_DATA --MIN_MAPPING_QUALITY 1"
        self.assertEquals(self.config.gatkCmd('CountCovariates'), s)

if __name__ == '__main__':
    unittest.main()    