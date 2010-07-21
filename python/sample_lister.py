import os
import subprocess

#this will create a list of collaborator ids or a list of sample names--it will find cleaned bams even if they're not stored in the sample_set file system. I'm mainly using it as a module import. 
class SampleSet:
    def __init__(self, projectname, setname, path):
        self.setname = setname
        self.projectname = projectname
        self.path = path
    def sampslist(self):
        '''finds and lists all samples in a set'''
        try:
            searchpath="ls /humgen/gsa-firehose/firehose/firehose_output/trunk/Sample_Set/" + self.setname +"/ -I CleanBam -I UnifiedGenotyper -I MergeBam"  
            raw_samps= subprocess.Popen([searchpath], shell=True, stdout=subprocess.PIPE).communicate()[0]
        except IOError:
            print( "Can't make sample list. Those files are not where they ought to be, or Sample_Set is not valid.")
        samps = raw_samps.split("\n" + self.projectname+ "_")
        samplelist = raw_samps.split("\n")[0:len(samps)]
        samps[0] = samps[0].split(self.projectname+"_")[len(samps[0].split(self.projectname+"_"))-1]
        samps[len(samps)-1] = samps[len(samps)-1].split("\n")[0]
        return [samps, samplelist]
    def bamlist(self, samplist, write=True):
        '''finds and lists all cleaned bams in a sample set'''
        if (write == True):
            try:
                if os.path.exists(self.path + "bamsfor" + self.setname + ".list"):
                    os.remove(self.path + "bamsfor" + self.setname + ".list")
                listfile = open(self.path + "bamsfor" + self.setname + ".list", "a")
                for samp in samplist:
                    searcher="find /humgen/gsa-firehose/firehose/firehose_output/trunk/Sample/" + repr(samp)  +"/ -name \*cleaned.bam"
                    raw_samp= subprocess.Popen([searcher], shell=True, stdout=subprocess.PIPE).communicate()[0]
                    listfile.write(raw_samp)
                listfile.close()
                print (listfile.name)
            except IOError:
                print( "can't make .bam list.Those files are not where they ought to be, or Sample_Set is not valid")
        else:
            for samp in samplist:
                searcher="find /humgen/gsa-firehose/firehose/firehose_output/trunk/Sample/" + samp  +"/ -name \*cleaned.bam"
                raw_samp= subprocess.Popen([searcher], shell=True, stdout=subprocess.PIPE).communicate()[0]
                print(raw_samp)
    def bedlist(self, samplist, write=True):
        '''finds and lists all beds for a sample set'''
        if (write == True):
            try:
                if os.path.exists(self.path + "bedsfor" + self.setname + ".list"):
                    os.remove(self.path + "bedsfor" + self.setname + ".list")
                listfile = open(self.path + "bedsfor" + self.setname + ".list", "a")
                for samp in samplist:
                    searcher="find /humgen/gsa-firehose/firehose/firehose_output/trunk/Sample/" + repr(samp)  +"/ -name \*.bed"
                    raw_samp= subprocess.Popen([searcher], shell=True, stdout=subprocess.PIPE).communicate()[0]
                    listfile.write(raw_samp)
                listfile.close()
                print (listfile.name)
            except IOError:
                print( "can't make .bed list.Those files are not where they ought to be, or Sample_Set is not valid")
        else:
            for samp in samplist:
                searcher="find /humgen/gsa-firehose/firehose/firehose_output/trunk/Sample/" + samp  +"/ -name \*.bed"
                raw_samp= subprocess.Popen([searcher], shell=True, stdout=subprocess.PIPE).communicate()[0]
                print(raw_samp)
'''next two lines are example usage
#pfizer5=SampleSet("T2D_Altshuler_Pfizer_Plate_5", "T2D_Altshuler_Pfizer", "humgen/gsa-hphome1/corin/oneoffs/pfizer5/")
#pfizer5.bamlist(pfizer5.sampslist()[0], write=False)'''
