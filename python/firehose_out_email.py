#this script produces the output to go in emails 
import subprocess
import os
import re
import sys
import getopt
import sample_lister

try:
        opts, args = getopt.getopt(sys.argv[1:], "dp:s:")
except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

opts=dict(opts)
givenname=opts['-s']
projname=opts['-p']
if '-d' in opts:
	dirname=opts['-d']
else:
	dirname=""

class SampSet(sample_lister.SampleSet):
    def __init__(self, projectname, sampset, pathname):
        self.sampset=sampset
        self.pathname=pathname
        self.projectname=projectname
        sample_lister.SampleSet.__init__(self, projectname, sampset, pathname)
    def evalout(self):
        '''This produced the output that needs to go in the emails'''
        filename = "/humgen/gsa-firehose/firehose/firehose_output/trunk/Sample_Set/" + self.sampset +"/UnifiedGenotyper/"+ self.sampset+".filtered.eval"
        evalfile = open(filename, "r").read()
        annotations= ["all", "novel", "known", "snp_at_known_non_snps", "filtered"]
        variant=dict(zip(annotations, ('','','','','')))
        ratio=dict(zip(annotations, ('','','','','')))
        bpre=re.compile("all,summary,variant_counts +n bases covered +(\d+)")
        size=repr(bpre.search(evalfile).group(1))
        bamsearch="find /humgen/gsa-firehose/firehose/firehose_output/trunk/Sample_Set/"+self.sampset+"/* -name \*.bam | grep -v unfixed"
        bams = subprocess.Popen([bamsearch], shell=True, stdout=subprocess.PIPE).communicate()[0]
        sampno=bams.count("bam")
        bedsearch="find /humgen/gsa-firehose/firehose/firehose_output/trunk/Sample_Set/"+self.sampset+"/* -name \*filtered_indels.bed"
        beds = subprocess.Popen([bedsearch], shell=True, stdout=subprocess.PIPE).communicate()[0]
        vcf="/humgen/gsa-firehose/firehose/firehose_output/trunk/Sample_Set/"+self.sampset+"/UnifiedGenotyper/"+self.sampset+'.maf.annotated.vcf'
        for a in annotations:
            anregexv  = re.compile(a + ",summary,variant_counts +variants +(\d+)")
            variant[a] = repr(anregexv.search(evalfile).group(1))
            anregexr = re.compile(a + ",summary,transitions_transversions +ratio +(\d+.\d+|Infinity)")
            ratio[a] = repr(anregexr.search(evalfile).group(1))
        out1="Samples processed:"+repr(sampno)+"\n\n Target size: \t" +size+"  bp \n\n\t\t\t\t\t Variants \t\t Ti/TV \n (true positives)\t All \t\t " +variant["all"]+ " \t\t " + ratio["all"] +" \n \t\t\t Known \t\t " +variant["known"]+ " \t\t " + ratio['known']+" \n \t\t\t Novel \t\t " +variant["novel"]+" \t\t " + ratio['novel']+  " \n*************************************************************************\n (false  \tSNPS at known indels \t " +variant["snp_at_known_non_snps"]+"\t\t\t " + ratio['snp_at_known_non_snps']+ " \n positives) \t\t filtered \t " +variant["filtered"]+" \t\t " + ratio['filtered']
        out2="\n\n\nSNP calls:"+vcf+"\n\nIndel-realigned Bam files:\n"+bams+"\nIndel calls:\n"+beds
        if self.pathname == '':
            print(out1+out2)
        else:
            filename=self.pathname+self.sampset+".emailtxt"
            putthere=open(filename, "w")
            putthere.write(out1+out2)

target=SampSet(projname,givenname,dirname)
target.evalout()
#TODO: make this send the email when run
#TODO: make this find the list of bams, bed files, and annotated vcfs.
