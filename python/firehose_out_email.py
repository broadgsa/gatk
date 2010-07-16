#this script produces the output to go in emails 
import os
import re
import sample_lister

class SampSet(sample_lister.SampleSet):
    def __init__(self, sampset, pathname):
        self.sampset=sampset
        self.pathname=pathname
        sample_lister.SampleSet.__init__(self, sampset, pathname)
    def evalout(self):
        '''This produced the output that needs to go in the emails'''
        filename = "/humgen/gsa-firehose/firehose/firehose_output/trunk/Sample_Set/" + self.sampset +"/UnifiedGenotyper/"+ self.sampset+".filtered.eval"
        evalfile = open(filename, "r").read()
        annotations= ["all", "novel", "known", "snp_at_known_non_snps", "filtered"]
        variant=dict(zip(annotations, ('','','','','')))
        ratio=dict(zip(annotations, ('','','','','')))
        bpre=re.compile("all,summary,variant_counts +n bases covered +(\d+)")
        size=repr(bpre.search(evalfile).group(1))
        for a in annotations:
            anregexv  = re.compile(a + ",summary,variant_counts +variants +(\d+)")
            variant[a] = repr(anregexv.search(evalfile).group(1))
            anregexr = re.compile(a + ",summary,transitions_transversions +ratio +(\d+.\d+)")
            ratio[a] = repr(anregexr.search(evalfile).group(1)
        print("Samples processed:\n\n Target size: \t" +size+"  bp \n\n\t\t\t\t\t Variants \t\t Ti/TV \n (true positives)\t All \t\t " +variant["all"]+ " \t\t " + ratio["all"] +" \n \t\t\t Known \t\t " +variant["known"]+ " \t\t " + ratio['known']+" \n \t\t\t Novel \t\t " +variant["novel"]+" \t\t " + ratio['novel']+  " \n*************************************************************************\n (false  \tSNPS at known indels \t " +variant["snp_at_known_non_snps"]+"\t\t\t " + ratio['snp_at_known_non_snps']+ " \n positives) \t\t filtered \t " +variant["filtered"]+" \t\t " + ratio['filtered'] )


#EOMI=SampSet("EOMI_Kathiresan_NHGRI", "test")
#EOMI.evalout() <-this and the line above are examples

#TODO: make this send the email when run with a setname as input
#TODO: make this find the list of bams, bed files, and annotated vcfs.
