#!/usr/bin/env python

"""Starts from TGTC list of fastb, qualb, and qltout files and creates merged SAM/BAM files"""

import os, commands, sys
from farm_commands import cmd

# Global list of failed conversions to be reported at end of execution
failed_list = []

class FailedConversion:
    def __init__(self, lane, error_str):
        self.lane = lane
        self.error_str = error_str
        print self
    def __str__(self):
        return self.error_str+" while attempting conversion of "+str(self.lane)

class lane:
    bam_ref_list = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.bam.ref_list";    
    pipe_dirs_head = "/slxa/pipelineOutput/farm"
    pipe_dirs = os.listdir(pipe_dirs_head)
    pipe_dirs_dict = dict( [(dir.split("_")[1], dir) for dir in pipe_dirs if "_" in dir] )

    def __init__(self, line):
        fields = line.split(" ")
        if len(fields) == 5:
            self.fdl, self.fastb, self.qualb, self.qltout, group = fields
        else:
            self.fdl, self.fastb, self.qualb, self.qltout = fields
            self.qltout = self.qltout.rstrip()

        self.fqltout = os.path.splitext(self.qltout)[0]+".fais_filtered.qltout"
        self.flowcell, self.date, self.lane = self.fdl.split(".")
        self.flowlane = self.flowcell+"."+self.lane
        self.path_head = os.path.splitext(self.qltout)[0]
        self.sam = self.path_head+".sam"
        self.bam = self.path_head+".bam"
        self.head = ".".join(os.path.basename(self.path_head).split(".")[:2]) # should amount to FLOWCELL.LANE

        #self.refdict_metrics()

    def pipe_dir(self):
        dir = lane.pipe_dirs_dict[self.flowcell]
        return lane.pipe_dirs_head+"/"+dir

    def set_default_refdict(self):
        self.refdict = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.dict";
        #print "WARNING: Using standard Hg18 reference dictionary (refdict) for",self

    def set_default_metrics(self):
        self.metrics = "/broad/1KG/legacy_data/tgtc2sam/unknown.metrics"
        #print "WARNING: Using default unkonwn.metrics file for",self

    def refdict_metrics(self):
        "Ensures that a reference dictionary and metrics are available"
        try:
            pipe_lane_dir = self.pipe_dir()
        except KeyError:
            self.set_default_refdict()
            self.set_default_metrics()
            return True

        self.metrics = pipe_lane_dir+"/"+str(self.flowlane)+".metrics"
        if os.path.exists(self.metrics):
            line = commands.getoutput("grep '^reference=' "+pipe_lane_dir+"/"+self.flowlane+".metrics")
            if line != "":
                reference_fasta = line.split("=")[1]
                self.refdict = os.path.splitext(reference_fasta)[0]+'.dict'
                if "PhiX" in self.refdict:
                    failed_list.append(FailedConversion(self, "PhiX dictionary missing"))
                    return False
            else:
                self.set_default_refdict()
            return True # We were able to generate refdict and metrics files
        else:
            self.set_default_refdict()
            self.set_default_metrics()
            return True # We were able to generate refdict file and substituted unknown.metrics as the metrics file

            #failed_list.append(FailedConversion(self, "Could not find metrics file"))
            #return False 
                                      
    def __str__(self):
        return self.fdl

def usage():
    print "Usage: tgtc_sam.py <command> <command_args>..."
    print
    print "  Works from tgtc files which contain fastb, qualb, qltout filenames corresponding to one lane of data"

    print "  Commands and command arguments:"
    print "    convert tgtc_list - given a tgtc file converts files to SAM and BAM"
    print "      tgtc - file with lines of fastb, qualb, qltout filenames"
    print "      [clean] - optional 3rd argument that will create SAM/BAM files"
    print "                regardless of whether they already exist"
    print
    print "    merge tgtc output_file - merges converted BAM files"
    print "      tgtc - file with lines of fastb, qualb, qltout filenames"
    print "      output_file - path and name of merged BAM file"
    print 
    print "    lists list_of_tgtc_lists - run over a list of tgtc lists given in a file"
    print "      list_of_tgtc_files - file with tgtc filenames in it"
    print
    print "  Global options:"
    print "    print_only: do not run commands, just print what they would be "
    sys.exit(1)

def make_lane_list(tgtc_file):
    all_lanes = [lane(l) for l in open(tgtc_file) if not l.startswith("#")]
    lanes = []
    for ln in all_lanes:
        if ln.refdict_metrics(): # if we were able to generate refdict and metrics
            # append this lane to the lanes list that we will work with
            lanes.append(ln)
    unpaired = [l for l in lanes if ".end1." not in l.qualb and ".end2." not in l.qualb]
    paired1 = [l for l in lanes if ".end1." in l.qualb]
    paired2 = [l for l in lanes if ".end2." in l.qualb]
    
    print "Unpaired lanes:",len(unpaired)
    print "Paired 1 lanes:",len(paired1)
    print "Paired 2 lanes:",len(paired2)
    assert len(paired1) == len(paired2)

    return [(u, None) for u in unpaired] + zip(paired1, paired2)

def convert_lanes(tgtc_file, recreate_all=False, just_print_commands=False):
    lanes_done = 0
    lane_list = make_lane_list(tgtc_file)
    print "Beginning processing of these lanes (unpaired and paired):"
    print "\n".join(["%s %s" % (p1,p2) for p1,p2 in lane_list])
 
    for p1, p2 in lane_list:
        if recreate_all or not os.path.exists(p1.sam) or not os.path.getsize(p1.sam):

            # Filer reads with > 4 errors
            cmd_str = "\"FilterAlignsILTStyle QLTIN="+p1.qltout+" QLTOUT="+p1.fqltout

            if p2 == None:
                # This is an UNPAIRED LANE
                print "Processing unpaired lane",p1
                cmd_str += (" && Qltout2SAM"+
                            " HEAD="+p1.head+
                            " SAM="+p1.sam+
                            " QLTOUT="+p1.fqltout+
                            " FASTB="+p1.fastb+
                            " QUALB="+p1.qualb+
                            " REFDICT="+p1.refdict+
                            " METRICS="+p1.metrics+
                            " ALIGNER=PairFinder"+
                            " TMPDIR=/broad/hptmp/andrewk")
            else:
                # This is a PAIRED LANE
                print "Processing paired lanes",p1,p2
                if p1.flowcell != p2.flowcell or p1.lane != p2.lane:
                    print "### WARNING ### paired lanes ",p1,p2,"do not match! Skipping."
                    continue
                # Filter 2nd Qltout file too
                cmd_str += " && FilterAlignsILTStyle QLTIN="+p2.qltout+" QLTOUT="+p2.fqltout
                cmd_str += (" && QltoutPair2SAM"+
                            " HEAD="+p1.head+
                            " SAM="+p1.sam+
                            " QLTOUT1="+p1.fqltout+" QLTOUT2="+p2.fqltout+
                            " FASTB1="+p1.fastb+" FASTB2="+p2.fastb+
                            " QUALB1="+p1.qualb+" QUALB2="+p2.qualb+
                            " REFDICT="+p1.refdict+
                            " METRICS="+p1.metrics+
                            " ALIGNER=PairFinder"+
                            " TMPDIR=/broad/hptmp/andrewk")
                
            # Make the BAM file
            cmd_str += " && /seq/dirseq/samtools/current/samtools import "+lane.bam_ref_list+" "+p1.sam+" "+p1.bam+"\""
            # Now remove one or both SAM files if all previous commands returned 0 for success
            cmd_str += " && rm -f "+p1.sam
            if p2 != None: cmd_str += " && rm -f "+p2.sam

            status = cmd(cmd_str, "long", output_head=p1.path_head, just_print_commands=just_print_commands)
            if status == 0:
                lanes_done += 1
            else:
                failed_list.append("Unkown failure")
            #if lanes_done + len(failed_list) > 0:
                #break

    failed = len(failed_list)
    lanes_total = lanes_done+failed
    percent_done = lanes_done / lanes_total if lanes_done != 0 else 0
    print "Lane conversion success: %d/%d (%.1f%%)" % (lanes_done, lanes_total, percent_done)

def merge_bams(tgtc_file, outfile, just_print_commands=False):
    input_lane_list = make_lane_list(tgtc_file)
    #lane_list = [lp for lp in lane_list if lp[0].flowcell =="30LLT"]
    missing_files = False
    lane_list = []
    for p1,p2 in input_lane_list:
        if "adapted" not in p1.qualb or (p2 != None and "adapted" not in p2.qualb):
            print "NOT ADAPTED", p1.qualb, p2.qualb
        if not os.path.exists(p1.bam):
            #print "Missing file:",p1.bam
            missing_files = True
        else:
            lane_list.append( (p1,p2) )

    #return
    if not missing_files:
        if not os.path.exists(outfile):
            print len(lane_list),"lanes / lane pairs to merge"
            if not os.path.exists(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))

            cmd("java -cp /home/radon01/depristo/dev/workspace/Libraries/bin"+
                " edu.mit.broad.picard.sam.MergeSamFiles"+
                " "+" ".join(["I="+p1.bam for p1,p2 in lane_list])+
                " O="+outfile+
                " TMP_DIR=/broad/hptmp/andrewk/"+
                " && /seq/dirseq/samtools/current/samtools index "+outfile+ 
                " && /xchip/tcga/gbm/visualization/sam/IGVAlignmentProcessor/bin/run_linux-64.sh "+outfile+" "+os.path.dirname(outfile)+" hg18 -w 100",
                "long", outfile, just_print_commands=just_print_commands)

            # Original samtools merge command below that was used for rev. 1 merging

            #cmd("/seq/dirseq/samtools/current/samtools merge "+outfile+
            #    " "+" ".join([p1.bam for p1,p2 in lane_list])+
            #    " && /seq/dirseq/samtools/current/samtools index "+outfile+
            #    " && /xchip/tcga/gbm/visualization/sam/IGVAlignmentProcessor/bin/run_linux-64.sh "+outfile+" "+os.path.dirname(outfile)+" hg18 -w 100",
            #    "long", outfile) 

    else:
        # Convert the lanes that had missing outfiles
        #convert_lanes(tgtc_file)
        pass

def process_tgtc_lists(tgtc_list_file):
    "Works through a list of tgtc files and merges them to produce the output"
    lists_done = 0
    
    for line in open(tgtc_list_file):
        if line[0] == '#':
            continue
        fields = line.split()
        if len(fields) != 2:
            continue
        head, tgtc_list = fields
        bam = head+".bam"
        if True: #not os.path.exists(bam):
            print "Processing tgtc_list_file",tgtc_list,"with goal of producing",bam
            merge_bams(tgtc_list, bam)
            print 
        else:
            print bam,"already exists\n"

        lists_done += 1
        #if lists_done > 0:
        #    break

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
    operation = sys.argv[1]

    # Check for print commands only option
    if "print_only" in sys.argv:
        just_print_commands = True
        sys.argv.remove("print_only")
    else:
        just_print_commands = False

    if operation == "convert":
        if len(sys.argv) < 3:
            usage()
        tgtc_file = sys.argv[2]
        # Run convert_lanes and recreate_all == True if clean option was given as 3rd arg
        convert_lanes(tgtc_file, len(sys.argv) > 3 and sys.argv[3] == "clean", just_print_commands=just_print_commands)
    elif operation == "merge":
        if len(sys.argv) < 4:
            usage()
        tgtc_file = sys.argv[2]
        outfile = sys.argv[3]
        merge_bams(tgtc_file, outfile, just_print_commands=just_print_commands)
    elif operation == "lists":
        tgtc_list_file = "/broad/1KG/legacy_data/tgtc2sam/tgtc_lists_2_legacy_bams"
        if len(sys.argv) >= 3:
            tgtc_list_file = sys.argv[2]
        process_tgtc_lists(tgtc_list_file)
    else:
        print "Didn't understand operation value: "+operation
        usage()
