#!/usr/bin/env python

import sys
import os

bam_file = sys.argv[1]
fingerprint_file = sys.argv[2]
project = sys.argv[3]
directory = bam_file.rsplit("/",1)[0]+"/"
sample_id = bam_file.rsplit("/",1)[1].rsplit(".",1)[0]
is_metrics = directory+sample_id+".insert_size_metrics"
his_metrics = directory+sample_id+".hybrid_selection_metrics"
ali_metrics = directory+sample_id+".alignment_summary_metrics"
os.system("zip -j "+project+"_"+sample_id+"_sequencing_metrics"+" "+is_metrics+" "+his_metrics+" "+ali_metrics+" "+fingerprint_file)
