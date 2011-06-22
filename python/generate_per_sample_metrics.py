#
# Reads in selected Picard metrics, generating an R-compatible TSV suitable for pre-QC analysis.
#
# To run:
#   /humgen/gsa-hpprojects/software/bin/jython2.5.2/jython \
#     -J-classpath $STING_HOME/dist/sam-1.47.869.jar:$STING_HOME/dist/picard-1.47.869.jar:$STING_HOME/dist/picard-private-parts-1941.jar \
#     $STING_HOME/python/generate_per_sample_metrics.py <bam.list> > <output_metrics_file.tsv>
#
# To add a new metric:
#   - If the metric file is new to Picard, add the relevant parser to the picard-private jar
#     (see http://www.broadinstitute.org/gsa/wiki/index.php/Adding_and_updating_dependencies for details).
#   - Add the field name to the header array.
#   - Add the field data to the statement printing the data array.
#
from java.lang import *
from java.io import File,FileReader
from net.sf.picard.metrics import MetricsFile

import os,string,sys

def median(l):
    return sorted(l)[(len(l)+1)/2]
def mean(l):
    return float(sum(l))/len(l)

def get_metrics(filename):
    if not os.path.exists(filename):
        return None
    file_reader = FileReader(filename)
    metrics_file = MetricsFile()
    metrics_file.read(file_reader)
    metrics = metrics_file.getMetrics()
    file_reader.close()
    return metrics

if len(sys.argv) != 2:
    print 'USAGE: %s <pipeline_file.yaml>'
    sys.exit(1)
if not os.path.exists(sys.argv[1]):
    print 'BAM list %s not found' % sys.argv[1]
    sys.exit(1)

bam_list_filename = sys.argv[1]

header = ['sample','HAPLOTYPES_CONFIDENTLY_MATCHING.MIN','HAPLOTYPES_CONFIDENTLY_MATCHING.MAX','HAPLOTYPES_CONFIDENTLY_MATCHING.MEDIAN',
          'BAIT_SET','GENOME_SIZE','PCT_SELECTED_BASES','MEAN_TARGET_COVERAGE','ZERO_CVG_TARGETS_PCT','FOLD_80_BASE_PENALTY','PCT_TARGET_BASES_2X',
          'PCT_TARGET_BASES_10X','PCT_TARGET_BASES_20X','PCT_TARGET_BASES_30X','HS_LIBRARY_SIZE','PCT_PF_READS_ALIGNED','PF_HQ_ERROR_RATE',
          'PF_INDEL_RATE','MEAN_READ_LENGTH','BAD_CYCLES','STRAND_BALANCE','PCT_CHIMERAS','PCT_ADAPTER','MEDIAN_INSERT_SIZE','TOTAL_SNPS']
data = ['%s'] * len(header)

print string.join(header,'\t')

# get a representative BAM file for each sample, to use as a base path.  Note that this assumes every sample corresponds to the same base path.
bam_list = open(bam_list_filename,'r')
samples = dict()

for bam_filename in bam_list:
    bam_filename = bam_filename.strip()
    if bam_filename == '':
        continue
    bam_filename_tokens = bam_filename.split('/')
    sample_id = bam_filename_tokens[len(bam_filename_tokens)-3]
    samples[sample_id] = bam_filename
bam_list.close()

for sample_id,filename in samples.items():
    basepath = filename[:filename.rindex('.bam')]
    
    fingerprinting_summary_metrics = get_metrics('%s.%s' % (basepath,'fingerprinting_summary_metrics'))

    if fingerprinting_summary_metrics != None:
        haplotypes_confidently_matching = [metric.HAPLOTYPES_CONFIDENTLY_MATCHING for metric in fingerprinting_summary_metrics]
        min_haplotypes_confidently_matching = str(min(haplotypes_confidently_matching))
        max_haplotypes_confidently_matching = str(max(haplotypes_confidently_matching))
        median_haplotypes_confidently_matching = str(median(haplotypes_confidently_matching))
    else:
        min_haplotypes_confidently_matching = 'NA'
        max_haplotypes_confidently_matching = 'NA'
        median_haplotypes_confidently_matching = 'NA'

    insert_size_metrics = get_metrics('%s.%s' % (basepath,'insert_size_metrics'))

    if insert_size_metrics != None:
        median_insert_size = insert_size_metrics[0].MEDIAN_INSERT_SIZE
    else:
        median_insert_size = 'NA'


    hybrid_selection_metrics = get_metrics('%s.%s' % (basepath,'hybrid_selection_metrics'))[0]
    alignment_summary_metrics = get_metrics('%s.%s' % (basepath,'alignment_summary_metrics'))[0]
    dbsnp_matches = get_metrics('%s.%s' % (basepath,'dbsnp_matches'))[0]

    print string.join(data,'\t')%(sample_id,min_haplotypes_confidently_matching,max_haplotypes_confidently_matching,median_haplotypes_confidently_matching,
                                  hybrid_selection_metrics.BAIT_SET,hybrid_selection_metrics.GENOME_SIZE,hybrid_selection_metrics.PCT_SELECTED_BASES,
                                  hybrid_selection_metrics.MEAN_TARGET_COVERAGE,hybrid_selection_metrics.ZERO_CVG_TARGETS_PCT,hybrid_selection_metrics.FOLD_80_BASE_PENALTY,
                                  hybrid_selection_metrics.PCT_TARGET_BASES_2X,hybrid_selection_metrics.PCT_TARGET_BASES_10X,hybrid_selection_metrics.PCT_TARGET_BASES_20X,
                                  hybrid_selection_metrics.PCT_TARGET_BASES_30X,hybrid_selection_metrics.HS_LIBRARY_SIZE,alignment_summary_metrics.PCT_PF_READS_ALIGNED,
                                  alignment_summary_metrics.PF_HQ_ERROR_RATE,alignment_summary_metrics.PF_INDEL_RATE,alignment_summary_metrics.MEAN_READ_LENGTH,
                                  alignment_summary_metrics.BAD_CYCLES,alignment_summary_metrics.STRAND_BALANCE,alignment_summary_metrics.PCT_CHIMERAS,
                                  alignment_summary_metrics.PCT_ADAPTER,median_insert_size,dbsnp_matches.TOTAL_SNPS)
