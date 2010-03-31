import os
import JobDispatcher

PRINT_ONLY = True
QUEUES = ["gsa","short","long"]
LIMITS = dict([["gsa",50],["short",20],["long",50]])
ACTION = "space"
TIME_DIFF = "0:1:0" # 1 hour
RESULTS_DIR = "/humgen/gsa-hpprojects/dev/kiran/wholeExomeGeneCoverage/scratch/"
GATK_JAR = "/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar"
MEMORY = "2g"
WALKER = "CoverageStatistics"
REFERENCE = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
BAM_FILES = RESULTS_DIR+"bam_files_use.list"
INTERVAL_LIST = "/seq/references/HybSelOligos/whole_exome_agilent_designed_120/whole_exome_agilent_designed_120.targets.interval_list"
MAX_BASES_PER_JOB = 100000
ARGUMENTS = "-mmq 20 -mbq 10 -omitBaseOutput -dels -l INFO"

dispatcher = JobDispatcher.GATKDispatcher(GATK_JAR,MEMORY,WALKER,ARGUMENTS,RESULTS_DIR,REFERENCE,
                                          BAM_FILES,INTERVAL_LIST,QUEUES,LIMITS,PRINT_ONLY,ACTION,TIME_DIFF)

dispatcher.setProject("GSA_WholeExome_CoverageByExon")
dispatcher.dispatchByInterval(MAX_BASES_PER_JOB)
