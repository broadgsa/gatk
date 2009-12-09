#!/usr/bin/env python
# imports
import subprocess
import os
import sys
import re
import time
import math
import codecs
## script keys and triggers
remove_all_files = False # be very careful with this one
make_directories = True
copy_vcf_files = True
annotate_vcf_files = True
run_snp_selector = False
make_hard_threshold_vcf = False
snp_select_threshold_vcf = False
generate_hapmap_info = False
generate_threshold_hapmap_info = False
false_negatives_by_snp_selector = False
false_negatives_by_snp_selector_on_threshold = False
plot_info_field_metrics = False
plot_info_field_metrics_threshold = False
make_best_effort_vcf = False
snp_select_best_effort_vcf = False
plot_best_effort_metrics = False


## global stuff
DEBUG = True
override_pilot_bamfile = True
override_production_bamfile = True
matlab_dir = "/humgen/gsa-scr1/projects/FHS/scripts"
pipeline_samples_info_path = "/humgen/gsa-hphome1/flannick/pfizer/pspipeline/meta/samples.tsv"
samples_pipeline_path = "/humgen/gsa-hphome1/flannick/pfizer/pspipeline/output/samples/"
home_dir = "/humgen/gsa-scr1/projects/FHS/"
annotations_data_base = home_dir + "oneOffAnalyses/annotationEffectiveness/data/"
figures_base = home_dir + "oneOffAnalyses/annotationEffectiveness/figures/"
production_calls_dir = home_dir + "production/raw_production_calls_12_07/"
pilot_calls_dir = home_dir + "pilot/calls/raw_pilot_calls_12_07/"
hapmap_pool_info_dir = home_dir+"pilot/project_info/pilot_pool_information/"
pilot_standard_project_name = "framinghampilot"
pilot_clipped_project_name = "pilot_clipped"
pilot_no_gatk_project_name = "pilot_no_gatk"
FHS_project_name = "FHS"
FHS_round_2_project_name = "FHS_round_2"
all_projects = [pilot_standard_project_name,pilot_clipped_project_name,pilot_no_gatk_project_name,FHS_project_name,FHS_round_2_project_name]
pilot_projects = [pilot_standard_project_name,pilot_clipped_project_name,pilot_no_gatk_project_name]
production_projects = [FHS_project_name,FHS_round_2_project_name]

gatk_cmd_base = "java -jar /humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -L /humgen/gsa-scr1/projects/FHS/interval_lists/FHS_exons_only.interval_list"
dbsnp_129_path = "/humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod"

vcf_header = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]

def subSystem(cmd):
    if ( DEBUG ):
        print(cmd)
    else:
        os.system(cmd)

def safeRemove(cmd):
    if ( not ( cmd.startswith("rm -rf "+annotations_data_base) or cmd.startswith("rm -rf "+production_calls_dir) or cmd.startswith("rm -rf "+pilot_calls_dir) or cmd.startswith("rm -rf "+figures_base) ) ):
        print("Unsafe removal attempt: "+cmd)
    else:
        subSystem(cmd)

def productionBamfilePicard(proj, samp):
    for line in open("/humgen/gsa-hphome1/flannick/pfizer/pspipeline/meta/framingham.flowcell_lanes.tsv").readlines():
        spline = line.strip().split()
        if ( spline[1] == proj and spline[2] == samp ):
            return spline[5]+spline[3]+"."+spline[4]+".aligned.bam"
    print("no return for: "+proj+"_"+samp)

def annotationProjectPath(proj):
    return annotations_data_base+proj+"/"

def annotationPath(proj,pool):
    return annotationProjectPath(proj)+pool+"/"

def annotationProjectFile(proj, pool, appender):
    return annotationPath(proj,pool)+pool+appender

def figureProjectPath(proj):
    return figures_base+proj+"/"

def figurePath(proj, samp):
    return figureProjectPath(proj)+samp+"/"

def pipelineProjectPath(proj):
    return samples_pipeline_path+proj+"/"

def pipelinePath(proj, samp):
    return pipelineProjectPath(proj)+samp+"/"

def pipelineBam(proj, samp):
    if ( proj == "framinghampilot" and override_pilot_bamfile ):
        if ( samp == "CEPH1" ):
            return "/seq/picard/302JDAAXX/C1-252_2009-04-30_2009-11-13/1/Solexa-10930/302JDAAXX.1.aligned.bam -I /seq/picard/302JDAAXX/C1-252_2009-04-30_2009-11-13/5/Solexa-10933/302JDAAXX.5.aligned.bam"
        if ( samp == "CEPH2" ):
            return "/seq/picard/302JDAAXX/C1-252_2009-04-30_2009-11-13/2/Solexa-10931/302JDAAXX.2.aligned.bam -I /seq/picard/302JDAAXX/C1-252_2009-04-30_2009-11-13/8/Solexa-10934/302JDAAXX.8.aligned.bam"
        if ( samp == "CEPH3" ):
            return "/seq/picard/302JDAAXX/C1-252_2009-04-30_2009-11-13/3/Solexa-10932/302JDAAXX.3.aligned.bam -I /seq/picard/302JDAAXX/C1-252_2009-04-30_2009-11-13/6/Solexa-10936/302JDAAXX.6.aligned.bam -I /seq/picard/302JDAAXX/C1-252_2009-04-30_2009-11-13/7/Solexa-10935/302JDAAXX.7.aligned.bam"

    if ( proj in production_projects and override_production_bamfile ):
        return productionBamfilePicard(proj,samp)

    return pipelinePath(proj,samp)+proj+"."+samp+".bam"

def pipelineVCF(proj,samp):
    return pipelinePath(proj,samp)+samp+".vcf"

def homePath(proj, samp):
    if ( proj in pilot_projects ):
        return pilot_calls_dir+proj+"/"
    else:
        return production_calls_dir+proj+"/"

def homeProjectFile(proj, samp, appender):
    return homePath(proj,samp)+samp+appender

def homeVCF(proj, samp, extension):
    return homeProjectFile(proj,samp,extension)+".vcf"

def hapmapVCF(proj, samp):
    return "/humgen/gsa-scr1/projects/FHS/pilot/analysis/"+samp.lower()+"_hapmap_snp_only.vcf"

def poolSize(proj):
    if ( proj in pilot_projects ):
        return 40
    else:
        return 28

def poolBinding(proj, samp):
    if ( proj in production_projects ):
        raise Exception, "Production projects do not have hapmap pool bindings"
    else:
        return hapmap_pool_info_dir+samp+".pool.path"

def poolSampleNames(proj,samp):
    if ( proj in production_projects ):
        raise Exception, "Production projects do not have hapmap pool bindings"
    else:
        return hapmap_pool_info_dir+samp+".pool"


## working stuff

#working_projects = [pilot_standard_project_name]
working_projects = [FHS_round_2_project_name]
working_info_fields = ["syzy_SB","syzy_DP","syzy_NMMR","HRun"]
working_quality_scores = ["0","1","2","3","4","5"]
working_samples = set()
projects_to_samples = []

# create the project --> list of samples dictionary

for project in working_projects:
    sample_list = os.listdir(pipelineProjectPath(project))
    projects_to_samples.append([project, sample_list])
    for s in sample_list:
        working_samples.add(s)

projects_to_samples = dict(projects_to_samples)

# create the sample --> pool name dictionary

samples_to_pool_name = []

for line in open(pipeline_samples_info_path).readlines():
    spline = line.strip().split()
    if(spline[0] in working_samples):
        keyval = [spline[0], spline[1]]
        samples_to_pool_name.append(keyval)

samples_to_pool_name = dict(samples_to_pool_name)

# remove everything
if ( remove_all_files ):
    for project in working_projects:
        safeRemove("rm -rf "+homePath(project,"") )
        safeRemove("rm -rf "+annotationProjectPath(project))
        safeRemove("rm -rf "+figureProjectPath(project))

# make directories if necessary
if ( make_directories ):
    for project in working_projects:
        subSystem("mkdir "+homePath(project,""))
        subSystem("mkdir "+annotationProjectPath(project))
        subSystem("mkdir "+figureProjectPath(project))
        for sample in projects_to_samples[project]:
            subSystem("mkdir "+annotationPath(project, samples_to_pool_name[sample]))
            subSystem("mkdir "+figurePath(project,samples_to_pool_name[sample]))

# copy the vcf files from the pipeline into the /FHS/ directory
if ( copy_vcf_files ):
    for project in working_projects:
        sample_list = projects_to_samples[project]
        for sample in sample_list:
            pool = samples_to_pool_name[sample]
            sample_home_vcf = homeVCF(project,samples_to_pool_name[sample],"_raw")
            pipeline_vcf = pipelineVCF(project,sample)
            subSystem("cat "+pipeline_vcf+" | sed 's/"+sample+"/"+pool+"/g' > "+sample_home_vcf)

# annotate the vcf files using VariantAnnotator
if ( annotate_vcf_files ):
    prior_step_extension = "_raw"
    this_step_extension = "_annotated"
    gatk_annotate_base = gatk_cmd_base + " -T VariantAnnotator -exp -D "+dbsnp_129_path
    for project in working_projects:
        sample_list = projects_to_samples[project]
        for sample in sample_list:
            pool = samples_to_pool_name[sample]
            inputVCF = homeVCF(project, pool, prior_step_extension)
            outputVCF = homeVCF(project, pool, this_step_extension)
            bamfile = pipelineBam(project, sample)
            gatk_args = " -I "+bamfile+" -B variant,VCF,"+inputVCF+" -vcf "+outputVCF
            subSystem("bsub -q gsa "+gatk_annotate_base + gatk_args)

# definition for next section (and another one further down)

def runSnpSelector(project_list,prior_step_annotation,this_step_annotation,info_field_list):
    snp_selector_base = "python /humgen/gsa-scr1/chartl/sting/python/snpSelector.py -p 10 --plottable"
    for project in project_list:
        sample_list = projects_to_samples[project]
        for sample in sample_list:
            pool = samples_to_pool_name[sample]
            if ( project in pilot_projects ):
                truth_arg = " -t "+hapmapVCF(project,sample)
            else:
                truth_arg = " -titv=3.6"
            for info_field in info_field_list:
                input_vcf = homeVCF(project,pool,prior_step_extension)
                output_vcf = homeVCF(project,pool,this_step_annotation+"_"+info_field)
                log_output = annotationProjectFile(project,pool,this_step_annotation+"_"+info_field+".log")
                snp_selector_args = truth_arg+" -o "+output_vcf + " -l "+log_output + " -f "+info_field+" "+input_vcf
                subSystem(snp_selector_base+snp_selector_args)

# run snp selector to analyze info fields
if ( run_snp_selector ):
    prior_step_extension = "_annotated"
    this_step_extension = "_snpsel"
    runSnpSelector(working_projects,prior_step_extension, this_step_extension, working_info_fields)

# definition for the next stage

def parseInfoField(info):
    info_list = info.split(";")
    info_dict = []
    for info_field in info_list:
        keyval = info_field.split("=")
        key = keyval[0]
        val = float(keyval[1])
        info_dict.append([key, val])
    return dict(info_dict)


def vcfLinePassesThreshold(line, fields, vals,greaterthan):
    if( line.startswith('#') ):
        return True
    else:
        # parse the line
        general_fields = line.split()
        info_dict = parseInfoField(general_fields[vcf_header.index("INFO")])
        header_set = set(vcf_header)
        # compare with thresholds
        for j in range(len(fields)):
            if ( fields[j] in header_set ):
                if ( float(general_fields[vcf_header.index(fields[j])]) > vals[j] ):
                    if ( not greaterthan[j] ):
                        return False
                    else:
                        continue
                else:
                    if ( greaterthan[j] ):
                        return False
                    else:
                        continue
            elif ( fields[j] in info_dict ):
                if ( float(info_dict[fields[j]]) > vals[j] ):
                    if ( not greaterthan[j] ):
                        return False
                    else:
                        continue
                else:
                    if ( greaterthan[j] ):
                        return False
                    else:
                        continue
            else:
                print("Field not found "+fields[j])
        return True

def makeThresholdVCF(prev_vcf,this_vcf,fields,vals,greater):
    if ( DEBUG ):
        filtered_lines = 0
    else:
        out_vcf = open(this_vcf,'w')
        print("open for writing: "+this_vcf)
        wrotelines = 0
        for line in open(prev_vcf).readlines():
            if ( vcfLinePassesThreshold(line,fields,vals,greater) ):
                if ( DEBUG ):
                    filtered_lines = filtered_lines + 1
                else:                    
                    out_vcf.write(line)
                    wrotelines=wrotelines+1
        if ( DEBUG ):
            print(this_vcf+" would have filtered "+str(filtered_lines)+" variants")
        else:
            out_vcf.close()

# make hard-threshold vcf files ( to see if a more strenuous set with less noise helps snp selector determine an ordering )
if ( make_hard_threshold_vcf ):
    prior_step_extension = "_annotated"
    this_step_extension = "_hard_threshold"
    threshold_fields = ["QUAL", "syzy_NMMR","syzy_DP"]
    threshold_values = [15, 0.4,280]
    threshold_greaterthan = [True, False,True]
    for project in working_projects:
        sample_list = projects_to_samples[project]
        for sample in sample_list:
            pool = samples_to_pool_name[sample]
            prev_vcf = homeVCF(project,pool,prior_step_extension)
            this_vcf = homeVCF(project,pool,this_step_extension)
            makeThresholdVCF(prev_vcf,this_vcf,threshold_fields,threshold_values,greaterthan)

# run snp selector on hard-threshold vcf files
if ( snp_select_threshold_vcf ):
    prior_step_extension = "_hard_threshold"
    this_step_extension = "_hard_snpsel"
    runSnpSelector(working_projects, prior_step_extension, this_step_extension, working_info_fields)

# definitions for next few steps

def hapmapInfoFile(proj,pool,info,qscore):
    return annotationProjectFile(project,pool,"_"+info+"_filter_at_q_"+qscore+"_hapmap_info.txt")

def generateHapmapInfo(extension, qscorelist):
    gatk_hapmap_info = gatk_cmd_base+" -T HapmapPoolAllelicInfo"
    for project in working_projects:
        sample_list = projects_to_samples[project]
        for sample in sample_list:
            pool = samples_to_pool_name[sample]
            for infofield in working_info_fields:
                inputVCF = homeVCF(project,pool,prior_step_extension+"_"+infofield)
                inputBam = pipelineBam(project, sample)
                ps = str(poolSize(project))
                for qscore in qscorelist:
                    outputFile = hapmapInfoFile(project, pool, infofield, qscore)
                    gatk_args = " -ps "+ps+" -I "+inputBam+" -of "+outputFile+" -B "+pool+",VCF,"+inputVCF+" -B "+poolBindings(project,sample)+" -samples "+poolSampleNames(project,sample)+" -q "+qscore
                    subSystem("bsub -q gsa "+gatk_hapmap_info+gatk_args)

# generate hapmap (false postive/false negative) info files on filtered vcfs
if ( generate_hapmap_info ):
    prior_step_extension = "_snpsel"
    generateHapmapInfo(prior_step_extension, working_quality_scores)

if ( generate_threshold_hapmap_info ):
    prior_step_extension = "_hard_snpsel"
    generateHapmapInfo(prior_step_extension, working_quality_scores)

# definitions for next section
def falsenegSnpSelect(invcf, false_neg_vcf, qscore):
    tmp_vcf = "tmp.vcf"
    tmp_recal_vcf = "tmp2.vcf"
    tmp_log = "tmp.log"
    snpsel_base = "python /humgen/gsa-scr1/chartl/sting/python/snpSelector.py -p 30 --plottable -l "+tmp_log+" -o "+tmp_recal_vcf+" --FNOutputVCF="+false_neg_vcf+" "+tmp_vcf
    cmd = "cat "+invcf+" | awk '{if ( substr($1,0,1) == '#' || $5 > "+qscore+" ) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' | sed 's/ /\t/g' > "+tmp.vcf
    subSystem(cmd)
    subSystem(snpsel_base)
    cmd = "rm "+tmp_vcf
    subSystem(cmd)
    cmd = "rm "+tmp_recal_vcf
    subSystem(cmd)
    cmd = "rm "+tmp_log
    subSystem(cmd)

def falseNegsBySnpSelector(extension):
    for project in working_projects:
        sample_list = projects_to_samples[project]
        for sample in sample_list:
            pool = samples_to_pool_name[sample]
            for infofield in working_info_fields:
                input_vcf = homeVCF(project,pool,extension+"_"+infofield)
                for qscore in working_quality_scores:
                    false_neg_vcf = annotationProjectFile(project, pool, "_"+infofield+"_q_"+qscore+"_false_negs.vcf")
                    falsenegSnpSelect(input_vcf,false_neg_vcf,qscore)

# generate false negative sites in a vcf file via snp selector
if ( false_negatives_by_snp_selector):
    prior_step_extension = "_snpsel"
    falseNegsBySnpSelector(prior_step_extension)

# generate false negative sites ina  vcf file on the thresholded vcfs
if ( false_negatives_by_snp_selector_on_threshold ):
    prior_step_extension = "_hard_snpsel"
    falseNegsBySnpSelector(prior_step_extension)

# definition for matlab section
def makeMatlabMetricPlot(extension):
        for project in working_projects:
            sample_list = projects_to_samples[project]
            for sample in sample_list:
                pool = samples_to_pool_name[sample]
                for infofield in working_info_fields:
                    log_output = annotationProjectFile(project,pool,extension+"_"+infofield+".log")
                    figureName = figurePath(project,pool)+pool+extension
                    cmd = "matlab -r \"cd "+matlab_dir+" , snpSelectorPlots('"+log_output+"','"+figureName+"') , exit\""
                    subSystem(cmd)

# make matlab plots of the metrics (generated from the log files)
if ( plot_info_field_metrics ):
    prior_step_extension = '_snpsel'
    makeMatlabMetricPlot(prior_step_extension)

# make matlab plots of the metrics from the hard-thresholded log files
if ( plot_info_field_metrics_threshold ):
    prior_step_extension = "_hard_snpsel"
    makeMatlabMetricPlot(prior_step_extension)

# definitions for best effort section

def bestEffortVCF(proj):
    bevcf_path = homeVCF(proj,proj,"_best_effort")
    #bevcf = open(bevcf_path,'w')
    #bevcf.write("##source=best_effort_vcf")
    #bevcf.write("##format=VCRv3.2")
    #bevcf.write("#"+"\t".join(vcf_header)+"\t"+proj+"_combined_calls")
    #bevcf.close()
    return bevcf_path

def keyFromFields(fields):
    return ":".join([fields[0],fields[1]])

def inVCFDict(vdict, fields):
    return ( keyFromFields(fields) in vdict )

def addToVCFDict(vdict, fields):
    vdict[keyFromFields(fields)] = fields
    return vdict

def updateInfoField(key, info1, info2):
    if ( key == "HRun" ): # this one doesn't change
        return [key, str(info1)]
    else:
        return [key, str( (float(info1) + float(info2))/2 )]

def incorporateNewInfo(vdict, infodict, fields):
    addNPools = True
    poskey = keyFromFields(fields)
    vline = vdict[poskey]
    vinfodict = parseInfoField(vline[vcf_header.index("INFO")])
    newInfo = []
    all_keys = set(vinfodict.keys()).union(set(infodict.keys()))
    for infokey in all_keys:
        if ( infokey in vinfodict and infokey in infodict ):
            newInfo.append("=".join(updateInfoField(infokey, vinfodict[infokey], infodict[infokey]) ) )
        elif ( infokey in infodict ):
            newInfo.append("=".join([infokey, str(infodict[infokey])]) )
        elif ( infokey == "nPools" ):
            newPools = str( int(vinfodict[infokey]) + 1 )
            newInfo.append("=".join([infokey,newPools]))
            addNPools = False
        else:
            continue
    if ( addNPools ):
        newInfo.append("nPools=1")
    vline[vcf_header.index("INFO")] = ";".join(newInfo)
    vdict[poskey] = vline
    return vdict

def incorporateQuality(vdict, fields):
    poskey = keyFromFields(fields)
    vfields = vdict[poskey]
    qual_index = vcf_header.index("QUAL")
    vfields[qual_index] = str( float(vfields[qual_index]) + float(fields[qual_index]) )
    vdict[poskey] = vfields
    return vdict

def bestEffortMerge(call_dict,vcf_to_merge_path):
    for line in open(vcf_to_merge_path).readlines():
        if ( line.startswith("#") ):
            continue
        else:
            vcf_fields = line.strip().split()
            if ( not inVCFDict(call_dict, vcf_fields) ):
                call_dict = addToVCFDict(call_dict,vcf_fields)
            else:
                vcf_info_dict = parseInfoField(vcf_fields[vcf_header.index("INFO")])
                call_dict = incorporateNewInfo(call_dict,vcf_info_dict,vcf_fields)
                call_dict = incorporateQuality(call_dict,vcf_fields)
    return call_dict

def dictToFile(vdict, filepath,source,samplename):
    file = open(filepath,'w')
    file.write("##source="+source)
    file.write("##version=VCRv3.2")
    file.write("#"+"\t".join(vcf_header)+"\t"+samplename)
    if ( vdict ):
        for key in vdict.keys():
            file.write("\t".join(vdict[key])+"\n")
        file.close()
    else:
        print("dict had no keys for "+filepath)

# make the best effort VCF
if ( make_best_effort_vcf ):
    best_thresh_fields = ["QUAL"]
    best_thresh_val = [7]
    best_thresh_gt = [True]
    prev_extension = "_annotated"
    for project in working_projects:
        best_effort_vcf = bestEffortVCF(project)
        best_effort_dict = dict()
        sample_list = projects_to_samples[project]
        for sample in sample_list:
            pool = samples_to_pool_name[sample]
            in_vcf = homeVCF(project,pool,prev_extension)
            outvcf = homeVCF(project,pool,"_temporary")
            makeThresholdVCF(in_vcf,outvcf,best_thresh_fields,best_thresh_val,best_thresh_gt)
            print("dict about to update, size is currently: "+str(len(best_effort_dict)))
            best_effort_dict = bestEffortMerge(best_effort_dict,outvcf)
            print("updated size: "+str(len(best_effort_dict)))
            subSystem("rm "+outvcf)
        dictToFile(best_effort_dict,best_effort_vcf,"make_best_effort_vcf",project+"_best_effort")

if ( snp_select_best_effort_vcf ):
    snp_selector_base = "python /humgen/gsa-scr1/chartl/sting/python/snpSelector.py -p 20 --plottable -f "+",".join(working_info_fields)
    for project in working_project:
        best_effort_vcf = bestEffortVCF(project)
        snp_recal_vcf = homeVCF(project,project,"_best_effort_selected")
        log = annotationProjectFile(project,project,"_best_effort_selected.log")
