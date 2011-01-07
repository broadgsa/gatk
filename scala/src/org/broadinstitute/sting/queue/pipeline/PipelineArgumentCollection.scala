package org.broadinstitute.sting.queue.pipeline

import org.broadinstitute.sting.commandline.{Input, Argument, Hidden}
import java.io.File


class PipelineArgumentCollection {

  @Argument(doc="Number of cleaning jobs", shortName="cleaningJobs", required=false)
  var cleaningJobs: Int = 1

  @Argument(doc="the YAML file specifying inputs, interval lists, reference sequence, etc.", shortName="Y", required=false)
  var yamlFile: File = _

  @Input(doc="path to trigger track (for UnifiedGenotyper)", shortName="trigger", required=false)
  var trigger: File = _

  @Input(doc="path to refseqTable (for GenomicAnnotator)", shortName="refseqTable",required=false)
  var refseqTable: File = _

  @Input(doc="path to Picard FixMateInformation.jar.  See http://picard.sourceforge.net/ .", required=false)
  var picardFixMatesJar: File = new java.io.File("/seq/software/picard/current/bin/FixMateInformation.jar")

  @Input(doc="path to GATK jar", shortName="gatk")
  var gatkJar: File = _

  @Input(doc="target Ti/Tv ratio for recalibration", shortName="titv", required=true)
  var target_titv: Float = _

  @Input(doc="per-sample downsampling level",shortName="dcov",required=false)
  var downsampling_coverage = 300

  @Input(doc="level of parallelism for UnifiedGenotyper", shortName="snpScatter", required=false)
  var num_snp_scatter_jobs = 20

  @Input(doc="level of parallelism for IndelGenotyperV2", shortName="indelScatter", required=false)
  var num_indel_scatter_jobs = 5

  @Input(doc="Skip indel-cleaning for BAM files (for testing only)", shortName="skipCleaning", required=false)
  var skip_cleaning = false

  @Input(doc="List of samples and bams (in the form sample_id k1:v1,k2:v2 cleaned:/path/to/cleaned.bam,recalibrated:/path/to/recal.bam,unreacalibrated:/path/to/unrecal.bam). Mutually exclusive with YAML",
    required=false, shortName="pBams")
  var projectBams: File = _

  @Input(doc="The project name. Mutually exclusive with YAML.", required = false, shortName="pName")
  var projectName: String = _

  @Input(doc="The reference file. Mutually exclusive with YAML.", required=false, shortName="pRef")
  var projectRef: File = _

  @Input(doc="The project interval list. Mutually exclusive with YAML.", required=false, shortName="pInt")
  var projectIntervals: File = _

  @Input(doc="The project dbsnp. Mutually exclusive with YAML",required=false, shortName="pDB")
  var projectDBSNP: File = _

  //@Input(doc="ADPR script")
  //var adprScript: File = _

  //@Input(doc="Sequencing maching name (for use by adpr)")
  //var machine: String = _

  //@Input(doc="Sequencing experiement type (for use by adpr)--Whole_Exome, Whole_Genome, or Hybrid_Selection")
  //var protocol: String = _

  def verifyArguments = {
    // if no yaml file we require all the mutually exclusive arguments
    if ( yamlFile == null ) {
      if ( projectBams == null ) throw new IllegalArgumentException("No YAML file provided; and no project bam list. Pipeline requires either YAML file, or full project specs.")

      if ( projectName == null ) throw new IllegalArgumentException("No YAML file provided; and no project name. Pipeline requires either YAML file, or full project specs.")

      if ( projectRef == null) throw new IllegalArgumentException("No YAML file provided; and no project reference. Pipeline requires either YAML file, or full project specs.")

      if ( projectIntervals == null ) throw new IllegalArgumentException("No YAML file provided; and no project intervals. Pipeline requires either YAML file, or full project specs.")

      if ( projectDBSNP == null ) throw new IllegalArgumentException("No YAML file provided; and no project DBSNP. Pipeline requires either YAML file, or full project specs.")
    } else {
      if ( projectBams != null || projectName != null || projectRef != null || projectIntervals != null || projectDBSNP != null )
        throw new IllegalArgumentException("YAML file provided, along with other project spec arguments. YAML file is mutually exclusive with project specs.")
    }
  }

}