/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.extensions.cancer

import java.io.File
import org.broadinstitute.gatk.utils.commandline.Argument
import org.broadinstitute.gatk.utils.commandline.Gather
import org.broadinstitute.gatk.utils.commandline.Input
import org.broadinstitute.gatk.utils.commandline.Output
import org.broadinstitute.gatk.queue.function.scattergather.ScatterGatherableFunction
import org.broadinstitute.gatk.queue.extensions.gatk.{TaggedFile, VcfGatherFunction, LocusScatterFunction}
import org.broadinstitute.gatk.utils.commandline.ArgumentTypeDescriptor.isCompressed

class MuTect extends org.broadinstitute.gatk.queue.extensions.gatk.CommandLineGATK with ScatterGatherableFunction {
  analysisName = "MuTect"
  analysis_type = "MuTect"
  scatterClass = classOf[LocusScatterFunction]

  /** used for debugging, basically exit as soon as we get the reads */
  @Argument(fullName="noop", shortName="", doc="used for debugging, basically exit as soon as we get the reads", required=false, exclusiveOf="", validation="")
  var noop: Boolean = _

  /** add many additional columns of statistics to the output file */
  @Argument(fullName="enable_extended_output", shortName="", doc="add many additional columns of statistics to the output file", required=false, exclusiveOf="", validation="")
  var enable_extended_output: Boolean = _

  /** used when running the caller on a normal (as if it were a tumor) to detect artifacts */
  @Argument(fullName="artifact_detection_mode", shortName="", doc="used when running the caller on a normal (as if it were a tumor) to detect artifacts", required=false, exclusiveOf="", validation="")
  var artifact_detection_mode: Boolean = _

  /** name to use for tumor in output files */
  @Argument(fullName="tumor_sample_name", shortName="", doc="name to use for tumor in output files", required=false, exclusiveOf="", validation="")
  var tumor_sample_name: String = _

  /** if the tumor bam contains multiple samples, only use read groups with SM equal to this value */
  @Argument(fullName="bam_tumor_sample_name", shortName="", doc="if the tumor bam contains multiple samples, only use read groups with SM equal to this value", required=false, exclusiveOf="", validation="")
  var bam_tumor_sample_name: String = _

  /** name to use for normal in output files */
  @Argument(fullName="normal_sample_name", shortName="", doc="name to use for normal in output files", required=false, exclusiveOf="", validation="")
  var normal_sample_name: String = _

  /** force output for each site */
  @Argument(fullName="force_output", shortName="", doc="force output for each site", required=false, exclusiveOf="", validation="")
  var force_output: Boolean = _

  /** force output for all alleles at each site */
  @Argument(fullName="force_alleles", shortName="", doc="force output for all alleles at each site", required=false, exclusiveOf="", validation="")
  var force_alleles: Boolean = _

  /** only emit passing calls */
  @Argument(fullName="only_passing_calls", shortName="", doc="only emit passing calls", required=false, exclusiveOf="", validation="")
  var only_passing_calls: Boolean = _

  /** Initial LOD threshold for calling tumor variant */
  @Argument(fullName="initial_tumor_lod", shortName="", doc="Initial LOD threshold for calling tumor variant", required=false, exclusiveOf="", validation="")
  var initial_tumor_lod: Option[Float] = None

  /** Format string for initial_tumor_lod */
  @Argument(fullName="initial_tumor_lodFormat", shortName="", doc="Format string for initial_tumor_lod", required=false, exclusiveOf="", validation="")
  var initial_tumor_lodFormat: String = "%s"

  /** LOD threshold for calling tumor variant */
  @Argument(fullName="tumor_lod", shortName="", doc="LOD threshold for calling tumor variant", required=false, exclusiveOf="", validation="")
  var tumor_lod: Option[Float] = None

  /** Format string for tumor_lod */
  @Argument(fullName="tumor_lodFormat", shortName="", doc="Format string for tumor_lod", required=false, exclusiveOf="", validation="")
  var tumor_lodFormat: String = "%s"

  /** estimate of fraction (0-1) of physical contamination with other unrelated samples */
  @Argument(fullName="fraction_contamination", shortName="", doc="estimate of fraction (0-1) of physical contamination with other unrelated samples", required=false, exclusiveOf="", validation="")
  var fraction_contamination: Option[Float] = None

  /** Format string for fraction_contamination */
  @Argument(fullName="fraction_contaminationFormat", shortName="", doc="Format string for fraction_contamination", required=false, exclusiveOf="", validation="")
  var fraction_contaminationFormat: String = "%s"

  /** minimum fraction of cells which are presumed to have a mutation, used to handle non-clonality and contamination */
  @Argument(fullName="minimum_mutation_cell_fraction", shortName="", doc="minimum fraction of cells which are presumed to have a mutation, used to handle non-clonality and contamination", required=false, exclusiveOf="", validation="")
  var minimum_mutation_cell_fraction: Option[Float] = None

  /** Format string for minimum_mutation_cell_fraction */
  @Argument(fullName="minimum_mutation_cell_fractionFormat", shortName="", doc="Format string for minimum_mutation_cell_fraction", required=false, exclusiveOf="", validation="")
  var minimum_mutation_cell_fractionFormat: String = "%s"

  /** LOD threshold for calling normal non-germline */
  @Argument(fullName="normal_lod", shortName="", doc="LOD threshold for calling normal non-germline", required=false, exclusiveOf="", validation="")
  var normal_lod: Option[Float] = None

  /** Format string for normal_lod */
  @Argument(fullName="normal_lodFormat", shortName="", doc="Format string for normal_lod", required=false, exclusiveOf="", validation="")
  var normal_lodFormat: String = "%s"

  /** LOD threshold for calling normal non-variant */
  @Argument(fullName="normal_artifact_lod", shortName="", doc="LOD threshold for calling normal non-variant", required=false, exclusiveOf="", validation="")
  var normal_artifact_lod: Option[Float] = None

  /** Format string for normal_artifact_lod */
  @Argument(fullName="normal_artifact_lodFormat", shortName="", doc="Format string for normal_artifact_lod", required=false, exclusiveOf="", validation="")
  var normal_artifact_lodFormat: String = "%s"

  /** LOD threshold for calling strand bias */
  @Argument(fullName="strand_artifact_lod", shortName="", doc="LOD threshold for calling strand bias", required=false, exclusiveOf="", validation="")
  var strand_artifact_lod: Option[Float] = None

  /** Format string for strand_artifact_lod */
  @Argument(fullName="strand_artifact_lodFormat", shortName="", doc="Format string for strand_artifact_lod", required=false, exclusiveOf="", validation="")
  var strand_artifact_lodFormat: String = "%s"

  /** power threshold for calling strand bias */
  @Argument(fullName="strand_artifact_power_threshold", shortName="", doc="power threshold for calling strand bias", required=false, exclusiveOf="", validation="")
  var strand_artifact_power_threshold: Option[Float] = None

  /** Format string for strand_artifact_power_threshold */
  @Argument(fullName="strand_artifact_power_thresholdFormat", shortName="", doc="Format string for strand_artifact_power_threshold", required=false, exclusiveOf="", validation="")
  var strand_artifact_power_thresholdFormat: String = "%s"

  /** LOD threshold for calling normal non-variant at dbsnp sites */
  @Argument(fullName="dbsnp_normal_lod", shortName="", doc="LOD threshold for calling normal non-variant at dbsnp sites", required=false, exclusiveOf="", validation="")
  var dbsnp_normal_lod: Option[Float] = None

  /** Format string for dbsnp_normal_lod */
  @Argument(fullName="dbsnp_normal_lodFormat", shortName="", doc="Format string for dbsnp_normal_lod", required=false, exclusiveOf="", validation="")
  var dbsnp_normal_lodFormat: String = "%s"

  /** Power threshold for normal to determine germline vs variant */
  @Argument(fullName="somatic_classification_normal_power_threshold", shortName="", doc="Power threshold for normal to determine germline vs variant", required=false, exclusiveOf="", validation="")
  var somatic_classification_normal_power_threshold: Option[Float] = None

  /** Format string for somatic_classification_normal_power_threshold */
  @Argument(fullName="somatic_classification_normal_power_thresholdFormat", shortName="", doc="Format string for somatic_classification_normal_power_threshold", required=false, exclusiveOf="", validation="")
  var somatic_classification_normal_power_thresholdFormat: String = "%s"

  /** minimum allele fraction to be considered in normal, useful for normal sample contaminated with tumor */
  @Argument(fullName="minimum_normal_allele_fraction", shortName="", doc="minimum allele fraction to be considered in normal, useful for normal sample contaminated with tumor", required=false, exclusiveOf="", validation="")
  var minimum_normal_allele_fraction: Option[Float] = None

  /** Format string for minimum_normal_allele_fraction */
  @Argument(fullName="minimum_normal_allele_fractionFormat", shortName="", doc="Format string for minimum_normal_allele_fraction", required=false, exclusiveOf="", validation="")
  var minimum_normal_allele_fractionFormat: String = "%s"

  /** for computational efficiency, reject sites with allelic fraction below this threshold */
  @Argument(fullName="tumor_f_pretest", shortName="", doc="for computational efficiency, reject sites with allelic fraction below this threshold", required=false, exclusiveOf="", validation="")
  var tumor_f_pretest: Option[Float] = None

  /** Format string for tumor_f_pretest */
  @Argument(fullName="tumor_f_pretestFormat", shortName="", doc="Format string for tumor_f_pretest", required=false, exclusiveOf="", validation="")
  var tumor_f_pretestFormat: String = "%s"

  /** threshold for minimum base quality score */
  @Argument(fullName="min_qscore", shortName="", doc="threshold for minimum base quality score", required=false, exclusiveOf="", validation="")
  var min_qscore: Option[Int] = None

  /** how many gapped events (ins/del) are allowed in proximity to this candidate */
  @Argument(fullName="gap_events_threshold", shortName="", doc="how many gapped events (ins/del) are allowed in proximity to this candidate", required=false, exclusiveOf="", validation="")
  var gap_events_threshold: Option[Int] = None

  /** if this fraction or more of the bases in a read are soft/hard clipped, do not use this read for mutation calling */
  @Argument(fullName="heavily_clipped_read_fraction", shortName="", doc="if this fraction or more of the bases in a read are soft/hard clipped, do not use this read for mutation calling", required=false, exclusiveOf="", validation="")
  var heavily_clipped_read_fraction: Option[Float] = None

  /** Format string for heavily_clipped_read_fraction */
  @Argument(fullName="heavily_clipped_read_fractionFormat", shortName="", doc="Format string for heavily_clipped_read_fraction", required=false, exclusiveOf="", validation="")
  var heavily_clipped_read_fractionFormat: String = "%s"

  /** pvalue threshold for fishers exact test of clipping bias in mutant reads vs ref reads */
  @Argument(fullName="clipping_bias_pvalue_threshold", shortName="", doc="pvalue threshold for fishers exact test of clipping bias in mutant reads vs ref reads", required=false, exclusiveOf="", validation="")
  var clipping_bias_pvalue_threshold: Option[Float] = None

  /** Format string for clipping_bias_pvalue_threshold */
  @Argument(fullName="clipping_bias_pvalue_thresholdFormat", shortName="", doc="Format string for clipping_bias_pvalue_threshold", required=false, exclusiveOf="", validation="")
  var clipping_bias_pvalue_thresholdFormat: String = "%s"

  /** threshold for determining if there is relatedness between the alt and ref allele read piles */
  @Argument(fullName="fraction_mapq0_threshold", shortName="", doc="threshold for determining if there is relatedness between the alt and ref allele read piles", required=false, exclusiveOf="", validation="")
  var fraction_mapq0_threshold: Option[Float] = None

  /** Format string for fraction_mapq0_threshold */
  @Argument(fullName="fraction_mapq0_thresholdFormat", shortName="", doc="Format string for fraction_mapq0_threshold", required=false, exclusiveOf="", validation="")
  var fraction_mapq0_thresholdFormat: String = "%s"

  /** threshold for clustered read position artifact median */
  @Argument(fullName="pir_median_threshold", shortName="", doc="threshold for clustered read position artifact median", required=false, exclusiveOf="", validation="")
  var pir_median_threshold: Option[Double] = None

  /** Format string for pir_median_threshold */
  @Argument(fullName="pir_median_thresholdFormat", shortName="", doc="Format string for pir_median_threshold", required=false, exclusiveOf="", validation="")
  var pir_median_thresholdFormat: String = "%s"

  /** threshold for clustered read position artifact MAD */
  @Argument(fullName="pir_mad_threshold", shortName="", doc="threshold for clustered read position artifact MAD", required=false, exclusiveOf="", validation="")
  var pir_mad_threshold: Option[Double] = None

  /** Format string for pir_mad_threshold */
  @Argument(fullName="pir_mad_thresholdFormat", shortName="", doc="Format string for pir_mad_threshold", required=false, exclusiveOf="", validation="")
  var pir_mad_thresholdFormat: String = "%s"

  /** required minimum value for tumor alt allele maximum mapping quality score */
  @Argument(fullName="required_maximum_alt_allele_mapping_quality_score", shortName="", doc="required minimum value for tumor alt allele maximum mapping quality score", required=false, exclusiveOf="", validation="")
  var required_maximum_alt_allele_mapping_quality_score: Option[Int] = None

  /** threshold for maximum alternate allele counts in normal */
  @Argument(fullName="max_alt_alleles_in_normal_count", shortName="", doc="threshold for maximum alternate allele counts in normal", required=false, exclusiveOf="", validation="")
  var max_alt_alleles_in_normal_count: Option[Int] = None

  /** threshold for maximum alternate allele quality score sum in normal */
  @Argument(fullName="max_alt_alleles_in_normal_qscore_sum", shortName="", doc="threshold for maximum alternate allele quality score sum in normal", required=false, exclusiveOf="", validation="")
  var max_alt_alleles_in_normal_qscore_sum: Option[Int] = None

  /** threshold for maximum alternate allele fraction in normal */
  @Argument(fullName="max_alt_allele_in_normal_fraction", shortName="", doc="threshold for maximum alternate allele fraction in normal", required=false, exclusiveOf="", validation="")
  var max_alt_allele_in_normal_fraction: Option[Double] = None

  /** Format string for max_alt_allele_in_normal_fraction */
  @Argument(fullName="max_alt_allele_in_normal_fractionFormat", shortName="", doc="Format string for max_alt_allele_in_normal_fraction", required=false, exclusiveOf="", validation="")
  var max_alt_allele_in_normal_fractionFormat: String = "%s"

  /** Phred scale quality score constant to use in power calculations */
  @Argument(fullName="power_constant_qscore", shortName="", doc="Phred scale quality score constant to use in power calculations", required=false, exclusiveOf="", validation="")
  var power_constant_qscore: Option[Int] = None

  /** Absolute Copy Number Data, as defined by Absolute, to use in power calculations */
  @Argument(fullName="absolute_copy_number_data", shortName="", doc="Absolute Copy Number Data, as defined by Absolute, to use in power calculations", required=false, exclusiveOf="", validation="")
  var absolute_copy_number_data: File = _

  /** Allelic fraction constant to use in power calculations */
  @Argument(fullName="power_constant_af", shortName="", doc="Allelic fraction constant to use in power calculations", required=false, exclusiveOf="", validation="")
  var power_constant_af: Option[Double] = None

  /** Format string for power_constant_af */
  @Argument(fullName="power_constant_afFormat", shortName="", doc="Format string for power_constant_af", required=false, exclusiveOf="", validation="")
  var power_constant_afFormat: String = "%s"

  /** Call-stats output */
  @Output(fullName="out", shortName="o", doc="Call-stats output", required=false, exclusiveOf="", validation="")
  @Gather(classOf[org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction])
  var out: File = _

  /**
   * Short name of out
   * @return Short name of out
   */
  def o = this.out

  /**
   * Short name of out
   * @param value Short name of out
   */
  def o_=(value: File) { this.out = value }

  /** VCF output of mutation candidates */
  @Output(fullName="vcf", shortName="vcf", doc="VCF output of mutation candidates", required=false, exclusiveOf="", validation="")
  @Gather(classOf[VcfGatherFunction])
  var vcf: File = _

  /** Automatically generated index for vcf */
  @Output(fullName="vcfIndex", shortName="", doc="Automatically generated index for vcf", required=false, exclusiveOf="", validation="")
  @Gather(enabled=false)
  private var vcfIndex: File = _

  /** VCF file of DBSNP information */
  @Input(fullName="dbsnp", shortName="dbsnp", doc="VCF file of DBSNP information", required=false, exclusiveOf="", validation="")
  var dbsnp: Seq[File] = Nil

  /** Dependencies on any indexes of dbsnp */
  @Input(fullName="dbsnpIndexes", shortName="", doc="Dependencies on any indexes of dbsnp", required=false, exclusiveOf="", validation="")
  private var dbsnpIndexes: Seq[File] = Nil

  /** VCF file of COSMIC sites */
  @Input(fullName="cosmic", shortName="cosmic", doc="VCF file of COSMIC sites", required=false, exclusiveOf="", validation="")
  var cosmic: Seq[File] = Nil

  /** Dependencies on any indexes of cosmic */
  @Input(fullName="cosmicIndexes", shortName="", doc="Dependencies on any indexes of cosmic", required=false, exclusiveOf="", validation="")
  private var cosmicIndexes: Seq[File] = Nil

  /** VCF file of sites observed in normal */
  @Input(fullName="normal_panel", shortName="normal_panel", doc="VCF file of sites observed in normal", required=false, exclusiveOf="", validation="")
  var normal_panel: Seq[File] = Nil

  /** Dependencies on any indexes of normal_panel */
  @Input(fullName="normal_panelIndexes", shortName="", doc="Dependencies on any indexes of normal_panel", required=false, exclusiveOf="", validation="")
  private var normal_panelIndexes: Seq[File] = Nil

  /** write out coverage in WIGGLE format to this file */
  @Output(fullName="coverage_file", shortName="cov", doc="write out coverage in WIGGLE format to this file", required=false, exclusiveOf="", validation="")
  @Gather(classOf[org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction])
  var coverage_file: File = _

  /**
   * Short name of coverage_file
   * @return Short name of coverage_file
   */
  def cov = this.coverage_file

  /**
   * Short name of coverage_file
   * @param value Short name of coverage_file
   */
  def cov_=(value: File) { this.coverage_file = value }

  /** write out 20x of Q20 coverage in WIGGLE format to this file */
  @Output(fullName="coverage_20_q20_file", shortName="cov_q20", doc="write out 20x of Q20 coverage in WIGGLE format to this file", required=false, exclusiveOf="", validation="")
  @Gather(classOf[org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction])
  var coverage_20_q20_file: File = _

  /**
   * Short name of coverage_20_q20_file
   * @return Short name of coverage_20_q20_file
   */
  def cov_q20 = this.coverage_20_q20_file

  /**
   * Short name of coverage_20_q20_file
   * @param value Short name of coverage_20_q20_file
   */
  def cov_q20_=(value: File) { this.coverage_20_q20_file = value }

  /** write out power in WIGGLE format to this file */
  @Output(fullName="power_file", shortName="pow", doc="write out power in WIGGLE format to this file", required=false, exclusiveOf="", validation="")
  @Gather(classOf[org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction])
  var power_file: File = _

  /**
   * Short name of power_file
   * @return Short name of power_file
   */
  def pow = this.power_file

  /**
   * Short name of power_file
   * @param value Short name of power_file
   */
  def pow_=(value: File) { this.power_file = value }

  /** write out tumor read depth in WIGGLE format to this file */
  @Output(fullName="tumor_depth_file", shortName="tdf", doc="write out tumor read depth in WIGGLE format to this file", required=false, exclusiveOf="", validation="")
  @Gather(classOf[org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction])
  var tumor_depth_file: File = _

  /**
   * Short name of tumor_depth_file
   * @return Short name of tumor_depth_file
   */
  def tdf = this.tumor_depth_file

  /**
   * Short name of tumor_depth_file
   * @param value Short name of tumor_depth_file
   */
  def tdf_=(value: File) { this.tumor_depth_file = value }

  /** write out normal read depth in WIGGLE format to this file */
  @Output(fullName="normal_depth_file", shortName="ndf", doc="write out normal read depth in WIGGLE format to this file", required=false, exclusiveOf="", validation="")
  @Gather(classOf[org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction])
  var normal_depth_file: File = _

  /**
   * Short name of normal_depth_file
   * @return Short name of normal_depth_file
   */
  def ndf = this.normal_depth_file

  /**
   * Short name of normal_depth_file
   * @param value Short name of normal_depth_file
   */
  def ndf_=(value: File) { this.normal_depth_file = value }

  /** if a read has mismatching number of bases and base qualities, filter out the read instead of blowing up. */
  @Argument(fullName="filter_mismatching_base_and_quals", shortName="filterMBQ", doc="if a read has mismatching number of bases and base qualities, filter out the read instead of blowing up.", required=false, exclusiveOf="", validation="")
  var filter_mismatching_base_and_quals: Boolean = _

  /**
   * Short name of filter_mismatching_base_and_quals
   * @return Short name of filter_mismatching_base_and_quals
   */
  def filterMBQ = this.filter_mismatching_base_and_quals

  /**
   * Short name of filter_mismatching_base_and_quals
   * @param value Short name of filter_mismatching_base_and_quals
   */
  def filterMBQ_=(value: Boolean) { this.filter_mismatching_base_and_quals = value }

  override def freezeFieldValues() {
    super.freezeFieldValues()
    if (vcf != null && !org.broadinstitute.gatk.utils.io.IOUtils.isSpecialFile(vcf))
      if (!org.broadinstitute.gatk.utils.commandline.ArgumentTypeDescriptor.isCompressed(vcf.getPath))
        vcfIndex = new File(vcf.getPath + ".idx")
    dbsnpIndexes ++= dbsnp.filter(orig => orig != null).map(orig => new File(orig.getPath + ".idx"))
    cosmicIndexes ++= cosmic.filter(orig => orig != null).map(orig => new File(orig.getPath + ".idx"))
    normal_panelIndexes ++= normal_panel.filter(orig => orig != null).map(orig => new File(orig.getPath + ".idx"))
  }

  override def commandLine = super.commandLine + conditional(noop, "--noop", escape=true, format="%s") + conditional(enable_extended_output, "--enable_extended_output", escape=true, format="%s") + conditional(artifact_detection_mode, "--artifact_detection_mode", escape=true, format="%s") + optional("--tumor_sample_name", tumor_sample_name, spaceSeparated=true, escape=true, format="%s") + optional("--bam_tumor_sample_name", bam_tumor_sample_name, spaceSeparated=true, escape=true, format="%s") + optional("--normal_sample_name", normal_sample_name, spaceSeparated=true, escape=true, format="%s") + conditional(force_output, "--force_output", escape=true, format="%s") + conditional(force_alleles, "--force_alleles", escape=true, format="%s") + conditional(only_passing_calls, "--only_passing_calls", escape=true, format="%s") + optional("--initial_tumor_lod", initial_tumor_lod, spaceSeparated=true, escape=true, format=initial_tumor_lodFormat) + optional("--tumor_lod", tumor_lod, spaceSeparated=true, escape=true, format=tumor_lodFormat) + optional("--fraction_contamination", fraction_contamination, spaceSeparated=true, escape=true, format=fraction_contaminationFormat) + optional("--minimum_mutation_cell_fraction", minimum_mutation_cell_fraction, spaceSeparated=true, escape=true, format=minimum_mutation_cell_fractionFormat) + optional("--normal_lod", normal_lod, spaceSeparated=true, escape=true, format=normal_lodFormat) + optional("--normal_artifact_lod", normal_artifact_lod, spaceSeparated=true, escape=true, format=normal_artifact_lodFormat) + optional("--strand_artifact_lod", strand_artifact_lod, spaceSeparated=true, escape=true, format=strand_artifact_lodFormat) + optional("--strand_artifact_power_threshold", strand_artifact_power_threshold, spaceSeparated=true, escape=true, format=strand_artifact_power_thresholdFormat) + optional("--dbsnp_normal_lod", dbsnp_normal_lod, spaceSeparated=true, escape=true, format=dbsnp_normal_lodFormat) + optional("--somatic_classification_normal_power_threshold", somatic_classification_normal_power_threshold, spaceSeparated=true, escape=true, format=somatic_classification_normal_power_thresholdFormat) + optional("--minimum_normal_allele_fraction", minimum_normal_allele_fraction, spaceSeparated=true, escape=true, format=minimum_normal_allele_fractionFormat) + optional("--tumor_f_pretest", tumor_f_pretest, spaceSeparated=true, escape=true, format=tumor_f_pretestFormat) + optional("--min_qscore", min_qscore, spaceSeparated=true, escape=true, format="%s") + optional("--gap_events_threshold", gap_events_threshold, spaceSeparated=true, escape=true, format="%s") + optional("--heavily_clipped_read_fraction", heavily_clipped_read_fraction, spaceSeparated=true, escape=true, format=heavily_clipped_read_fractionFormat) + optional("--clipping_bias_pvalue_threshold", clipping_bias_pvalue_threshold, spaceSeparated=true, escape=true, format=clipping_bias_pvalue_thresholdFormat) + optional("--fraction_mapq0_threshold", fraction_mapq0_threshold, spaceSeparated=true, escape=true, format=fraction_mapq0_thresholdFormat) + optional("--pir_median_threshold", pir_median_threshold, spaceSeparated=true, escape=true, format=pir_median_thresholdFormat) + optional("--pir_mad_threshold", pir_mad_threshold, spaceSeparated=true, escape=true, format=pir_mad_thresholdFormat) + optional("--required_maximum_alt_allele_mapping_quality_score", required_maximum_alt_allele_mapping_quality_score, spaceSeparated=true, escape=true, format="%s") + optional("--max_alt_alleles_in_normal_count", max_alt_alleles_in_normal_count, spaceSeparated=true, escape=true, format="%s") + optional("--max_alt_alleles_in_normal_qscore_sum", max_alt_alleles_in_normal_qscore_sum, spaceSeparated=true, escape=true, format="%s") + optional("--max_alt_allele_in_normal_fraction", max_alt_allele_in_normal_fraction, spaceSeparated=true, escape=true, format=max_alt_allele_in_normal_fractionFormat) + optional("--power_constant_qscore", power_constant_qscore, spaceSeparated=true, escape=true, format="%s") + optional("--absolute_copy_number_data", absolute_copy_number_data, spaceSeparated=true, escape=true, format="%s") + optional("--power_constant_af", power_constant_af, spaceSeparated=true, escape=true, format=power_constant_afFormat) + optional("-o", out, spaceSeparated=true, escape=true, format="%s") + optional("-vcf", vcf, spaceSeparated=true, escape=true, format="%s") + conditional(no_cmdline_in_header, "-no_cmdline_in_header", escape=true, format="%s") + conditional(sites_only, "-sites_only", escape=true, format="%s") + conditional(bcf, "-bcf", escape=true, format="%s") + repeat("-dbsnp", dbsnp, formatPrefix=TaggedFile.formatCommandLineParameter, spaceSeparated=true, escape=true, format="%s") + repeat("-cosmic", cosmic, formatPrefix=TaggedFile.formatCommandLineParameter, spaceSeparated=true, escape=true, format="%s") + repeat("-normal_panel", normal_panel, formatPrefix=TaggedFile.formatCommandLineParameter, spaceSeparated=true, escape=true, format="%s") + optional("-cov", coverage_file, spaceSeparated=true, escape=true, format="%s") + optional("-cov_q20", coverage_20_q20_file, spaceSeparated=true, escape=true, format="%s") + optional("-pow", power_file, spaceSeparated=true, escape=true, format="%s") + optional("-tdf", tumor_depth_file, spaceSeparated=true, escape=true, format="%s") + optional("-ndf", normal_depth_file, spaceSeparated=true, escape=true, format="%s") + conditional(filter_mismatching_base_and_quals, "-filterMBQ", escape=true, format="%s")
}
