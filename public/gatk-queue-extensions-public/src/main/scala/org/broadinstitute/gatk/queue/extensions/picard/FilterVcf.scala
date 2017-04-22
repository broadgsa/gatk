/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.queue.extensions.picard

import java.io.File

import htsjdk.samtools.ValidationStringency
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import org.broadinstitute.gatk.utils.commandline.{Argument, Input, Output}

/**
  * Reads a VCF/VCF.gz/BCF and and filters it quickly based on some simple hard filters
  */
class FilterVcf extends JavaCommandLineFunction {
  analysisName = "FilterVcf"
  javaMainClass = "picard.vcf.FilterVcf"

  @Input(doc = "The input VCF files to filter.", shortName = "input", fullName = "input_vcf_file", required = true)
  var input: File = _

  @Output(doc = "The output VCF which will have it's FILTER field updated", required = false)
  var output: File = _


  @Argument(doc = "The minimum allele balance acceptable before filtering a site. Allele balance is calculated for heterozygotes as " +
    "the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygous " +
    "genotypes at the same locus are measured independently. The locus is filtered if any allele balance is below the limit.", required = false)
  var minAb: Option[Double] = _

  @Argument(doc = "The minimum sequencing depth supporting a genotype before the genotype will be filtered out.", required = false)
  var minDp: Option[Double] = _

  @Argument(doc = "The minimum genotype quality that must be achieved for a sample otherwise the genotype will be filtered out.", required = false)
  var minGQ: Option[Double] = _

  @Argument(doc = "The maximum phred scaled fisher strand value before a site will be filtered out.", required = false)
  var maxFs: Option[Double] = _

  @Argument(doc = "The minimum QD value to accept or otherwise filter out the variant.", required = false)
  var minQd: Option[Double] = _

  var validationStringency = ValidationStringency.SILENT
  var compressionLevel: Option[Int] = None
  var createIndex: Option[Boolean] = None
  var maxRecordsInRam: Option[Int] = None
  var assumeSorted: Option[Boolean] = None

  override def commandLine = super.commandLine +
    required("INPUT=", input, spaceSeparated = false) +
    required("TMP_DIR=" + jobTempDir) +
    optional("OUTPUT=", output, spaceSeparated = false) +
    optional("MIN_AB=", minAb, spaceSeparated = false) +
    optional("MIN_DP=", minDp, spaceSeparated = false) +
    optional("MIN_GQ=", minGQ, spaceSeparated = false) +
    optional("MIN_FS=", maxFs, spaceSeparated = false) +
    optional("MIN_QD=", minQd, spaceSeparated = false) +
    optional("COMPRESSION_LEVEL=", compressionLevel, spaceSeparated = false) +
    optional("VALIDATION_STRINGENCY=", validationStringency, spaceSeparated = false) +
    optional("MAX_RECORDS_IN_RAM=", maxRecordsInRam, spaceSeparated = false) +
    optional("ASSUME_SORTED=", assumeSorted, spaceSeparated = false) +
    optional("CREATE_INDEX=", createIndex, spaceSeparated = false)

}
