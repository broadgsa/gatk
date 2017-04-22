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
 * "Reads a VCF/VCF.gz/BCF and removes all genotype information from it while retaining all site level information,
 * including annotations based on genotypes (e.g. AN, AF). Output can be any support variant format including .vcf, .vcf.gz or .bcf."
 */
class MakeSitesOnlyVcf extends JavaCommandLineFunction {
  analysisName = "MakeSitesOnlyVcf"
  javaMainClass = "picard.vcf.MakeSitesOnlyVcf"

  @Input(doc = "The input VCF files to analyze.", shortName = "input", fullName = "input_vcf_file", required = true)
  var input: File = _

  @Output(doc = "The output VCF which will not have any genotypes, but will keep the INFO field intact.", required = false)
  var output: File = _

  @Argument(shortName = "S", doc = "Optionally one or more samples to retain when building the \'sites-only\' VCF.", required = false)
  var samples: List[String] = _

  var validationStringency = ValidationStringency.SILENT
  var compressionLevel: Option[Int] = None
  var createIndex: Option[Boolean] = None
  var maxRecordsInRam: Option[Int] = None
  var assumeSorted: Option[Boolean] = None

  override def commandLine = super.commandLine +
    required("INPUT=", input, spaceSeparated = false) +
    required("TMP_DIR=" + jobTempDir) +
    optional("OUTPUT=", output, spaceSeparated = false) +
    repeat("SAMPLE=", samples, spaceSeparated = false) +
    optional("COMPRESSION_LEVEL=", compressionLevel, spaceSeparated = false) +
    optional("VALIDATION_STRINGENCY=", validationStringency, spaceSeparated = false) +
    optional("MAX_RECORDS_IN_RAM=", maxRecordsInRam, spaceSeparated = false) +
    optional("ASSUME_SORTED=", assumeSorted, spaceSeparated = false) +
    optional("CREATE_INDEX=", createIndex, spaceSeparated = false)

}
