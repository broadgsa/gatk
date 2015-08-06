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

package org.broadinstitute.gatk.queue.extensions.picard

/*
* Copyright (c) 2012 The Broad Institute
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

import org.broadinstitute.gatk.utils.commandline._

import java.io.File

class FastqToSam extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction /*with PicardBamFunction*/ {
  analysisName = "FastqToSam"
  javaMainClass = "picard.sam.FastqToSam"

  @Input(shortName = "fq1", fullName = "input_fq_file1", required = true, doc = "Input Fastq file to extract reads from (single-end fastq or, if paired, first end of the pair fastq)")
  var fastq: File = _

  @Input(shortName = "fq2", fullName = "input_fq_file2", required = false, doc = "Input Fastq file to extract reads from (if paired, second end of the pair fastq).")
  var secondEndFastQ: File = _

  @Output(shortName = "bam", fullName = "output_bam_file", required = true, doc = "Output bam file .")
  var bam: File = _

  @Argument(shortName = "SM", fullName = "SM", required = false, doc = "SM")
  var SM: String = "SM"

  @Argument(shortName = "LIB", fullName = "LIB", required = false, doc = "LIB")
  var LIB: String = "LIB"

  @Argument(shortName = "PU", fullName = "PU", required = false, doc = "PU")
  var PU: String = "PU"

  @Argument(shortName = "RG", fullName = "RG", required = false, doc = "RG")
  var RG: String = "RG"

  @Argument(shortName = "PL", fullName = "PL", required = false, doc = "PL")
  var PL: String = "illumina"

  @Argument(shortName = "CN", fullName = "CN", required = false, doc = "CN")
  var CN: String = "CN"


//  override def inputBams = Seq(fastq)
//  override def outputBam = bam
//  this.sortOrder = null
  val createIndex:Boolean = true
  override def commandLine = super.commandLine +
    required("FASTQ=" + fastq) +
    optional("FASTQ2=", secondEndFastQ, spaceSeparated=false) +
    required("OUTPUT=" + bam) +
    optional("READ_GROUP_NAME=", RG, spaceSeparated=false) +
    required("SAMPLE_NAME=" + SM) +
    optional("LIBRARY_NAME=", LIB, spaceSeparated=false) +
    optional("PLATFORM_UNIT=", PU, spaceSeparated=false) +
    optional("PLATFORM=", PL, spaceSeparated=false) +
    optional("CREATE_INDEX=", createIndex, spaceSeparated=false) +
    optional("SEQUENCING_CENTER=", CN, spaceSeparated=false)
}