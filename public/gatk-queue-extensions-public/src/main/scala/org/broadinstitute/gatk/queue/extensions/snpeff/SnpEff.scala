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

package org.broadinstitute.gatk.queue.extensions.snpeff

import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import java.io.File
import org.broadinstitute.gatk.utils.commandline.{Argument, Output, Input}

/**
 * Basic snpEff support.
 * See: http://www.broadinstitute.org/gatk/guide/article?id=50
 */
class SnpEff extends JavaCommandLineFunction {
  javaMainClass = "ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff"

  @Input(doc="snp vcf4 file")
  var inVcf: File = _

  @Input(doc="config file with path to data dir", required=false)
  var config: File = _

  @Argument(doc="genome version")
  var genomeVersion: String = _

  @Argument(doc="verbose", required=false)
  var verbose = true

  @Argument(doc="onlyCoding", required=false)
  var onlyCoding = true

  @Output(doc="snp eff output")
  var outVcf: File = _

  override def commandLine = super.commandLine +
                             required("eff") +
                             conditional(verbose, "-v") +
                             required("-onlyCoding", onlyCoding.toString) +
                             optional("-c", config) +
                             required("-i", "vcf") +
                             required("-o", "vcf") +
                             required(genomeVersion) +
                             required(inVcf) +
                             required(">", escape=false) +
                             required(outVcf)
}
