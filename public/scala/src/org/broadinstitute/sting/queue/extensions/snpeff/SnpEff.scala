/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.extensions.snpeff

import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import java.io.File
import org.broadinstitute.sting.commandline.{Argument, Output, Input}

/**
 * Basic snpEff support.
 * See: http://www.broadinstitute.org/gsa/wiki/index.php/Adding_Genomic_Annotations_Using_SnpEff_and_VariantAnnotator
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

  @Output(doc="snp eff output")
  var outVcf: File = _

  override def commandLine = Array(
    super.commandLine,
    " eff",
    if (verbose) " -v" else "",
    optional(" -c ", config),
    " -i vcf -o vcf %s %s > %s".format(genomeVersion, inVcf, outVcf)
  ).mkString
}
