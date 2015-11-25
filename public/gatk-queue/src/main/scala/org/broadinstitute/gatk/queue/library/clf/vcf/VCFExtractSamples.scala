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

package org.broadinstitute.gatk.queue.library.clf.vcf

import java.io.File
import collection.JavaConversions._
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.utils.text.XReadLines

class VCFExtractSamples(inVCF: File, outVCF: File, samples: List[String]) extends CommandLineFunction {
  @Input(doc="input VCF from which to extract samples") var inputVCF : File = inVCF
  @Output(doc="output VCF to write extracted samples to") var outputVCF : File = outVCF
  @Argument(doc="List of samples to extract from the VCF") var sampleList : List[String] = samples

  var sampleGrep : String = _

  def this(in: File, out: File, samples: File) = this(in,out, (new XReadLines(samples)).readLines.toList)

  override def freezeFieldValues = {
    this.logger.warn("Note: Using VCFExtractSamples invalidates AC/AF/AN annotations. This is an explicit warning.")
    sampleGrep = "'" + sampleList.reduceLeft(_ + "|" + _) + "'"
    super.freezeFieldValues
  }

  def commandLine = {

    var first : String = "head -n 500 %s | grep \\\\#\\\\# > %s".format(inputVCF.getAbsolutePath,outputVCF.getAbsolutePath)
    var second : String = "head -n 500 %s | grep \\\\#CHR | tr '\\t' '\\n' | awk '{print ++count\"\\t\"$1}' ".format(inputVCF.getAbsolutePath)
    second += "| egrep %s | awk '{print $1}' | tr '\\n' ',' | xargs -i cut -f1-9,\\{\\} %s | grep -v \\\\#\\\\# >> %s".format(sampleGrep,inputVCF.getAbsolutePath,outputVCF.getAbsolutePath)

    first+" ; "+second
  }

}