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
import org.broadinstitute.gatk.utils.commandline.{Argument, Output, Input}
import org.broadinstitute.gatk.queue.function.CommandLineFunction

class VCFExtractIntervals(inVCF: File, outList: File, passOnly: Boolean ) extends CommandLineFunction {
  def this(vcf: File, list: File ) = this(vcf,list,false)
  def this(vcf: File) = this(vcf, new File(vcf.getAbsolutePath.replace("vcf","intervals.list")),false)
  def this(vcf: File, removeFilters: Boolean) = this(vcf, new File(vcf.getAbsolutePath.replace("vcf","intervals.list")), removeFilters)

  @Input(doc="The VCF from which to extract an interval list") var inputVCF : File = inVCF
  @Output(doc="The file to write the interval list to") var outputList : File = outList
  @Argument(doc="Whether to use all sites, or only the unfiltered sites") var usePFOnly: Boolean = passOnly

  def commandLine = {
    if ( usePFOnly ) "grep PASS %s | awk '{print $1\":\"$2}' | uniq > %s".format(inputVCF.getAbsolutePath,outputList.getAbsolutePath)
    else "grep -v \\\\# %s | awk '{print $1\":\"$2}' | uniq > %s".format(inputVCF.getAbsolutePath,outputList.getAbsolutePath)
  }

}