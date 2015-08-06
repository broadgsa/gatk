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

package org.broadinstitute.gatk.queue.qscripts.lib

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.gatk.{RodBind, VariantsToTable}
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._

class Vcf2Table extends QScript {
  @Argument(shortName="vcf",doc="VCF file",required=true) var vcf : File = _
  @Argument(shortName="f",doc="Info fields to extract",required=false) var fields : java.util.List[String] = new java.util.ArrayList[String]
  @Argument(shortName="o",doc="Output file",required=true) var output : File = _
  @Argument(shortName="useFilters",doc="Use filtered sites?",required=false) var useFilters : Boolean = false
  @Argument(shortName="r",doc="Reference file") var ref : File = _
  @Argument(shortName="i",doc="Intervals",required=false) var ints : java.util.List[File] = new java.util.ArrayList[File]
  @Argument(shortName="g",doc="gatk jar",required=true) var gatk: File = _


  def script = {
    var vcf2table : VariantsToTable = new VariantsToTable
    vcf2table.rodBind :+= new RodBind("variant","vcf",vcf)
    vcf2table.reference_sequence = ref
    vcf2table.intervals = ints.toList
    vcf2table.raw = useFilters
    vcf2table.out = output
    vcf2table.F = fields.toList
    vcf2table.jarFile = gatk
    add(vcf2table)

  }
}