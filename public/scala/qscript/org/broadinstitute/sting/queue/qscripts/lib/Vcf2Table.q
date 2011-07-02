package org.broadinstitute.sting.queue.qscripts.lib

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