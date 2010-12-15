import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.QScript

import org.broadinstitute.sting.queue.library.ipf.vcf.VCFInfoToTable
import collection.JavaConversions._

class Vcf2Table extends QScript {
  @Argument(shortName="vcf",doc="VCF file",required=true) var vcf : File = _
  @Argument(shortName="f",doc="Info fields to extract",required=false) var fields : java.util.List[String] = new java.util.ArrayList[String]
  @Argument(shortName="o",doc="Output file",required=true) var output : File = _
  @Argument(shortName="useFilters",doc="Use filtered sites?",required=false) var useFilters : Boolean = false
  @Hidden @Argument(shortName="pass",doc="set the hack filter string to this value",required=false) var filterString : String = "PASS"
  @Hidden @Argument(shortName="notfound",doc="set the hack no entry string to this value",required=false) var keyNotFound : String = "NA"

  def script = {
    var vcf2table : VCFInfoToTable = new VCFInfoToTable(vcf,output,fields,useFilters)
    vcf2table.PF_KEY = filterString
    vcf2table.NO_KEY = keyNotFound
    add(vcf2table)
    this.functions.foreach(u => logger.debug("added: %s%n".format(u.toString)))
  }
}