import org.broadinstitute.sting.queue.library.ipf.vcf.{VCFExtractIntervals, VCFExtractSamples, VCFSimpleMerge, VCFExtractSites}
import org.broadinstitute.sting.queue.library.ipf.SortByRef
import org.broadinstitute.sting.queue.library.ipf.intervals.ExpandIntervals
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._

// todo -- should the argument collection on which this runs be generated at compile-time into extensions??
// todo -- maybe a compile-time generated enum of available library functions? (ipf of course)
class QTools extends QScript {
  @Argument(doc="Tool to run",shortName="T", required=true) var qtool : String = _
  @Argument(doc="input VCF",shortName="ivcf",required=false) var inVCF : File = _
  @Argument(doc="input VCF files",shortName="vcfs",required=false) var inVCFs : String = _
  @Argument(doc="output file",shortName="out",required=true) var output : File = _
  @Argument(doc="reference file",shortName="ref",required=false) var ref : File = _
  @Argument(doc="The samples to extract",shortName="sm",required=false) var samples : String = _
  @Argument(doc="Keep filtered sites when merging or extracting?",shortName="kf",required=false) var keepFilters : Boolean = false
  @Argument(doc="Input interval list (not used with VCF tools)",shortName="il",required=false) var intervalList : File = _
  @Argument(doc="interval list expand start",shortName="il_start",required=false) var ilStart : Int = 1
  @Argument(doc="interval list expand size",shortName="il_size",required=false) var ilSize : Int = 50
  // todo -- additional arguments or argument collection

  def script = {
    if ( qtool.equals("VCFExtractSites") ) {
      runVCFExtractSites
    }

    if ( qtool.equals("VCFSimpleMerge") ) {
      runVCFSimpleMerge
    }

    if ( qtool.equals("VCFExtractSamples") ) {
      runVCFExtractSamples
    }

    if ( qtool.equals("VCFExtractIntervals") ) {
      runVCFExtractIntervals
    }

    if ( qtool.equals("SortByRef") ) {
      runSortByRef
    }

    if ( qtool.equals("ExpandTargets") ) {
      runExpandTargets
    }
  }

  def runVCFExtractSites = {
    var ves : VCFExtractSites = new VCFExtractSites(inVCF,output)
    add(ves)
  }

  def runVCFSimpleMerge = {
    var vsm : VCFSimpleMerge = new VCFSimpleMerge
    vsm.vcfs = inVCFs.split(",").toList.map(new File(_))
    vsm.outVCF = output
    vsm.fai = new File(ref.getAbsolutePath+".fai")

    add(vsm)
  }

  def runVCFExtractSamples = {
    var ves : VCFExtractSamples = new VCFExtractSamples(inVCF,output,samples.split(",").toList)
    add(ves)
  }

  def runVCFExtractIntervals = {
    var vei : VCFExtractIntervals = new VCFExtractIntervals(inVCF,output,keepFilters)
    add(vei)
  }

  def runSortByRef = {
    var sbr : SortByRef = new SortByRef(inVCF,new File(ref.getAbsolutePath+".fai"),output)
    add(sbr)
  }

  def runExpandTargets = {
    var ets : ExpandIntervals = new ExpandIntervals(intervalList,ilStart,ilSize,output,ref,"INTERVALS","INTERVALS")
    add(ets)
  }
}