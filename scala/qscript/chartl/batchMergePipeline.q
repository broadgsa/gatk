import java.io.{FileReader, BufferedReader}
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk.CommandLineGATK
import org.broadinstitute.sting.queue.pipeline.{ProjectManagement, BamProcessing, VariantCalling}
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.utils.yaml.YamlUtils

class batchMergePipeline extends QScript {
  batchMerge =>

  @Argument(doc="VCF list",shortName="vcfs") var vcfList: File = _
  @Argument(doc="bam list",shortName="bams") var bamList: File = _
  @Argument(doc="sting dir",shortName="sting") var stingDir: String = _
  @Argument(doc="reference file",shortName="ref") var ref: File = _
  @Argument(doc="batched output",shortName="batch") var batchOut: File = _

  def addAll( cmds : List[CommandLineFunction]) : Unit = { cmds.foreach( c => add(c)) }

  def script = {

    var vcfs : List[File] = extractFileEntries(vcfList)
    var bams : List[File] = extractFileEntries(bamList)
    var pmLib = new ProjectManagement(stingDir)
    addAll(pmLib.MergeBatches(vcfs,bams,batchOut,ref,20))
  }

  def extractFileEntries(in: File): List[File] = {
    return (new XReadLines(in)).readLines.toList.map( new File(_) )
  }
}