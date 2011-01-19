import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * A pipeline for Queue that runs a custom walker outside of the GATK jar.
 * NOTE: This code is an unsupported example for soliciting feedback on how to improve Queue.
 * Future syntax will simplify running the GATK so please expect the syntax below to change significantly.
 */
class DoC extends QScript {
  // The full packaged jar should be used.
  // You can build this jar via 'ant package' and then find it under
  // 'Sting/dist/packages/GenomeAnalysisTK-*/GenomeAnalysisTK.jar'
  @Input(doc="The path to the packaged GenomeAnalysisTK.jar file.", shortName="gatk")
  var gatkJar: File = null

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = null

  // NOTE: Do not initialize List, Set, or Option to null
  // as you won't be able to update the collection.
  // By default set:
  //   List[T] = Nil
  //   Set[T] = Set.empty[T]
  //   Option[T] = None
  @Input(doc="One or more bam files.", shortName="I")
  var bamFiles: List[File] = Nil

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervalsString: List[String] = List("2:87000001-90000000")

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait DepthOfCoverageArguments extends CommandLineGATK {
    this.jarFile = DoC.this.gatkJar
    this.reference_sequence = DoC.this.referenceFile
    this.intervalsString = DoC.this.intervalsString
    this.memoryLimit = Some(8)
  }


  def script = {
    // Create the four function that we can run.
    val doc = new DepthOfCoverage with DepthOfCoverageArguments
    doc.downsampling_type = Some(DownsampleType.NONE)
    doc.omitLocusTable = true
    doc.omitIntervals = true
    doc.omitSampleSummary = true        

    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.
    //doc.scatterCount = 3
    doc.input_file = DoC.this.bamFiles
    doc.out = new File("doc-all.out")

    add(doc)
  }

}
