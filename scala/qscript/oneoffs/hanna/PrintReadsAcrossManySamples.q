import java.io.PrintWriter
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException

/**
 * A pipeline for Queue that runs a custom walker outside of the GATK jar.
 * NOTE: This code is an unsupported example for soliciting feedback on how to improve Queue.
 * Future syntax will simplify running the GATK so please expect the syntax below to change significantly.
 */
class PrintReadsAcrossManySamples extends QScript {
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

  @Input(doc="Name of the test case", fullName="test_case",required=false)
  var testCaseName: String = "."

  @Input(doc="Max number of bam files to process", fullName="max_bams",required=false)
  var maxBams = 1

  @Input(doc="Step size",fullName="step_size",required=false)
  var stepSize = 1

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait PrintReadsAcrossManySamplesArguments extends CommandLineGATK {
    this.jarFile = PrintReadsAcrossManySamples.this.gatkJar
    this.reference_sequence = PrintReadsAcrossManySamples.this.referenceFile
    this.memoryLimit = 8
  }


  def script = {
    if(bamFiles.size != 1)
      throw new ReviewedStingException("-I argument must consist of exactly one file containing a list of BAM files.");

    var lines: List[String] = List[String]()
    for(line <- scala.io.Source.fromFile(bamFiles(0)).getLines) {
      lines = lines ::: List(line)
    }

    for(numBams <- 1 to math.min(maxBams,lines.size) by stepSize) {
      val dir = new File(testCaseName + "/%03d_bams".format(numBams))
      dir.mkdir()

      val file = new File(dir,"bams.list")
      val writer = new PrintWriter(file)
      
      for(bamIndex <- 0 to numBams-1)
        writer.println(lines(bamIndex))

      writer.close()

      // Create the function that we can run.
      val printreads = new PrintReads with PrintReadsAcrossManySamplesArguments

      printreads.jobOutputFile = new File(dir, "PrintReads.out")      
      printreads.input_file = List(file)
      printreads.reference_sequence = referenceFile
      printreads.out = new File("/dev/null")

      add(printreads)
    }

  }
}
