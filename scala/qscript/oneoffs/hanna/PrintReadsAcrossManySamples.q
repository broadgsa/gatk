import java.io.PrintWriter
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException

/**
 * A pipeline for Queue that runs a custom walker outside of the GATK jar.
 * NOTE: This code is an unsupported example for soliciting feedback on how to improve Queue.
 * Future syntax will simplify running the GATK so please expect the syntax below to change significantly.
 */
class PrintReadsAcrossManySamples extends QScript {
  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = null

  // NOTE: Do not initialize List, Set, or Option to null
  // as you won't be able to update the collection.
  // By default set:
  //   List[T] = Nil
  //   Set[T] = Set.empty[T]
  //   Option[T] = None
  @Input(doc="One or more bam files.", shortName="I")
  var bamFile: File = _

  @Input(doc="Name of the test case", fullName="test_case",required=false)
  var testCaseName: String = "."

  @Input(doc="Max number of bam files to process", fullName="max_bams",required=false)
  var maxBams = 1

  @Input(doc="Step size",fullName="step_size",required=false)
  var stepSize = 1


  def script = {
    var lines = scala.io.Source.fromFile(bamFile).getLines.map(new File(_)).toList

    for(numBams <- 1 to (maxBams min lines.size) by stepSize) {
      val basePath = "%s/%03d_bams".format(testCaseName,numBams)

      val writeBamList = new ListWriterFunction
      writeBamList.inputFiles = lines.take(numBams)
      writeBamList.listFile = new File(basePath+"/bams.list")
      add(writeBamList)

      // Create the function that we can run.
      val printreads = new PrintReads

      printreads.jobOutputFile = new File(basePath+"/PrintReads.out")
      printreads.input_file = List(writeBamList.listFile)
      printreads.reference_sequence = referenceFile
      printreads.out = new File("/dev/null")
      printreads.memoryLimit = 8

      add(printreads)
    }

  }
}
