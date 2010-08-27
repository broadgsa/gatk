import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * A pipeline for Queue that runs a custom walker outside of the GATK jar.
 * NOTE: This code is an unsupported example for soliciting feedback on how to improve Queue.
 * Future syntax will simplify running the GATK so please expect the syntax below to change significantly.
 */
class ExampleCustomWalker extends QScript {
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

  /**
   * In script, you create and then add() functions to the pipeline.
   */
  def script = {
    val myClasses = "myClassDir"

    val customWalker = new CommandLineGATK {
      // Set the name of your walker, for example this will be passed as -T MyCustomWalker
      this.analysis_type = "MyCustomWalker"
      // NOTE: At this time, you still need to specify the GATK jar or the pipeline won't validate.
      this.jarFile = gatkJar
      override def javaExecutable = "org.broadinstitute.sting.gatk.CommandLineGATK"
      override def javaOpts = "%s -cp %s:%s".format(super.javaOpts, gatkJar, myClasses)
    }

    customWalker.reference_sequence = referenceFile
    customWalker.input_file = bamFiles

    // Add the newly created function to the pipeline.
    add(customWalker)
  }
}
