package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * A pipeline for Queue that runs a custom walker on the classpath.
 * NOTE: This code is an unsupported example for soliciting feedback on how to improve Queue.
 * Future syntax will simplify running the GATK so please expect the syntax below to change significantly.
 */
class ExampleCustomWalker extends QScript {
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
    val customWalker = new CommandLineGATK {
      // Set the name of your walker, for example this will be passed as -T MyCustomWalker
      this.analysis_type = "MyCustomWalker"

      // If your walker is already on the classpath you shouldn't need to do anything else

      // If your walker is in a GATK jar that is for some reason NOT on the classpath
      // nor referenced in the Queue.jar's, specify the jar file here
      //this.jarFile = "myGATK.jar"

      // If your walker needs a custom classpath, specify it here
      //this.javaClasspath = List("myClasses")
    }

    customWalker.reference_sequence = referenceFile
    customWalker.input_file = bamFiles

    // Add the newly created function to the pipeline.
    add(customWalker)
  }
}
