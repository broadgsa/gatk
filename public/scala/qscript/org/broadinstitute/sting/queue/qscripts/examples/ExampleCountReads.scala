package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * An introductory pipeline for Queue.
 * Runs the GATK CountReads individually and across a set of bams.
 * All bams must have the same reference.
 */
class ExampleCountReads extends QScript {
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

    // Run CountReads for all bams jointly.

    // Create a new CountReads from the Queue GATK Extensions.
    // The names of walkers are the same as you would use for '-T <WalkerName>'
    val jointCountReads = new CountReads

    // Each field in the extensions is based off of the full form of the arguments.
    // To get the list of arguments and their descriptions run
    // java -jar <path to GenomeAnalysisTK.jar> -T <WalkerName> -help
    jointCountReads.reference_sequence = referenceFile

    // GATK inputs that take more than one file will have a singular name which
    // matches the full form of the argument, but will actually be a scala List[]
    jointCountReads.input_file = bamFiles

    // Add the newly created function to the pipeline.
    add(jointCountReads)

    // If there is more than one BAM, also run CountReads once for each bam.
    if (bamFiles.size > 1) {
      for (bamFile <- bamFiles) {
        val singleCountReads = new CountReads
        singleCountReads.reference_sequence = referenceFile
        // ':+' is the scala List append operator
        singleCountReads.input_file :+= bamFile
        add(singleCountReads)
      }
    }
  }
}
