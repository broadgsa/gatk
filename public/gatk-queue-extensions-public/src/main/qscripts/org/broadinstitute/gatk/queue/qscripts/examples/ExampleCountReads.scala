/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.qscripts.examples

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._

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
  def script() {

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

    // Set the memory limit. Also acts as a memory request on LSF and GridEngine.
    jointCountReads.memoryLimit = 1

    // Add the newly created function to the pipeline.
    add(jointCountReads)

    // If there is more than one BAM, also run CountReads once for each bam.
    if (bamFiles.size > 1) {
      for (bamFile <- bamFiles) {
        val singleCountReads = new CountReads
        singleCountReads.reference_sequence = referenceFile
        // ':+' is the scala List append operator
        singleCountReads.input_file :+= bamFile
        singleCountReads.memoryLimit = 1
        add(singleCountReads)
      }
    }
  }
}
