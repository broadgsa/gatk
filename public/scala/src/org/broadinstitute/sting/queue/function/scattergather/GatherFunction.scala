package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.{Input, Output}
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.utils.io.IOUtils
import collection.JavaConversions._

/**
 * Base class for Gather command line functions.
 */
trait GatherFunction extends QFunction {
  var originalFunction: ScatterGatherableFunction = _

  @Output(doc="The original output of the scattered function")
  var originalOutput: File = _

  @Input(doc="Parts to gather back into the original output")
  var gatherParts: List[File] = Nil

  @Input(doc="Other log files that will be gathered before this output", required=false)
  var originalLogFiles: List[File] = Nil

  /**
   * Called to initialize the gather function values after all other values have been setup for this function.
   */
  def init() {}

  /**
   * Waits for gather parts to propagate over NFS or throws an exception.
   */
  protected def waitForGatherParts() {
    val missing = IOUtils.waitFor(gatherParts, 120)
    if (!missing.isEmpty)
      throw new QException("Unable to find gather inputs: " + missing.mkString(", "))
  }
}
