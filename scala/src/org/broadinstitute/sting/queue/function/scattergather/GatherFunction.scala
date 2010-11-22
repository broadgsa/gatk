package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.{Input, Output}
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.queue.util.IOUtils

/**
 * Base class for Gather command line functions.
 */
trait GatherFunction extends QFunction {
  @Input(doc="Parts to gather back into the original output")
  var gatherParts: List[File] = Nil

  @Output(doc="The original output of the scattered function")
  var originalOutput: File = _

  /**
   * Waits for gather parts to propagate over NFS or throws an exception.
   */
  protected def waitForGatherParts = {
    val missing = IOUtils.waitFor(gatherParts, 120)
    if (!missing.isEmpty)
      throw new QException("Unable to find gather inputs: " + missing.mkString(", "))
  }
}
