package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.{Input, Output}
import org.broadinstitute.sting.queue.function.QFunction

/**
 * Base class for Gather command line functions.
 */
trait GatherFunction extends QFunction {
  @Input(doc="Parts to gather back into the original output")
  var gatherParts: List[File] = Nil

  @Output(doc="The original output of the scattered function")
  var originalOutput: File = _
}
