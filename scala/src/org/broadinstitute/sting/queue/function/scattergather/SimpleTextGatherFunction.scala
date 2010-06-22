package org.broadinstitute.sting.queue.function.scattergather

import java.io.File

class SimpleTextGatherFunction extends GatherFunction {
  type GatherType = File

  // TODO: Write a text merging utility that takes into account headers.
  def commandLine = "mergeText.sh %s%s".format(originalOutput, repeat(" ", gatherParts))
}
