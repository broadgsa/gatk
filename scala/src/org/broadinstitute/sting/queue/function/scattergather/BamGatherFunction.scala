package org.broadinstitute.sting.queue.function.scattergather

import java.io.File

class BamGatherFunction extends GatherFunction {
  type GatherType = File

  def commandLine = "samtools merge %s%s".format(originalOutput, repeat(" ", gatherParts))
}
