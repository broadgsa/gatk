package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.commandline.Argument

/**
 * Merges a text file.
 * The script can be changed by setting rmdirScript.
 * By default uses mergeText.sh in Sting/shell.
 * The format of the call is <mergeTextScript> <file_output> <file_1> [.. <file_n>]
 */
class SimpleTextGatherFunction extends GatherFunction {
  @Argument(doc="merge text script")
  var mergeTextScript = "mergeText.sh"

  def commandLine = "%s %s%s".format(mergeTextScript, originalOutput, repeat(" ", gatherParts))
}
