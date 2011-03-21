package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.scattergather.GatherFunction
import org.broadinstitute.sting.queue.extensions.picard.PicardBamEmbeddedFunction

/**
 * Merges BAM files using Picards net.sf.picard.sam.MergeSamFiles.
 */
class BamGatherFunction extends GatherFunction with PicardBamEmbeddedFunction {
  this.mainClass = "net.sf.picard.sam.MergeSamFiles"
  this.assumeSorted = Some(true)
  protected def inputBams = gatherParts
  protected def outputBam = originalOutput

  override def init() {
    // Whatever the original function can handle, merging *should* do less.
    this.memoryLimit = originalFunction.memoryLimit
  }
}
