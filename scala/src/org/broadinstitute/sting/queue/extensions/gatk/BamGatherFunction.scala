package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.scattergather.GatherFunction
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction

/**
 * Merges BAM files using Picards MergeSamFiles.jar.
 * At the Broad the jar can be found at /seq/software/picard/current/bin/MergeSamFiles.jar.  Outside the broad see http://picard.sourceforge.net/")
 */
class BamGatherFunction extends GatherFunction with PicardBamJarFunction {
  this.assumeSorted = Some(true)
  protected def inputBams = gatherParts
  protected def outputBam = originalOutput
}
