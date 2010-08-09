package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.JarCommandLineFunction
import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.queue.function.scattergather.GatherFunction

/**
 * Merges BAM files using Picards MergeSampFiles.jar.
 * At the Broad the jar can be found at /seq/software/picard/current/bin/MergeSamFiles.jar.  Outside the broad see http://picard.sourceforge.net/")
 */
class BamGatherFunction extends GatherFunction with JarCommandLineFunction {
  @Argument(doc="Compression level 1-9", required=false)
  var compressionLevel: Option[Int] = None

  override def commandLine = super.commandLine + "%s%s%s".format(
    optional(" COMPRESSION_LEVEL=", compressionLevel), " AS=true VALIDATION_STRINGENCY=SILENT SO=coordinate OUTPUT=" + originalOutput, repeat(" INPUT=", gatherParts))
}
