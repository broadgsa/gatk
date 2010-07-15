package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.Input

class BamGatherFunction extends GatherFunction {
  type GatherType = File

  @Input(doc="Picard MergeSamFiles.jar.  At the Broad this can be found at /seq/software/picard/current/bin/MergeSamFiles.jar.  Outside the broad see http://picard.sourceforge.net/")
  var picardMergeSamFilesJar: String = _

  @Input(doc="Compression level 1-9", required=false)
  var picardMergeCompressionLevel: Option[Int] = None

  def commandLine = "java -jar %s%s%s%s".format(picardMergeSamFilesJar,
    optional(" COMPRESSION_LEVEL=", picardMergeCompressionLevel), " AS=true VALIDATION_STRINGENCY=SILENT SO=coordinate OUTPUT=" + originalOutput, repeat(" INPUT=", gatherParts))
}
