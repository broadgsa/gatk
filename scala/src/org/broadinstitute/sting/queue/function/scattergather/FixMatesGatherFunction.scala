package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.commandline.Input

class FixMatesGatherFunction extends GatherFunction {
  type GatherType = File

  @Input(doc="Picard FixMateInformation.jar.  At the Broad this can be found at /seq/software/picard/current/bin/FixMateInformation.jar.  Outside the broad see http://picard.sourceforge.net/")
  var picardFixMatesJar: String = _

  @Input(doc="Compression level 1-9", required=false)
  var picardMergeCompressionLevel: Option[Int] = None

  def commandLine = "java -Djava.io.tmpdir=/broad/shptmp/queue -jar %s%s%s%s".format(picardFixMatesJar,
    optional(" COMPRESSION_LEVEL=", picardMergeCompressionLevel), " VALIDATION_STRINGENCY=SILENT SO=coordinate OUTPUT=" + originalOutput, repeat(" INPUT=", gatherParts))
}
