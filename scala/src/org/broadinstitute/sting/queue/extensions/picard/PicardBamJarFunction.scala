package org.broadinstitute.sting.queue.extensions.picard

import org.broadinstitute.sting.queue.function.JarCommandLineFunction
import java.io.File

/**
 * Wraps a Picard jar that operates on BAM files.
 * See http://picard.sourceforge.net/ for more info.
 *
 * Since the jar files take slightly different arguments
 * some values are optional.
 */
trait PicardBamJarFunction extends JarCommandLineFunction {
  var validationStringency = "SILENT"
  var sortOrder = "coordinate"
  var compressionLevel: Option[Int] = None
  var maxRecordsInRam: Option[Int] = None
  var assumeSorted: Option[Boolean] = None

  protected def inputBams: List[File]
  protected def outputBam: File

  override def commandLine = super.commandLine + "%s%s%s%s%s%s%s%s".format(
    optional(" COMPRESSION_LEVEL=", compressionLevel), optional(" VALIDATION_STRINGENCY=", validationStringency),
    optional(" SO=", sortOrder), optional( " MAX_RECORDS_IN_RAM=", maxRecordsInRam), optional(" ASSUME_SORTED=", assumeSorted),
    " OUTPUT=" + outputBam, repeat(" INPUT=", inputBams), " TMP_DIR=" + jobTempDir)
}
