package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.commandline.{Input, Output}
import java.io.{PrintWriter, File}
import org.apache.commons.io.IOUtils

/**
 * Writes a list of inputs to an output file.
 * Custom formats can override addFile.
 */
class ListWriterFunction extends InProcessFunction {
  @Input(doc="input files") var inputFiles: List[File] = Nil
  @Output(doc="output file") var listFile: File = _

  def run {
    val writer = new PrintWriter(listFile)
    try {
      for (inputFile <- inputFiles)
        addFile(writer, inputFile)
    } finally {
      IOUtils.closeQuietly(writer)
    }
  }

  /**
   * Adds the inputFile to the output list.
   * @param writer Output file.
   * @param inputFile File to add to the output file.
   */
  def addFile(writer: PrintWriter, inputFile: File) {
    writer.println(inputFile.toString)
  }
}
