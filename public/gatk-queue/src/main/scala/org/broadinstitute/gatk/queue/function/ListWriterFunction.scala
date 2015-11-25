/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.function

import org.broadinstitute.gatk.utils.commandline.{Input, Output}
import java.io.{PrintWriter, File}
import org.apache.commons.io.IOUtils

/**
 * Writes a list of inputs to an output file.
 * Custom formats can override addFile.
 */
class ListWriterFunction extends InProcessFunction {
  analysisName = "WriteList"

  @Input(doc="input files") var inputFiles: Seq[File] = Nil
  @Output(doc="output file") var listFile: File = _

  def run() {
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
