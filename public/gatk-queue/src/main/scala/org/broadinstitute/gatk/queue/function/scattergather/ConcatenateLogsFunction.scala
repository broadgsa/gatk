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

package org.broadinstitute.gatk.queue.function.scattergather

import org.broadinstitute.gatk.queue.function.InProcessFunction
import org.broadinstitute.gatk.queue.QException
import org.broadinstitute.gatk.utils.commandline.Input
import org.apache.commons.io.FileUtils
import java.io.File
import collection.JavaConversions._

/**
 * Concatenate log files to the jobOutputFile.
 */
class ConcatenateLogsFunction extends InProcessFunction {
  analysisName = "Concat"

  @Input(doc="Parts to gather back into the original output")
  var logs: Seq[File] = Nil

  override def description = "%s: %s > %s".format(analysisName, logs, jobOutputFile)
  override def shortDescription = analysisName + ": " + jobOutputFile.getName

  def run() {
    val missing = org.broadinstitute.gatk.utils.io.IOUtils.waitFor(logs, 120)
    if (!missing.isEmpty)
      throw new QException("Unable to find log: " + missing.mkString(", "))
    logs.foreach(log => {
      FileUtils.copyFile(log, this.jobOutputStream)
      this.jobOutputStream.println()
    })
  }
}
