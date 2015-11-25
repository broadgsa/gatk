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

package org.broadinstitute.gatk.queue.engine.parallelshell

import java.io.PrintWriter

import org.broadinstitute.gatk.queue.util.Logging
import org.broadinstitute.gatk.utils.runtime.ProcessSettings
import scala.sys.process._

/**
 *
 */
class ThreadSafeProcessController extends Logging {

  private var process: Option[Process] = None

  /**
   * Construct a process logger writing the stdout and stderr of the
   * process controlled by this instance to the files specified in
   * the provided ProcessSettings instance.
   * @param processSettings specifiying which files to write to
   * @return a process logger which can be used by the `scala.sys.process`
   */
  private def getProcessLogger(processSettings: ProcessSettings): ProcessLogger = {

    val (stdOutFile, stdErrFile) = {

      val stdOutFile = processSettings.getStdoutSettings.getOutputFile

      if(processSettings.getStderrSettings.getOutputFile != null) {
        val stdErrFile = processSettings.getStderrSettings.getOutputFile
        (stdOutFile, stdErrFile)
      } else {
        (stdOutFile, stdOutFile)
      }

    }

    val stdOutPrintWriter = new PrintWriter(stdOutFile)
    val stdErrPrintWriter = new PrintWriter(stdErrFile)

    def printToWriter(printWriter: PrintWriter)(line: String): Unit = {
      printWriter.println(line)
      printWriter.flush()
    }

    val stringStdOutPrinterFunc = printToWriter(stdOutPrintWriter) _
    val stringStdErrPrinterFunc = printToWriter(stdErrPrintWriter) _

    val processLogger = ProcessLogger(
      stringStdOutPrinterFunc,
      stringStdErrPrinterFunc
      )

    processLogger
  }

  /**
   * Execute the process specified in process settings
   * @param processSettings specifying the commandline to run.
   * @return the exit status of the process.
   */
  def exec(processSettings: ProcessSettings): Int = {

    val commandLine: ProcessBuilder = processSettings.getCommand.mkString(" ")
    logger.debug("Trying to start process: " + commandLine)
    process = Some(commandLine.run(getProcessLogger(processSettings)))
    process.get.exitValue()

  }

  /**
   * Attempt to destroy the underlying process.
   */
  def tryDestroy(): Unit = {
    logger.debug("Trying to kill process")
    process.getOrElse {
      throw new IllegalStateException("Tried to kill unstarted job.")
    }.destroy()
  }

}
