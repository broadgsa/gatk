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