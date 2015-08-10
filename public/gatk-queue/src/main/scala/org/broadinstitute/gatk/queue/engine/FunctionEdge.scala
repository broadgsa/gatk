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

package org.broadinstitute.gatk.queue.engine

import org.broadinstitute.gatk.queue.function.QFunction
import java.io.{StringWriter, PrintWriter}
import org.broadinstitute.gatk.queue.util.Logging
import org.broadinstitute.gatk.utils.io.IOUtils
import org.apache.commons.io.FileUtils
import org.apache.commons.lang.StringUtils

/**
 * An edge in the QGraph that runs a QFunction.
 * The edge is created first to determine inter-node dependencies,
 * and then the runner is specified later when the time comes to
 * execute the function in the edge.
 */
class FunctionEdge(val function: QFunction, val inputs: QNode, val outputs: QNode) extends QEdge with Logging {
  var runner: JobRunner[_] =_

  /**
   * The depth of this edge in the graph.
   */
  var depth = -1

  val myRunInfo: JobRunInfo = JobRunInfo.default // purely for dryRun testing

  /**
   * When using reset status this variable tracks the old status
   */
  var resetFromStatus: RunnerStatus.Value = null

  /**
   * Initializes with the current status of the function.
   */
  private var currentStatus = {
    if (function.isFail)
      RunnerStatus.FAILED
    else if (function.isDone)
      RunnerStatus.DONE
    else
      RunnerStatus.PENDING
  }

  def start() {
    try {
      if (logger.isDebugEnabled) {
        logger.debug("Starting: " + function.commandDirectory + " > " + function.description)
      } else {
        logger.info("Starting: " + function.description)
      }
      logger.info("Output written to " + function.jobOutputFile)
      if (function.jobErrorFile != null)
        logger.info("Errors written to " + function.jobErrorFile)

      function.deleteLogs()
      function.deleteOutputs()
      function.mkOutputDirectories()

      runner.init()
      runner.start()
    } catch {
      case e: Throwable =>
        currentStatus = RunnerStatus.FAILED
        try {
          runner.cleanup()
          function.failOutputs.foreach(_.createNewFile())
          writeStackTrace(e)
        } catch {
          case _: Throwable => /* ignore errors in the exception handler */
        }
        logger.error("Error: " + function.description, e)
    }
  }

  /**
   * Returns the current status of the edge.
   */
  def status = {
    if (currentStatus == RunnerStatus.PENDING || currentStatus == RunnerStatus.RUNNING) {
      if (runner != null) {
        try {
          currentStatus = runner.status

          if (currentStatus == RunnerStatus.FAILED) {
            try {
              runner.cleanup()
              function.failOutputs.foreach(_.createNewFile())
            } catch {
              case _: Throwable => /* ignore errors in the error handler */
            }
            logger.error("Error: " + function.description)
            tailError()
          } else if (currentStatus == RunnerStatus.DONE) {
            try {
              runner.cleanup()
              function.doneOutputs.foreach(_.createNewFile())
            } catch {
              case _: Throwable => /* ignore errors in the done handler */
            }
            logger.info("Done: " + function.description)
          }
        } catch {
          case e: Throwable =>
            currentStatus = RunnerStatus.FAILED
            try {
              runner.cleanup()
              function.failOutputs.foreach(_.createNewFile())
              writeStackTrace(e)
            } catch {
              case _: Throwable => /* ignore errors in the exception handler */
            }
            logger.error("Error retrieving status: " + function.description, e)
        }
      }
    }

    currentStatus
  }

  /**
   * Explicitly sets the status of the runner to done..
   */
  def markAsDone() {
    currentStatus = RunnerStatus.DONE
  }

  /**
   * Marks this edge as skipped as it is not needed for the current run.
   */
  def markAsSkipped() {
    currentStatus = RunnerStatus.SKIPPED
  }

  /**
   * Resets the edge to pending status.
   */
  def resetToPending(cleanOutputs: Boolean) {
    if (resetFromStatus == null)
      resetFromStatus = currentStatus
    currentStatus = RunnerStatus.PENDING
    if (cleanOutputs)
      function.deleteOutputs()
    function.jobErrorLines = Nil
    runner = null
  }

  override def shortDescription = function.shortDescription

  /**
   * Returns the path to the file to use for logging errors.
   * @return the path to the file to use for logging errors.
   */
  private def functionErrorFile = if (function.jobErrorFile != null) function.jobErrorFile else function.jobOutputFile

  /**
   * Outputs the last lines of the error logs.
   */
  private def tailError() {
    val errorFile = functionErrorFile
    if (IOUtils.waitFor(errorFile, 120)) {
      val maxLines = 100
      val tailLines = IOUtils.tail(errorFile, maxLines)
      val nl = "%n".format()
      val summary = if (tailLines.size > maxLines) "Last %d lines".format(maxLines) else "Contents"
      this.function.jobErrorLines = collection.JavaConversions.collectionAsScalaIterable(tailLines).toSeq
      logger.error("%s of %s:%n%s".format(summary, errorFile, StringUtils.join(tailLines, nl)))
    } else {
      logger.error("Unable to access log file: %s".format(errorFile))
    }
  }

  /**
   * Writes the stack trace to the error file.
   */
  private def writeStackTrace(e: Throwable) {
    val stackTrace = new StringWriter
    val printWriter = new PrintWriter(stackTrace)
    printWriter.println(function.description)
    e.printStackTrace(printWriter)
    printWriter.close()
    FileUtils.writeStringToFile(functionErrorFile, stackTrace.toString)
  }

  def getRunInfo = {
    if ( runner == null ) myRunInfo else runner.getRunInfo
  }
}
