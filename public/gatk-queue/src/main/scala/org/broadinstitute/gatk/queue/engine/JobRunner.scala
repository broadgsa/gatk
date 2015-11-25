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

/**
 * Base interface for job runners.
 */
trait JobRunner[TFunction <: QFunction] {

  /**
   * Initializes this job.
   */
  def init() {
  }

  /**
   * Runs the function.
   * After the function returns the status of the function should
   * be RUNNING, FAILED, or DONE.
   * @param function Command to run.
   */
  def start()

  /**
   * Returns the current run status.
   * Must only be called AFTER start().
   * @return RUNNING, DONE, or FAILED.
   */
  def status: RunnerStatus.Value

  /**
   * Checks if the status has been unknown for an extended period of time.
   */
  def checkUnknownStatus() {}

  /**
   * Returns the function to be run.
   */
  def function: TFunction

  /**
   * Cleans up after the function is run.
   * For example removing all temporary files.
   */
  def cleanup() {
  }

  /**
   * Must be overloaded
   */
  val runInfo = JobRunInfo.default
  def getRunInfo = runInfo

  /**
   * Calls back to a hook that an expert user can setup to modify a job.
   * @param value Value to modify.
   */
  protected def updateJobRun(value: Any) {
    val updater = function.updateJobRun
    if (updater != null)
      if (updater.isDefinedAt(value))
        updater(value)
  }
}
