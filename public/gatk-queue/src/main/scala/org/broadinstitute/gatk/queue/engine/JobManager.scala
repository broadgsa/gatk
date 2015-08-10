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
 * Creates and stops JobRunners
 */
trait JobManager[TFunction <: QFunction, TRunner <: JobRunner[TFunction]] {
  def init() {}
  def exit() {}

  /** The class type of the runner.  Available at runtime even after erasure. */
  def functionType: Class[TFunction]

  /** The class type of the functions processed by the runner.  Available at runtime even after erasure. */
  def runnerType: Class[TRunner]

  /** Creates a new runner.
   * @param function Function for the runner.
   */
  def create(function: TFunction): TRunner

  /**
   * Updates the status on a list of functions.
   * @param runners Runners to update.
   * @return runners which were updated.
   */
  def updateStatus(runners: Set[TRunner]): Set[TRunner] = Set.empty

  /**
   * Stops a list of functions.
   * @param runners Runners to stop.
   */
  def tryStop(runners: Set[TRunner]) {}
}
