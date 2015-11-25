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

import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.engine.CommandLineJobManager

/**
 * Runs multiple jobs locally without blocking.
 * Use this with care as it might not be the most efficient way to run things.
 * However, for some scenarios, such as running multiple single threaded
 * programs concurrently it can be quite useful.
 * 
 * All this code is based on the normal shell runner in GATK Queue and all 
 * credits for everything except the concurrency part goes to the GATK team.
 * 
 * @author Johan Dahlberg
 *
 */
class ParallelShellJobManager extends CommandLineJobManager[ParallelShellJobRunner] {

  def runnerType = classOf[ParallelShellJobRunner]

  /**
   * Create new ParallelShellJobRunner
   * @param function Function for the runner.
   * @return a new ParallelShellJobRunner instance
   */
  def create(function: CommandLineFunction) =
    new ParallelShellJobRunner(function)

  /**
   * Update the status of the specified jobrunners.
   * @param runners Runners to update.
   * @return runners which were updated.
   */
  override def updateStatus(
    runners: Set[ParallelShellJobRunner]): Set[ParallelShellJobRunner] =
    runners.filter { runner => runner.updateJobStatus() }

  /**
   * Stop the specified runners.
   * @param runners Runners to stop.
   */
  override def tryStop(runners: Set[ParallelShellJobRunner]) =
    runners.foreach(_.tryStop())
}
