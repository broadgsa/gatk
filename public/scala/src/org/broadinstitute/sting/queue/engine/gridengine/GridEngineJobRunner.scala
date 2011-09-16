/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.engine.gridengine

import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.engine.drmaa.DrmaaJobRunner
import org.ggf.drmaa.Session

/**
 * Runs jobs on a Grid Engine compute cluster.
 */
class GridEngineJobRunner(session: Session, function: CommandLineFunction) extends DrmaaJobRunner(session, function) with Logging {
  // Grid Engine disallows certain characters from being in job names.
  // This replaces all illegal characters with underscores
  protected override val jobNameFilter = """[\n\t\r/:@\\*?]"""
  protected override val minRunnerPriority = -1023
  protected override val maxRunnerPriority = 0

  override protected def functionNativeSpec = {
    // Force the remote environment to inherit local environment settings
    var nativeSpec: String = "-V"

    // If a project name is set specify the project name
    if (function.jobProject != null)
      nativeSpec += " -P " + function.jobProject

    // If the job queue is set specify the job queue
    if (function.jobQueue != null)
      nativeSpec += " -q " + function.jobQueue

    // If the resident set size is requested pass on the memory request
    if (function.residentRequest.isDefined)
      nativeSpec += " -l mem_free=%dM".format(function.residentRequest.map(_ * 1024).get.ceil.toInt)

    // If the resident set size limit is defined specify the memory limit
    if (function.residentLimit.isDefined)
      nativeSpec += " -l h_rss=%dM".format(function.residentLimit.map(_ * 1024).get.ceil.toInt)

    // Pass on any job resource requests
    nativeSpec += function.jobResourceRequests.map(" -l " + _).mkString

    // Pass on any job environment names
    nativeSpec += function.jobEnvironmentNames.map(" -pe " + _).mkString

    // If the priority is set specify the priority
    val priority = functionPriority
    if (priority.isDefined)
      nativeSpec += " -p " + priority.get

    (nativeSpec + " " + super.functionNativeSpec).trim()
  }
}
