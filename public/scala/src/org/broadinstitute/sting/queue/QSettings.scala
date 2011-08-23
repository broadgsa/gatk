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

package org.broadinstitute.sting.queue

import java.io.File
import org.broadinstitute.sting.commandline.{ArgumentCollection, Argument}
import org.broadinstitute.sting.queue.util.{SystemUtils, EmailSettings}

/**
 * Default settings settable on the command line and passed to CommandLineFunctions.
 */
class QSettings {
  @Argument(fullName="job_name_prefix", shortName="jobPrefix", doc="Default name prefix for compute farm jobs.", required=false)
  var jobNamePrefix: String = QSettings.processNamePrefix

  @Argument(fullName="job_project", shortName="jobProject", doc="Default project for compute farm jobs.", required=false)
  var jobProject: String = _

  @Argument(fullName="job_queue", shortName="jobQueue", doc="Default queue for compute farm jobs.", required=false)
  var jobQueue: String = _

  @Argument(fullName="job_priority", shortName="jobPriority", doc="Default priority for jobs. Min = 0, Max = 100", required=false)
  var jobPriority: Option[Int] = None

  @Argument(fullName="job_native_arg", shortName="jobNative", doc="Native arguments to pass to the job runner.", required=false)
  var jobNativeArgs: List[String] = Nil

  @Argument(fullName="job_resource_request", shortName="jobResReq", doc="Resource requests to pass to the job runner.", required=false)
  var jobResourceRequests: List[String] = Nil

  @Argument(fullName="job_environment_name", shortName="jobEnv", doc="Environment names for the job runner.", required=false)
  var jobEnvironmentNames: List[String] = Nil

  @Argument(fullName="memory_limit", shortName="memLimit", doc="Default memory limit for jobs, in gigabytes.", required=false)
  var memoryLimit: Option[Double] = None

  @Argument(fullName="resident_memory_limit", shortName="resMemLimit", doc="Default resident memory limit for jobs, in gigabytes.", required=false)
  var residentLimit: Option[Double] = None

  @Argument(fullName="resident_memory_request", shortName="resMemReq", doc="Default resident memory request for jobs, in gigabytes.", required=false)
  var residentRequest: Option[Double] = None

  @Argument(fullName="run_directory", shortName="runDir", doc="Root directory to run functions from.", required=false)
  var runDirectory = new File(".")

  @Argument(fullName="temp_directory", shortName="tempDir", doc="Temp directory to pass to functions.", required=false)
  var tempDirectory = new File(System.getProperty("java.io.tmpdir"))

  @Argument(fullName="job_scatter_gather_directory", shortName="jobSGDir", doc="Default directory to place scatter gather output for compute farm jobs.", required=false)
  var jobScatterGatherDirectory: File = _

  @ArgumentCollection
  val emailSettings = new EmailSettings
}

/**
 * Default settings settable on the command line and passed to CommandLineFunctions.
 */
object QSettings {
  /** A semi-unique job prefix using the host name and the process id. */
  private val processNamePrefix = "Q-" + SystemUtils.pidAtHost
}
