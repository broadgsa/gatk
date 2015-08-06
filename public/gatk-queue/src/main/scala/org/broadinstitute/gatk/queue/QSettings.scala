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

package org.broadinstitute.gatk.queue

import java.io.File
import org.broadinstitute.gatk.utils.commandline.{ClassType, Argument}

/**
 * Default settings settable on the command line and passed to CommandLineFunctions.
 */
class QSettings {
  @Argument(fullName="run_name", shortName="runName", doc="A name for this run used for various status messages.", required=false)
  var runName: String = _

  @Argument(fullName="job_project", shortName="jobProject", doc="Default project for compute farm jobs.", required=false)
  var jobProject: String = _

  @Argument(fullName="job_queue", shortName="jobQueue", doc="Default queue for compute farm jobs.", required=false)
  var jobQueue: String = _

  @Argument(fullName="job_priority", shortName="jobPriority", doc="Default priority for jobs. Min = 0, Max = 100", required=false)
  @ClassType(classOf[Int])
  var jobPriority: Option[Int] = None

  @Argument(fullName="job_native_arg", shortName="jobNative", doc="Native arguments to pass to the job runner.", required=false)
  var jobNativeArgs: Seq[String] = Nil

  @Argument(fullName="job_resource_request", shortName="jobResReq", doc="Resource requests to pass to the job runner.", required=false)
  var jobResourceRequests: Seq[String] = Nil

  @Argument(fullName="job_environment_name", shortName="jobEnv", doc="Environment names for the job runner.", required=false)
  var jobEnvironmentNames: Seq[String] = Nil

  @Argument(fullName="memory_limit", shortName="memLimit", doc="Default memory limit for jobs, in gigabytes. If not set defaults to 2GB.", required=false)
  @ClassType(classOf[Double])
  var memoryLimit: Option[Double] = Some(2)

  @Argument(fullName="memory_limit_threshold", shortName="memLimitThresh", doc="After passing this threshold stop increasing memory limit for jobs, in gigabytes.", required=false)
  @ClassType(classOf[Double])
  var memoryLimitThreshold: Option[Double] = None

  @Argument(fullName="resident_memory_limit", shortName="resMemLimit", doc="Default resident memory limit for jobs, in gigabytes.", required=false)
  @ClassType(classOf[Double])
  var residentLimit: Option[Double] = None

  @Argument(fullName="resident_memory_request", shortName="resMemReq", doc="Default resident memory request for jobs, in gigabytes.", required=false)
  @ClassType(classOf[Double])
  var residentRequest: Option[Double] = None

  @Argument(fullName="resident_memory_request_parameter", shortName="resMemReqParam", doc="Parameter for resident memory requests. By default not requested.", required=false)
  var residentRequestParameter: String = _

  @Argument(fullName="job_walltime", shortName="wallTime", doc="Setting the required DRMAA walltime or LSF run limit.", required=false)
  @ClassType(classOf[Long])
  var jobWalltime: Option[Long] = None

  /** The name of the parallel environment (required for SGE, for example) */
  @Argument(fullName="job_parallel_env", shortName="jobParaEnv", doc="An SGE style parallel environment to use for jobs requesting more than 1 core.  Equivalent to submitting jobs with -pe ARG nt for jobs with nt > 1", required=false)
  var parallelEnvironmentName: String = "smp_pe" // Broad default

  @Argument(fullName="dontRequestMultipleCores", shortName="multiCoreJerk", doc="If provided, Queue will not request multiple processors for jobs using multiple processors.  Sometimes you eat the bear, sometimes the bear eats you.", required=false)
  var dontRequestMultipleCores: Boolean = false

  @Argument(fullName="disableDefaultJavaGCOptimizations", shortName="noGCOpt", doc="If provided, Queue will not ensure that java GC threads are limited and that the a minimum amount of time is spent in GC.")
  var disableDefaultJavaGCOptimizations = false

  @Argument(fullName="run_directory", shortName="runDir", doc="Root directory to run functions from.", required=false)
  var runDirectory = new File(".")

  @Argument(fullName="temp_directory", shortName="tempDir", doc="Temp directory to pass to functions.", required=false)
  var tempDirectory = new File(System.getProperty("java.io.tmpdir"))

  @Argument(fullName="job_scatter_gather_directory", shortName="jobSGDir", doc="Default directory to place scatter gather output for compute farm jobs.", required=false)
  var jobScatterGatherDirectory: File = _
  
  @Argument(fullName="log_directory", shortName="logDir", doc="Directory to write log files into.", required=false)
  var logDirectory: File = _

  /**
   * If set, use Broad-specific cluster settings in the GridEngine job runner. Activated via the -qsub-broad argument in QGraphSettings.
   */
  var useBroadClusterSettings: Boolean = false
}
