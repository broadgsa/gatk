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

  @Argument(fullName="job_queue", shortName="jobQueue", doc="Default queue for compute farm jobs.", required=false)
  var jobQueue: String = _

  @Argument(fullName="job_project", shortName="jobProject", doc="Default project for compute farm jobs.", required=false)
  var jobProject: String = "Queue"

  @Argument(fullName="job_scatter_gather_directory", shortName="jobSGDir", doc="Default directory to place scatter gather output for compute farm jobs.", required=false)
  var jobScatterGatherDirectory: File = _

  @Argument(fullName="default_memory_limit", shortName="memLimit", doc="Default memory limit for jobs, in gigabytes.", required=false)
  var memoryLimit: Option[Int] = None

  @Argument(fullName="run_directory", shortName="runDir", doc="Root directory to run functions from.", required=false)
  var runDirectory = new File(".")

  @Argument(fullName="temp_directory", shortName="tempDir", doc="Temp directory to pass to functions.", required=false)
  var tempDirectory = new File(System.getProperty("java.io.tmpdir"))

  @Argument(fullName="mount_directory", shortName="mountDir", doc="Extra directory to automount via 'cd <dir>' before running functions.", required=false)
  var mountDirectories: Set[File] = Set.empty[File]

  @ArgumentCollection
  val emailSettings = new EmailSettings
}

/**
 * Default settings settable on the command line and passed to CommandLineFunctions.
 */
object QSettings {
  /** A semi-unique job prefix using the host name and the process id. */
  private val processNamePrefix = "Q-" + SystemUtils.pidAtHost.split('.')(0)
}
