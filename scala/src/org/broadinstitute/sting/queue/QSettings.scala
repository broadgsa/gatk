package org.broadinstitute.sting.queue

import org.broadinstitute.sting.commandline.Argument
import java.io.File
import java.lang.management.ManagementFactory

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
}

/**
 * Default settings settable on the command line and passed to CommandLineFunctions.
 */
object QSettings {
    /** A semi-unique job prefix using the host name and the process id. */
  private val processNamePrefix = "Q-" + {
    var prefix = ManagementFactory.getRuntimeMXBean.getName
    val index = prefix.indexOf(".")
    if (index >= 0)
      prefix = prefix.substring(0, index)
    prefix
  }
}
