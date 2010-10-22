package org.broadinstitute.sting.queue.engine

import java.io.File
import org.broadinstitute.sting.queue.QSettings
import org.broadinstitute.sting.commandline.{ArgumentCollection, Argument}
import org.broadinstitute.sting.queue.util.SystemUtils

/**
 * Command line options for a QGraph.
 */
class QGraphSettings {
  @ArgumentCollection
  val qSettings = new QSettings

  @Argument(fullName="bsub_all_jobs", shortName="bsub", doc="Use bsub to submit jobs", required=false)
  var bsubAllJobs = false

  @Argument(fullName="run_scripts", shortName="run", doc="Run QScripts.  Without this flag set only performs a dry run.", required=false)
  var run = false

  @Argument(fullName="dot_graph", shortName="dot", doc="Outputs the queue graph to a .dot file.  See: http://en.wikipedia.org/wiki/DOT_language", required=false)
  var dotFile: File = _

  @Argument(fullName="expanded_dot_graph", shortName="expandedDot", doc="Outputs the queue graph of scatter gather to a .dot file.  Otherwise overwrites the dot_graph", required=false)
  var expandedDotFile: File = _

  @Argument(fullName="start_from_scratch", shortName="startFromScratch", doc="Runs all command line functions even if the outputs were previously output successfully.", required=false)
  var startFromScratch = false

  @Argument(fullName="status",shortName="status",doc="Get status of jobs for the qscript",required=false)
  var getStatus = false

  @Argument(fullName="status_email_from", shortName="statusFrom", doc="Email address to send emails from upon completion or on error.", required=false)
  var statusEmailFrom: String = System.getProperty("user.name") + "@" + SystemUtils.domainName

  @Argument(fullName="status_email_to", shortName="statusTo", doc="Email address to send emails to upon completion or on error.", required=false)
  var statusEmailTo: List[String] = Nil

  @Argument(fullName="delete_intermediate_outputs", shortName="deleteIntermediates", doc="After a successful run delete the outputs of any Function marked as intermediate.", required=false)
  var deleteIntermediates = false

  @Argument(fullName="retry_failed", shortName="retry", doc="Retry the specified number of times after a command fails.  Defaults to no retries.", required=false)
  var retries = 0
}
