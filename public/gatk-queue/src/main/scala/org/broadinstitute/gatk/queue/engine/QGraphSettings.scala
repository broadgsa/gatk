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

import java.io.File
import org.broadinstitute.gatk.queue.QSettings
import org.broadinstitute.gatk.queue.util.{EmailSettings, SystemUtils}
import org.broadinstitute.gatk.utils.commandline.{ClassType, Advanced, ArgumentCollection, Argument}

/**
 * Command line options for a QGraph.
 */
class QGraphSettings {
  @Argument(fullName="run_scripts", shortName="run", doc="Run QScripts.  Without this flag set only performs a dry run.", required=false)
  var run = false

  @Argument(fullName="job_runner", shortName="jobRunner", doc="Use the specified job runner to dispatch command line jobs", required=false)
  var jobRunner: String = _

  @Argument(fullName="bsub", shortName="bsub", doc="Equivalent to -jobRunner Lsf706", required=false)
  var bsub = false

  @Argument(fullName="qsub", shortName="qsub", doc="Equivalent to -jobRunner GridEngine", required=false)
  var qsub = false

  @Argument(fullName="qsub-broad", shortName="qsub-broad", doc="Equivalent to -qsub, but uses GridEngine parameters specific to the Broad GridEngine cluster", required=false)
  var qsubBroad = false

  @Argument(fullName="status",shortName="status",doc="Get status of jobs for the qscript",required=false)
  var getStatus = false

  @Argument(fullName="retry_failed", shortName="retry", doc="Retry the specified number of times after a command fails.  Defaults to no retries.", required=false)
  var retries = 0

  @Argument(fullName="start_from_scratch", shortName="startFromScratch", doc="Runs all command line functions even if the outputs were previously output successfully.", required=false)
  var startFromScratch = false

  @Argument(fullName="keep_intermediate_outputs", shortName="keepIntermediates", doc="After a successful run keep the outputs of any Function marked as intermediate.", required=false)
  var keepIntermediates = false

  @Argument(fullName="status_email_to", shortName="statusTo", doc="Email address to send emails to upon completion or on error.", required=false)
  var statusEmailTo: Seq[String] = Nil

  @Argument(fullName="status_email_from", shortName="statusFrom", doc="Email address to send emails from upon completion or on error.", required=false)
  var statusEmailFrom: String = System.getProperty("user.name") + "@" + SystemUtils.mailName

  @Argument(fullName="graphviz", shortName="gv", doc="Outputs the queue graph to a Graphviz .gv file. See: http://www.graphviz.org/Documentation.php", required=false)
  var graphvizFile: File = _

  @Argument(fullName="graphviz_scatter_gather", shortName="gvsg", doc="Outputs the scatter/gather queue graph to a Graphviz .gv file.  Otherwise overwrites the --graphviz file.", required=false)
  var graphvizScatterGatherFile: File = _

  @Argument(fullName="jobReport", shortName="jobReport", doc="File where we will write the Queue job report", required=false)
  var jobReportFile: String = _

  @Advanced
  @Argument(fullName="disableJobReport", shortName="disableJobReport", doc="If provided, we will not create a job report", required=false)
  var disableJobReport: Boolean = false

  @Advanced
  @ClassType(classOf[Int])
  @Argument(fullName="maximumNumberOfJobsToRunConcurrently", shortName="maxConcurrentRun", doc="The maximum number of jobs to start at any given time. (Default is no limit)", required=false)
  var maximumNumberOfConcurrentJobs: Option[Int] = None

  @ArgumentCollection
  val emailSettings = new EmailSettings

  @ArgumentCollection
  val qSettings = new QSettings
}
