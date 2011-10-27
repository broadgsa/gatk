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

package org.broadinstitute.sting.queue.engine

import java.io.File
import org.broadinstitute.sting.queue.QSettings
import org.broadinstitute.sting.queue.util.SystemUtils
import org.broadinstitute.sting.commandline.{Advanced, ArgumentCollection, Argument}

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

  @Argument(fullName="status",shortName="status",doc="Get status of jobs for the qscript",required=false)
  var getStatus = false

  @Argument(fullName="retry_failed", shortName="retry", doc="Retry the specified number of times after a command fails.  Defaults to no retries.", required=false)
  var retries = 0

  @Argument(fullName="start_from_scratch", shortName="startFromScratch", doc="Runs all command line functions even if the outputs were previously output successfully.", required=false)
  var startFromScratch = false

  @Argument(fullName="keep_intermediate_outputs", shortName="keepIntermediates", doc="After a successful run keep the outputs of any Function marked as intermediate.", required=false)
  var keepIntermediates = false

  @Argument(fullName="status_email_to", shortName="statusTo", doc="Email address to send emails to upon completion or on error.", required=false)
  var statusEmailTo: List[String] = Nil

  @Argument(fullName="status_email_from", shortName="statusFrom", doc="Email address to send emails from upon completion or on error.", required=false)
  var statusEmailFrom: String = System.getProperty("user.name") + "@" + SystemUtils.mailName

  @Argument(fullName="dot_graph", shortName="dot", doc="Outputs the queue graph to a .dot file.  See: http://en.wikipedia.org/wiki/DOT_language", required=false)
  var dotFile: File = _

  @Argument(fullName="expanded_dot_graph", shortName="expandedDot", doc="Outputs the queue graph of scatter gather to a .dot file.  Otherwise overwrites the dot_graph", required=false)
  var expandedDotFile: File = _

  @Argument(fullName="jobReport", shortName="jobReport", doc="File where we will write the Queue job report", required=false)
  var jobReportFile: String = _

  @Advanced
  @Argument(fullName="disableJobReport", shortName="disabpleJobReport", doc="If provided, we will not create a job report", required=false)
  var disableJobReport: Boolean = false

  @ArgumentCollection
  val qSettings = new QSettings
}
