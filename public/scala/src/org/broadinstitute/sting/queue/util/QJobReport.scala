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

package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.gatk.report.{GATKReportTable, GATKReport}
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.queue.engine.JobRunInfo
import java.io.{FileOutputStream, PrintStream, File}
import org.broadinstitute.sting.utils.R.{RScriptLibrary, RScriptExecutor}
import org.broadinstitute.sting.utils.io.Resource

/**
 * A mixin to add Job info to the class
 */
trait QJobReport extends Logging {
  self: QFunction =>

  protected var reportGroup: String = null
  protected var reportFeatures: Map[String, String] = Map()
  protected var reportEnabled: Boolean = true

  def includeInReport = reportEnabled
  def enableReport() { reportEnabled = true }
  def disableReport() { reportEnabled = false }

  def setRunInfo(info: JobRunInfo) {
    //logger.info("info " + info)
    reportFeatures = Map(
      "iteration" -> 1,
      "analysisName" -> getReportGroup,
      "jobName" -> QJobReport.workAroundSameJobNames(this),
      "intermediate" -> self.isIntermediate,
      "exechosts" -> info.getExecHosts,
      "startTime" -> info.getStartTime.getTime,
      "doneTime" -> info.getDoneTime.getTime,
      "formattedStartTime" -> info.getFormattedStartTime,
      "formattedDoneTime" -> info.getFormattedDoneTime,
      "runtime" -> info.getRuntimeInMs).mapValues((x:Any) => if (x != null) x.toString else "null") ++ reportFeatures
    // note -- by adding reportFeatures second we override iteration
    // (or any other binding) with the user provided value
  }

  /** The report Group is the analysis name transform to only contain valid GATKReportTable characters */
  def getReportGroup = self.analysisName.replaceAll(GATKReportTable.INVALID_TABLE_NAME_REGEX, "_")
  def getReportFeatures = reportFeatures

  def getReportFeatureNames: List[String] = getReportFeatures.keys.toList
  def getReportFeature(key: String): String = {
    getReportFeatures.get(key) match {
      case Some(x) => x
      case None => throw new RuntimeException("Get called with key %s but no value was found".format(key))
    }
  }

  def getReportName: String = getReportFeature("jobName")

  def configureJobReport(features: Map[String, Any]) {
    this.reportFeatures = features.mapValues(_.toString)
  }

  // copy the QJobReport information -- todo : what's the best way to do this?
  override def copySettingsTo(function: QFunction) {
    self.copySettingsTo(function)
    function.reportFeatures = this.reportFeatures
  }
}

object QJobReport {
  val JOB_REPORT_QUEUE_SCRIPT = "queueJobReport.R"

  // todo -- fixme to have a unique name for Scatter/gather jobs as well
  var seenCounter = 1
  var seenNames = Set[String]()

  def printReport(jobsRaw: Map[QFunction, JobRunInfo], dest: File) {
    val jobs = jobsRaw.filter(_._2.isFilledIn).filter(_._1.includeInReport)
    jobs foreach {case (qf, info) => qf.setRunInfo(info)}
    val stream = new PrintStream(new FileOutputStream(dest))
    printJobLogging(jobs.keys.toList, stream)
    stream.close()
  }

  def plotReport(reportFile: File, pdfFile: File) {
    val executor = new RScriptExecutor
    executor.addLibrary(RScriptLibrary.GSALIB)
    executor.addScript(new Resource(JOB_REPORT_QUEUE_SCRIPT, classOf[QJobReport]))
    executor.addArgs(reportFile.getAbsolutePath, pdfFile.getAbsolutePath)
    executor.exec()
  }

  def workAroundSameJobNames(func: QFunction):String = {
    if ( seenNames.apply(func.jobName) ) {
      seenCounter += 1
      "%s_%d".format(func.jobName, seenCounter)
    } else {
      seenNames += func.jobName
      func.jobName
    }
  }

  /**
   * Prints the JobLogging logs to a GATKReport.  First splits up the
   * logs by group, and for each group generates a GATKReportTable
   */
  private def printJobLogging(logs: List[QFunction], stream: PrintStream) {
    // create the report
    val report: GATKReport = new GATKReport

    // create a table for each group of logs
    for ( (group, groupLogs) <- groupLogs(logs) ) {
      report.addTable(group, "Job logs for " + group)
      val table: GATKReportTable = report.getTable(group)
      table.addPrimaryKey("jobName", false)
      val keys = logKeys(groupLogs)

      // add the columns
      keys.foreach(table.addColumn(_, 0))
      for (log <- groupLogs) {
        for ( key <- keys )
          table.set(log.getReportName, key, log.getReportFeature(key))
      }
    }

    report.print(stream)
  }

  private def groupLogs(logs: List[QFunction]): Map[String, List[QFunction]] = {
    logs.groupBy(_.getReportGroup)
  }

  private def logKeys(logs: List[QFunction]): Set[String] = {
    // the keys should be the same for each log, but we will check that
    val keys = Set[String](logs(0).getReportFeatureNames : _*)

    for ( log <- logs )
      if ( keys.sameElements(Set(log.getReportFeatureNames)) )
        throw new UserException(("All JobLogging jobs in the same group must have the same set of features.  " +
          "We found one with %s and another with %s").format(keys, log.getReportFeatureNames))

    keys
  }
}
