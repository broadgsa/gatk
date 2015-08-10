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

package org.broadinstitute.gatk.queue.util

import java.io.{PrintStream, File}
import org.broadinstitute.gatk.utils.io.{Resource}
import org.broadinstitute.gatk.queue.engine.{JobRunInfo, QGraph}
import org.broadinstitute.gatk.queue.function.QFunction
import org.broadinstitute.gatk.utils.R.{RScriptLibrary, RScriptExecutor}
import org.broadinstitute.gatk.utils.report.GATKReportTable
import org.broadinstitute.gatk.utils.exceptions.UserException
import org.apache.commons.io.{FileUtils, IOUtils}
import org.broadinstitute.gatk.utils.report.{GATKReportTable, GATKReport}

/**
 * Writes out RunInfo to a GATKReport
 */
class QJobsReporter(val disabled: Boolean, val reportFile: File, val pdfFile: Option[File]) extends Logging {
  private val JOB_REPORT_QUEUE_SCRIPT = "queueJobReport.R"

  /**
   * Write out a job report based on the finished jobs graph
   * @param jobGraph
   * @param enabledPlotting if true, we will plot the report as well with the JOB_REPORT_QUEUE_SCRIPT
   */
  def write(jobGraph: QGraph, enabledPlotting: Boolean) {
    if ( ! disabled ) {
      logger.info("Writing JobLogging GATKReport to file " + reportFile)
      printReport(jobGraph.getFunctionsAndStatus, reportFile)

      if ( enabledPlotting )
        pdfFile match {
          case Some(file) =>
            logger.info("Plotting JobLogging GATKReport to file " + file)
            plotReport(reportFile, file)
          case None =>
        }
    }
  }

  private def printReport(jobsRaw: Map[QFunction, JobRunInfo], dest: File) {
    val jobs = jobsRaw.filter(_._2.isFilledIn).filter(_._1.includeInReport)
    jobs foreach {case (qf, info) => qf.setRunInfo(info)}
    val stream = new PrintStream(FileUtils.openOutputStream(dest))
    try {
      printJobLogging(jobs.keys.toSeq, stream)
    } finally {
      IOUtils.closeQuietly(stream)
    }
  }

  private def plotReport(reportFile: File, pdfFile: File) {
    val executor = new RScriptExecutor
    executor.addLibrary(RScriptLibrary.GSALIB)
    executor.addScript(new Resource(JOB_REPORT_QUEUE_SCRIPT, classOf[QJobReport]))
    executor.addArgs(reportFile.getAbsolutePath, pdfFile.getAbsolutePath)
    executor.exec()
  }

  /**
   * Prints the JobLogging logs to a GATKReport.  First splits up the
   * logs by group, and for each group generates a GATKReportTable
   */
  private def printJobLogging(logs: Seq[QFunction], stream: PrintStream) {
    // create the report
    val report: GATKReport = new GATKReport

    // create a table for each group of logs
    for ( (group, groupLogs) <- groupLogs(logs) ) {
      val keys = logKeys(groupLogs)
      report.addTable(group, "Job logs for " + group, keys.size)
      val table: GATKReportTable = report.getTable(group)

      // add the columns
      keys.foreach(table.addColumn(_))
      for (log <- groupLogs) {
        for ( key <- keys )
          table.set(log.getReportName, key, log.getReportFeature(key))
      }
    }

    report.print(stream)
  }

  private def groupLogs(logs: Seq[QFunction]): Map[String, Seq[QFunction]] = {
    logs.groupBy(_.getReportGroup)
  }

  private def logKeys(logs: Seq[QFunction]): Set[String] = {
    // the keys should be the same for each log, but we will check that
    val keys = Set[String](logs(0).getReportFeatureNames : _*)

    for ( log <- logs )
      if ( keys.sameElements(Set(log.getReportFeatureNames)) )
        throw new UserException(("All JobLogging jobs in the same group must have the same set of features.  " +
          "We found one with %s and another with %s").format(keys, log.getReportFeatureNames))

    keys
  }
}
