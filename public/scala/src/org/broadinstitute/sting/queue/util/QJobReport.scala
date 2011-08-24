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

/**
 * A mixin to add Job info to the class
 */
trait QJobReport extends QFunction {
  private var group: String = _
  private var features: Map[String, String] = null

  def getGroup = group
  def isEnabled = group != null
  def getFeatureNames: List[String] = features.keys.toList
  def getFeatures = features
  def get(key: String): String = {
    features.get(key) match {
      case Some(x) => x
      case None => throw new RuntimeException("Get called with key %s but no value was found".format(key))
    }
  }
  def getName: String = features.get("jobName").get

  private def addRunInfo(info: JobRunInfo) {
    features = features ++ Map(
      "analysisName" -> this.analysisName,
      "jobName" -> this.jobName,
      "intermediate" -> this.isIntermediate,
      "startTime" -> info.getStartTime,
      "doneTime" -> info.getDoneTime,
      "memUsedInGb" -> info.getMemoryUsedInGb,
      "runtime" -> info.getRuntimeInMs,
      "hostName" -> info.getHostname).mapValues(_.toString)
  }

  def setJobLogging(group: String) {
    this.group = group
  }

  def setJobLogging(group: String, features: Map[String, Any]) {
    this.group = group
    this.features = features.mapValues(_.toString)
  }
}

object QJobReport {
  def printReport(jobs: Map[QFunction, JobRunInfo], dest: File) {
    val jobLogs: List[QJobReport] = jobLoggingSublist(jobs.keys.toList)
    jobLogs.foreach((job: QJobReport) => job.addRunInfo(jobs.get(job).get))
    val stream = new PrintStream(new FileOutputStream(dest))
    printJobLogging(jobLogs, stream)
    stream.close()
  }

  private def jobLoggingSublist(l: List[QFunction]): List[QJobReport] = {
    def asJogLogging(qf: QFunction): QJobReport = qf match {
      case x: QJobReport => x
      case _ => null
    }

    l.map(asJogLogging).filter(_ != null)
  }

  /**
   * Prints the JobLogging logs to a GATKReport.  First splits up the
   * logs by group, and for each group generates a GATKReportTable
   */
  private def printJobLogging(logs: List[QJobReport], stream: PrintStream) {
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
          table.set(log.getName, key, log.get(key))
      }
    }

    report.print(stream)
  }

  private def groupLogs(logs: List[QJobReport]): Map[String, List[QJobReport]] = {
    logs.groupBy(_.getGroup)
  }

  private def logKeys(logs: List[QJobReport]): Set[String] = {
    // the keys should be the same for each log, but we will check that
    val keys = Set[String](logs(0).getFeatureNames : _*)

    for ( log <- logs )
      if ( keys.sameElements(Set(log.getFeatureNames)) )
        throw new UserException(("All JobLogging jobs in the same group must have the same set of features.  " +
          "We found one with %s and another with %s").format(keys, log.getFeatureNames))

    keys
  }
}
