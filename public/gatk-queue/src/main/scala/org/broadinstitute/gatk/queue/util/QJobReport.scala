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

import org.broadinstitute.gatk.queue.function.QFunction
import org.broadinstitute.gatk.queue.engine.JobRunInfo
import org.broadinstitute.gatk.utils.report.GATKReportTable

/**
 * A mixin to add Job info to the class
 */
trait QJobReport extends Logging {
  self: QFunction =>

  protected var reportGroup: String = null
  protected var userReportFeatures: Map[String, String] = Map()
  protected var reportFeatures: Map[String, String] = Map()
  protected var reportEnabled: Boolean = true

  def includeInReport = reportEnabled
  def enableReport() { reportEnabled = true }
  def disableReport() { reportEnabled = false }

  def setRunInfo(info: JobRunInfo) {
    //logger.info("info " + info)
    val runtimeFeatures = Map(
      "iteration" -> 1,
      "analysisName" -> getReportGroup,
      "jobName" -> QJobReport.workAroundSameJobNames(this),
      "intermediate" -> self.isIntermediate,
      "exechosts" -> info.getExecHosts,
      "startTime" -> info.getStartTime,
      "doneTime" -> info.getDoneTime,
      "formattedStartTime" -> info.getFormattedStartTime,
      "formattedDoneTime" -> info.getFormattedDoneTime,
      "runtime" -> info.getRuntimeInMs).mapValues((x:Any) => if (x != null) x.toString else "null")
    reportFeatures = runtimeFeatures ++ userReportFeatures
    // note -- by adding reportFeatures second we override iteration
    // (or any other binding) with the user provided value
  }

  /** The report Group is the analysis name transform to only contain valid GATKReportTable characters */
  def getReportGroup = self.analysisName.replaceAll(GATKReportTable.INVALID_TABLE_NAME_REGEX, "_")

  def getReportFeatureNames: Seq[String] = reportFeatures.keys.toSeq
  def getReportFeature(key: String): String = {
    reportFeatures.get(key) match {
      case Some(x) => x
      case None =>
        logger.warn("getReportFeature called with key %s but no value was found for group %s.  This can be caused by adding user-defined job features to a job with a generic name used elsewhere in the Queue script.  To fix the problem make sure that each group of commands with user-specific features has a unique analysisName".format(key, reportGroup))
        "NA"
    }
  }

  def getReportName: String = getReportFeature("jobName")

  def configureJobReport(features: Map[String, Any]) {
    this.userReportFeatures = features.mapValues(_.toString)
  }

  def addJobReportBinding(key: String, value: Any) {
    this.userReportFeatures += (key -> value.toString)
  }

  // copy the QJobReport information -- todo : what's the best way to do this?
  override def copySettingsTo(function: QFunction) {
    self.copySettingsTo(function)
    function.userReportFeatures = this.userReportFeatures
    function.reportFeatures = this.reportFeatures
  }
}

object QJobReport {
  // todo -- fixme to have a unique name for Scatter/gather jobs as well
  var seenCounter = 1
  var seenNames = Set[String]()

  def workAroundSameJobNames(func: QFunction):String = {
    if ( seenNames.apply(func.jobName) ) {
      seenCounter += 1
      "%s_%d".format(func.jobName, seenCounter)
    } else {
      seenNames += func.jobName
      func.jobName
    }
  }
}
