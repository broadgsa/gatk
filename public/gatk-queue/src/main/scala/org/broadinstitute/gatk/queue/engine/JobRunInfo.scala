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

import java.text.SimpleDateFormat

/**
 * Class containing tracked information about a job run.
 */
  // todo -- it might be nice to have the hostname
class JobRunInfo {
  /** constant date format */
  val formatter = new SimpleDateFormat("yy-MM-dd H:mm:ss:SSS");

  /** The start time with millisecond resolution of this job */
  var startTime: java.util.Date   = _
  /** The done time with millisecond resolution of this job */
  var doneTime: java.util.Date = _
  var exechosts: String = "localhost"

  def getStartTime: String = getTime(startTime)
  def getDoneTime: String = getTime(doneTime)
  def getFormattedStartTime = formatTime(startTime)
  def getFormattedDoneTime = formatTime(doneTime)

  /** Helper function that returns the time of the date */
  private def getTime(d: java.util.Date): String = if ( d != null ) d.getTime.toString else "null"

  /** Helper function that pretty prints the date */
  private def formatTime(d: java.util.Date): String = if ( d != null ) formatter.format(d) else "null"

  def getExecHosts = exechosts

  /**
   * Was any information set for this jobInfo?  JobInfo can be unset because
   * the job never ran or because it already completed.
   */
  def isFilledIn = startTime != null && doneTime != null

  /**
   * How long did the job run (in wall time)?  Returns -1 if this jobInfo isn't filled in
   */
  def getRuntimeInMs: Long = {
    if ( isFilledIn )
      doneTime.getTime - startTime.getTime
    else
      -1
  }

  override def toString: String =
    "started %s ended %s runtime %s".format(getFormattedStartTime, getFormattedDoneTime, getRuntimeInMs)
}

object JobRunInfo {
  def default: JobRunInfo = new JobRunInfo()
}