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

import java.util.Date
import java.text.SimpleDateFormat

/**
 * Class containing tracked information about a job run.
 */
  // todo -- it might be nice to have the hostname
class JobRunInfo {
  /** constant date format */
  val formatter = new SimpleDateFormat("yy-MM-dd H:mm:ss:SSS");

  /** The start time with millisecond resolution of this job */
  var startTime: Date = _
  /** The done time with millisecond resolution of this job */
  var doneTime: Date = _
  var exechosts: String = "localhost"

  def getStartTime = startTime
  def getDoneTime = doneTime
  def getFormattedStartTime = formatTime(getStartTime)
  def getFormattedDoneTime = formatTime(getDoneTime)

  /** Helper function that pretty prints the date */
  private def formatTime(d: Date) = if ( d != null ) formatter.format(d) else "null"

  def getExecHosts = exechosts

  /**
   * Was any information set for this jobInfo?  JobInfo can be unset because
   * the job never ran or because it already completed.
   */
  def isFilledIn = startTime != null

  /**
   * How long did the job run (in wall time)?  Returns -1 if this jobInfo isn't filled in
   */
  def getRuntimeInMs: Long = {
    if ( isFilledIn )
      getDoneTime.getTime - getStartTime.getTime
    else
      -1
  }

  override def toString: String =
    "started %s ended %s runtime %s".format(getFormattedStartTime, getFormattedDoneTime, getRuntimeInMs)
}

object JobRunInfo {
  def default: JobRunInfo = new JobRunInfo()
}