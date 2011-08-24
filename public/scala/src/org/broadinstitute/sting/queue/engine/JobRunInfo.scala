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

/**
 * Base class containing all of the information about a job run.
 */
class JobRunInfo {
  var startTime: Date = _
  var doneTime: Date = _
  var memUsedInGb: Int = -1
  var status: RunnerStatus.Value = RunnerStatus.DONE

  def getStatus = status

  def getStartTime = startTime
  def getDoneTime = doneTime
  def getFormattedStartTime = formatTime(getStartTime)
  def getFormattedDoneTime = formatTime(getDoneTime)

  val formatter = new SimpleDateFormat("dd.MM.yy/H:mm:ss:SSS");
  private def formatTime(d: Date) = if ( d != null ) formatter.format(d) else "null"

  def getMemoryUsedInGb = memUsedInGb
  def isFilledIn = startTime != null

  def getRuntimeInMs: Long = {
    if ( getDoneTime != null && getStartTime != null )
      getDoneTime.getTime - getStartTime.getTime
    else
      -1
  }

  override def toString: String =
    "started %s ended %s runtime %s using %d Gb memory".format(getFormattedStartTime, getFormattedDoneTime, getRuntimeInMs, getMemoryUsedInGb)
}

object JobRunInfo {
  def default: JobRunInfo = new JobRunInfo()
}