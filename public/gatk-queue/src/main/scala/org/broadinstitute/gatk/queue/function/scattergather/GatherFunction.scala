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

package org.broadinstitute.gatk.queue.function.scattergather

import java.io.File
import org.broadinstitute.gatk.utils.commandline.{Input, Output}
import org.broadinstitute.gatk.queue.function.QFunction
import org.broadinstitute.gatk.queue.QException
import org.broadinstitute.gatk.utils.io.IOUtils
import collection.JavaConversions._

/**
 * Base class for Gather command line functions.
 */
trait GatherFunction extends QFunction {
  analysisName = "Gather"

  var originalFunction: ScatterGatherableFunction = _

  @Output(doc="The original output of the scattered function")
  var originalOutput: File = _

  @Input(doc="Parts to gather back into the original output")
  var gatherParts: Seq[File] = Nil

  /**
   * Called to initialize the gather function values after all other values have been setup for this function.
   */
  def init() {}

  /**
   * Don't include this @Gather's log file when tracking .done.
   * The done files for original log file being produced will do.
   *
   * The logs from the scatter/gather jobs are concatenated together into the original log.
   * Without removing the logs a .done file would be created for the logs. If a SGF is switched
   * from scatterCount=1 to >1 then this Gather would be "missing" its logs and re-run.
   */
  override protected def statusPaths = outputs

  /**
   * Waits for gather parts to propagate over NFS or throws an exception.
   */
  protected def waitForGatherParts() {
    val missing = IOUtils.waitFor(gatherParts, 120)
    if (!missing.isEmpty)
      throw new QException("Unable to find gather inputs: " + missing.mkString(", "))
  }
}
