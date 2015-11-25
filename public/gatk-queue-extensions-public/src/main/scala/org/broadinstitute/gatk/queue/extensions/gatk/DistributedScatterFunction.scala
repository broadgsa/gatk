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

package org.broadinstitute.gatk.queue.extensions.gatk

import java.io.File
import org.broadinstitute.gatk.queue.function.InProcessFunction
import org.broadinstitute.gatk.queue.function.scattergather.{ScatterFunction, CloneFunction}
import org.broadinstitute.gatk.utils.commandline.Output
import util.Random

/**
 * An scatter function that uses the Distributed GATK.
 */
class DistributedScatterFunction extends ScatterFunction with InProcessFunction {
  @Output(doc="processingTracker")
  var processingTracker: File = _

  override def init() {
    this.processingTracker = new File(this.commandDirectory, "processingTracker.%8d".format(Random.nextInt(100000000)))
  }

  override def initCloneInputs(cloneFunction: CloneFunction, index: Int) {
    cloneFunction.setFieldValue("processingTracker", this.processingTracker)
  }

  def run() {
    /* doesn't actually need to run. */
  }
}
