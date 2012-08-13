/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.queue.function

/** A mixin that on retry increases the memory limit when certain text is found. */
trait RetryMemoryLimit extends CommandLineFunction {

  /** How to increase the memory. By default doubles the memory. */
  var retryMemoryFunction: (Double => Double) = (2 * _)

  /** Once the threshold is passed, no more memory will be added to memory limit. */
  var memoryLimitThreshold: Option[Double] = None

  /** Various strings to look for to determine we ran out of memory. */
  var memoryLimitErrorText = Seq("OutOfMemory", "you did not provide enough memory", "TERM_MEMLIMIT")

  override def freezeFieldValues() {
    super.freezeFieldValues()
    if (this.memoryLimitThreshold.isEmpty)
      this.memoryLimitThreshold = this.qSettings.memoryLimitThreshold
  }

  override def setupRetry() {
    super.setupRetry()
    if (this.memoryLimitThreshold.isDefined && this.memoryLimit.isDefined) {

      // NOTE: If we're already at or above the memoryLimit, don't do anything.
      if (this.memoryLimit.get < this.memoryLimitThreshold.get) {
          updateMemoryLimits()
      }

    } else {
      updateMemoryLimits()
    }
  }

  def updateMemoryLimits() {
    if (isMemoryError) {
      this.memoryLimit = this.memoryLimit.map(this.retryMemoryFunction)
      this.residentRequest = this.residentRequest.map(this.retryMemoryFunction)
      this.residentLimit = this.residentLimit.map(this.retryMemoryFunction)

      // Rebuffer the memory limit if the limit was set exactly to the request
      if (this.residentLimit == this.residentRequest)
        this.residentLimit = this.residentRequest.map(this.residentLimitBuffer)

      this match {
        case java: JavaCommandLineFunction =>
          java.javaMemoryLimit = java.javaMemoryLimit.map(this.retryMemoryFunction)
        case _ => /* ignore */
      }
    }
  }

  def isMemoryError = this.jobErrorLines.exists(line => this.memoryLimitErrorText.exists(error => line.contains(error)))
}
