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

package org.broadinstitute.gatk.queue.qscripts.examples

import org.broadinstitute.gatk.queue.function.RetryMemoryLimit
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._

class ExampleRetryMemoryLimit extends QScript {
  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _

  @Input(doc="Bam file to print.", shortName="I")
  var bamFile: File = _

  def script() {
    for (scatterCount <- 1 to 2) {
      val printReads = new PrintReads with RetryMemoryLimit
      // First run with 1m
      printReads.memoryLimit = .001
      // On retry run with 1g
      printReads.retryMemoryFunction = (d => d * 1000)
      printReads.reference_sequence = referenceFile
      printReads.input_file = Seq(bamFile)
      printReads.out = swapExt(bamFile, ".bam", ".out.scattered_%d.bam".format(scatterCount))
      printReads.scatterCount = scatterCount
      add(printReads)
    }
  }
}
