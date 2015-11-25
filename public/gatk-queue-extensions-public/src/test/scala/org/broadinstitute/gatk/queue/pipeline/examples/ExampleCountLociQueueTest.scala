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

package org.broadinstitute.gatk.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.gatk.queue.pipeline.{QueueTest, QueueTestSpec}
import org.broadinstitute.gatk.utils.BaseTest

class ExampleCountLociQueueTest {
  @Test(timeOut=36000000)
  def testCountLoci() {
    val testOut = "count.out"
    val spec = new QueueTestSpec
    spec.name = "countloci"
    spec.args = Array(
      " -S " + QueueTest.publicQScriptsPackageDir + "examples/ExampleCountLoci.scala",
      " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
      " -I " + BaseTest.publicTestDir + "exampleBAM.bam",
      " -o " + testOut).mkString
    spec.fileMD5s += testOut -> "ade93df31a6150321c1067e749cae9be"
    QueueTest.executeTest(spec)
  }
}